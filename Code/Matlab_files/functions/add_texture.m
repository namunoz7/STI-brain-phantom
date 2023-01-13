function nii_eig = add_texture(params)
NSTD=3;
smoothThresh=0.9; % region where the width of accepted values is calculated
smoothThresh2=0.05; % region where the width of accepted values is enforced
FWHM_V=1.2;

FinalSegment = load_nii(params.Segmentation_file); % this has the infromation
FinalSegm = FinalSegment.img;

load(params.ChiModulation_file); % this has the information on the reference susceptibility in each region
labels = struct2table(label);
R12chi = labels.R12chi;
R1only2chi = labels.R1only2chi;
R2star2chi = labels.R2star2chi;

disp('Wavelet filtering')
disp('Filtering R1 map')
R1 = load_nii(params.R1map_file);

N = size(FinalSegm);

eig = zeros([N, 3]);

voxel_size = R1.hdr.dime.pixdim(2:4);

tmp=Wavedec3Denoising(R1.img,25/1000,8,'db2','verysoft',double(R1.img~=0));
R1 = abs(real(tmp)); % same as for segmentation purposes

disp('Filtering R2* map')
R2star= load_nii(params.R2starmap_file);
tmp=Wavedec3Denoising(R2star.img,45/10,8,'db2','soft',double(R1~=0)); % same as for segmentation purposes
R2star = abs(real(tmp));

%% we also load some originally computed field maps so that we know where we trust or we don't trust the R2* values
HighGradMask=load_nii(params.highGradMask_file);
RealFreqMap=load_nii(params.rawField_file);

disp('Applying 3D gradient field')
[fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(RealFreqMap.img,[],0);
fieldgradient=(sos(fieldgradient,4));
fieldgradient(FinalSegm>=11)=0;

fieldgradient=fieldgradient.*single(HighGradMask.img);

disp('Truncate R2* map')
fieldgradient(and(fieldgradient~=0,fieldgradient>0.2))=0.2;
fieldgradient(and(fieldgradient~=0,fieldgradient<0.05))=0.05;
fieldgradient(fieldgradient~=0)=(fieldgradient(fieldgradient~=0)-0.05)/(0.2-0.05);
WeightingBetweenModels=cos(fieldgradient);

LG_mask_index=find(or(HighGradMask.img==0,HighGradMask.img==1));

disp('Start creating tensor eigen values')

for nn = 1:3
    if (nn == 1)
        chiref = labels.chiref1;
    elseif (nn == 2)
        chiref = labels.chiref2;
    else
        chiref = labels.chiref3;
    end
    Chimap3 = zeros(N);
    Chimap4 = zeros(N);
    ProbabilityAccumulated=zeros(N);
    disp([num2str(nn), ' chi tissue'])
    for n_label = 1:10
        disp(['... ', label(n_label).name])
        Mask=FinalSegm==n_label;
        indexes=find(FinalSegm==n_label);
        if sum(Mask(:))>1
            Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
            ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;
            [Chimap3] = simulate_first_tissues(Mask_smooth, LG_mask_index, indexes, R2star, R1, chiref(n_label), ...
                R2star2chi(n_label), R1only2chi(n_label), R12chi(n_label), Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD);
        end
    end
    disp(['... ', label(11).name])
    Mask=FinalSegm==11;
    indexes=find(FinalSegm==11);
    [Chimap3] = simulate_second_tissues(Mask, indexes, R2star, R1, chiref(11), R2star2chi(11), R1only2chi(11),...
        R12chi(11), Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD, ProbabilityAccumulated);
    
    disp(['... ', label(16).name])
    Mask=FinalSegm==16;
    indexes=find(FinalSegm==16);
    [Chimap3] = simulate_second_tissues(Mask, indexes, R2star, R1, chiref(16), R2star2chi(16), R1only2chi(16),...
        R12chi(16), Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD, ProbabilityAccumulated);
        
%     eig(:, :, :, nn) = (1-WeightingBetweenModels.^4).*Chimap4 + WeightingBetweenModels.^4.*Chimap3;
    eig(:, :, :, nn) = WeightingBetweenModels.^4.*Chimap3 ; % + (1-WeightingBetweenModels.^4).*Chimap4;
    
end
nii_eig = make_nii(eig, voxel_size);
save_nii(nii_eig,params.OutputLambda)
disp([params.OutputLambda, ' saved'])
end