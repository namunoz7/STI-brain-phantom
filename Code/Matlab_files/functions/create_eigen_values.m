function chi_lambda = create_eigen_values(params)
% params.R1map_file
% params.R2starmap_file
% params.M0map_file
% params.Segmentation_file
% params.ChiModulation_file
% params.Segmentation_file
% params.rawField_file
% params.highGradMask_file
% params.OutputChiModel_file

NSTD=3;
smoothThresh=0.9; % region where the width of accepted values is calculated
smoothThresh2=0.05; % region where the width of accepted values is enforced

load(params.ChiModulation_file); % this has the information on the reference susceptibility in each region
labels = struct2table(label);
R12chi = 0.5*labels.R12chi;
R1only2chi = 0.5*labels.R1only2chi;
R2star2chi = 0.5*labels.R2star2chi;

for k=1:length(label)
    R12chi(k)=label(k).R12chi;
    R1only2chi(k)=label(k).R1only2chi;
    R2star2chi(k)=label(k).R2star2chi;
end
FinalSegment = load_nii(params.Segmentation_file); % this has the infromation
FinalSegm = FinalSegment.img;
disp('Wavelet filtering')
disp('Filtering R1 map')
R1 = load_nii(params.R1map_file);
tmp=Wavedec3Denoising(R1.img,25/1000,8,'db2','verysoft',double(R1.img~=0));
R1.img = abs(real(tmp)); % same as for segmentation purposes
voxel_size = R1.hdr.dime.pixdim(2:4);

disp('Filtering R2* map')
R2star= load_nii(params.R2starmap_file);
tmp=Wavedec3Denoising(R2star.img,45/10,8,'db2','soft',double(R1.img~=0)); % same as for segmentation purposes
R2star.img = abs(real(tmp ));

%% we also load some originally computed field maps so that we know where we trust or we don't trust the R2* values

% M0=load_nii(params.M0map_file);

HighGradMask =load_nii(params.highGradMask_file);
RealFreqMap =load_nii(params.rawField_file);

WeightingBetweenModels = apply_gradient(HighGradMask, RealFreqMap, FinalSegm);

FWHM_V=1.2;
LG_mask_index=find(or(HighGradMask.img==0,HighGradMask.img==1));


ProbabilityAccumulated=0*zeros(size(FinalSegm));

disp('Start creating tensor eigen values')
chi_lambda = zeros([size(FinalSegm), 3]);
for nn = 1:3
    if (nn == 1)
        chiref = labels.chiref1;
    elseif (nn == 2)
        chiref = labels.chiref2;
    else
        chiref = labels.chiref3;
    end
    Chimap3=zeros(size(FinalSegm));% in this one the modulation on the regions of high gradients is only based on the R1 maps
    Chimap4=zeros(size(FinalSegm));% in this one the modulation is only based on the R1 maps
    for k=1:length(label)
        
        Mask=FinalSegm==k;
        indexes=find(FinalSegm==k);
        
        if sum(Mask(:))>1
            
            Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
            if k <= 10
                ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;
                [Chimap3, Chimap4] = simulate_first_tissues(Mask_smooth, LG_mask_index, indexes, R2star.img, R1.img, R2star2chi(k), R1only2chi(k),...
                    R12chi(k), chiref(k), Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD);
            else
                [Chimap3, Chimap4] = simulate_second_tissues(Mask, indexes, R2star.img, R1.img, R2star2chi(k), R1only2chi(k),...
                    R12chi(k), chiref(k), Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD, ProbabilityAccumulated);
            end
        end
    end
    
    Chimap5=WeightingBetweenModels.^4.*Chimap3+(1-WeightingBetweenModels.^4).*Chimap4;    
    chi_lambda(:, :, :, nn) = Chimap5;
    disp(['lambda ', num2str(nn)])
    
end
% nii_chi = make_nii(chi_lambda, voxel_size);
% save_nii(nii_chi,params.OutputLambda)
% disp([params.OutputLambda, ' saved'])
end

