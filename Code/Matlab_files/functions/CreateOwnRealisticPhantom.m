function CreateOwnRealisticPhantom(params)

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
% label_names = 

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

% M0=load_nii(params.M0map_file);

%% we also load some originally computed field maps so that we know where we trust or we don't trust the R2* values
HighGradMask=load_nii(params.highGradMask_file);
RealFreqMap=load_nii(params.rawField_file);

disp('Applying 3D gradient field')
[fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(RealFreqMap.img,[],0);
fieldgradient=(sos(fieldgradient,4));
fieldgradient(FinalSegm>=11)=0;

fieldgradient=fieldgradient.*single(HighGradMask.img);
%     %define acceptable threshold ???
% thresholdgrad=prctile(fieldgradient(fieldgradient~=0),[10 90]); % is this a good threshold, check on the fv of the fiedlgradient
% looking at the figure 0.05 seems to be the gradient close to the edges
% looking at the figure 0.2 seems to be a too large gradient

disp('Truncate R2* map')
fieldgradient(and(fieldgradient~=0,fieldgradient>0.2))=0.2;
fieldgradient(and(fieldgradient~=0,fieldgradient<0.05))=0.05;
fieldgradient(fieldgradient~=0)=(fieldgradient(fieldgradient~=0)-0.05)/(0.2-0.05);
WeightingBetweenModels=cos(fieldgradient);

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
        %    indexes=find(and(FinalSegm==k,HighGradMask.img==0));
        %    indexes_HG=find(and(FinalSegm==k,HighGradMask.img==1));
        
        if sum(Mask(:))>1
            
            Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
            % ensures the modulation is only applied in regions of low gradient
            if k <= 10
                ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;
                chitemp=...
                    R2star2chi(k)*(R2star.img(LG_mask_index)-mean(R2star.img(indexes)))+...
                    R12chi(k)*(R1.img(LG_mask_index)-mean(R1.img(indexes)));
                
                % trying to find a way to avoid creating outlyers because of errors on the R2* and R1 maps...
                modulation_std = std(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
                modulation_mean = mean(chitemp(Mask_smooth(LG_mask_index)>smoothThresh));
                % makes all the outliers have the average value of the segmented tissue
                chitemp(and(Mask_smooth(LG_mask_index)>smoothThresh2,abs(chitemp-modulation_mean)>NSTD*modulation_std))=modulation_mean;
                
                
                Chimap3(LG_mask_index)=Chimap3(LG_mask_index)+Mask_smooth(LG_mask_index)...
                    .*(chiref(k)+1*chitemp);
                
                % does the all R1 based method
                
                chitemp=...
                    R1only2chi(k)*(R1.img(LG_mask_index)-mean(R1.img(indexes)));
                
                % trying to find a way to avoid creating outlyers because of errors on the R2* and R1 maps...
                Chimap4(LG_mask_index)=Chimap4(LG_mask_index)+Mask_smooth(LG_mask_index)...
                    .*(chiref(k)+chitemp);
                
            else
                %    this is not quite correct because there will be some suseptibilitz
                %    missing in the regions arround vessels, dura and so on... hopefullz it
                %    is a small effect
                chitemp=R2star2chi(k)*(R2star.img(indexes)-mean(R2star.img(indexes)))+...
                    R12chi(k)*(R1.img(indexes)-mean(R1.img(indexes)));
                
                modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
                modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
                chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
                chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
                
                
                Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                    .*(chiref(k)+1*chitemp);
                
                % only based on R1
                chitemp=...
                    R1only2chi(k)*(R1.img(indexes)-mean(R1.img(indexes)));
                
                modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
                modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
                chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
                chitemp(and(Mask(indexes)>smoothThresh2,chitemp<modulation_mean+NSTD*modulation_std))=modulation_mean;
                Chimap4(indexes)=Chimap4(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
                    .*(chiref(k)+chitemp);
                
                
            end
        end
    end
    
    Chimap5=WeightingBetweenModels.^4.*Chimap3+(1-WeightingBetweenModels.^4).*Chimap4;    
    clear chitemp
    chi_lambda(:, :, :, nn) = Chimap5;
    disp(['lambda ', num2str(nn)])
    
end
nii_chi = make_nii(chi_lambda, voxel_size);
save_nii(nii_chi,params.OutputLambda)
disp([params.OutputLambda, ' saved'])


