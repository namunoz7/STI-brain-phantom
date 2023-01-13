function nii = generate_eigenvalues_sti(params)
% Function that generates the eigenvalues of the susceptibility tensor,
% according to the parameters, i.e., mean values of each region in the RC1
% challenge

NSTD=3;
smoothThresh=0.9; % region where the width of accepted values is calculated
smoothThresh2=0.05; % region where the width of accepted values is enforced
FWHM_V=1.2;

labels = readtable(params.ChiModulation_file);
% labels = struct2table(label);
FinalSegment = load_untouch_nii(params.Segmentation_file); % this has the infromation
FinalSegm = FinalSegment.img;

FAonlyChi = labels.FAonlyChi;

nii_FA = load_untouch_nii(params.FAmap_file);
N = size(FinalSegm);
% voxel_size = nii_FA.hdr.dime.pixdim(2:4);
% HighGradMask =load_nii(params.highGradMask_file);
% RealFreqMap =load_nii(params.rawField_file);

% WeightingBetweenModels = apply_gradient(HighGradMask, RealFreqMap, FinalSegm);
% LG_mask_index=find(or(HighGradMask.img==0,HighGradMask.img==1));
chi_lambda = zeros([N, 3]);


disp('Simulating eigenvalues')
for nn = 1:3
    Chimap3 = zeros(N);
    ProbabilityAccumulated = zeros(N);
    if (nn == 1)
        disp('Lambda 1')
        chiref = labels.chiref1;
    elseif (nn == 2)
        disp('Lambda 2')
        chiref = labels.chiref2;
    else
        disp('Lambda 3')
        chiref = labels.chiref3;
    end
    for n_label = 1:size(labels,1)
        disp(['... Simulating ', labels.name{n_label}])
%         tmp(FinalSegm == n_label) = chiref(n_label);
        Mask=FinalSegm==n_label;
        indexes=find(FinalSegm==n_label);
        Mask_smooth=real(smooth3D(double(Mask),FWHM_V,[1,1,1]));
%         Chimap3 = add_fa_tissue(Mask_smooth, LG_mask_index, indexes, FA.img, FAonlyChi(n_label), chiref(n_label), Chimap3, smoothThresh, smoothThresh2, NSTD);
        Chimap3 = add_fa_tissue(Mask_smooth, indexes, nii_FA.img, FAonlyChi(n_label), chiref(n_label), Chimap3, smoothThresh, smoothThresh2, NSTD);
        ProbabilityAccumulated=ProbabilityAccumulated+Mask_smooth;        
    end
    
%     disp(['... Simulating ', labels.name{11}])
%     Mask=FinalSegm==11;
%     indexes=find(FinalSegm==11);    
%     tissue=FAonlyChi(11)*(nii_FA.img(indexes)-mean(nii_FA.img(indexes)));    
%     % trying to find a way to avoid creating outlyers because of errors on the R2* and R1 maps...
%     modulation_std = std(tissue(Mask(indexes)>smoothThresh));
%     modulation_mean = mean(tissue(Mask(indexes)>smoothThresh));
%     % makes all the outliers have the average value of the segmented tissue
%     tissue(and(Mask(indexes)>smoothThresh2,abs(tissue-modulation_mean)>NSTD*modulation_std))=modulation_mean;
%     Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes) - ProbabilityAccumulated(indexes)).*(chiref(11)+tissue);
%     
%     disp(['... Simulating ', labels.name{16}])
%     Mask=FinalSegm==16;
%     indexes=find(FinalSegm==16);
%     tissue=FAonlyChi(16)*(nii_FA.img(indexes)-mean(nii_FA.img(indexes)));    
%     % trying to find a way to avoid creating outlyers because of errors on the R2* and R1 maps...
%     modulation_std = std(tissue(Mask(indexes)>smoothThresh));
%     modulation_mean = mean(tissue(Mask(indexes)>smoothThresh));
%     % makes all the outliers have the average value of the segmented tissue
%     tissue(and(Mask(indexes)>smoothThresh2,abs(tissue-modulation_mean)>NSTD*modulation_std))=modulation_mean;
%     Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes) - ProbabilityAccumulated(indexes)).*(chiref(16)+tissue);
    
%     Chimap5=WeightingBetweenModels.^4.*Chimap3;
%     chi_lambda(:, :, :, nn) = Chimap5;
    chi_lambda(:, :, :, nn) = Chimap3;
end
nii_FA.img = chi_lambda;
nii_FA.hdr.dime.dim(5) = 3;
nii_FA.hdr.dime.dim(1) = 4;
save_untouch_nii(nii_FA,params.OutputLambda)
disp([params.OutputLambda, ' saved'])

nii = nii_FA;
end