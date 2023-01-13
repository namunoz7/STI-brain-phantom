clearvars; clc
cd ..
addpath(genpath('functions'))
diary generate_phantom.txt
% cd ../
principal_folder = pwd;

IMAGES_FOLDER = fullfile('..', '..', 'Imagenes', 'Phantom_real_data');
OUT_FOLDER = fullfile(IMAGES_FOLDER, 'Phantom_tensor');
DTI_FOLDER = fullfile(IMAGES_FOLDER, 'Diffusion_data', 'DTI_reg_files');
STI_FOLDER = fullfile(IMAGES_FOLDER, 'Susceptibility_data');

%% Generate eigenvalues

ModelParams.FAmap_file = fullfile(DTI_FOLDER, 'dti_reg_FA.nii.gz');
ModelParams.Segmentation_file = fullfile(IMAGES_FOLDER, 'Masks', 'atlas_mask.nii.gz');
% ModelParams.rawField_file = fullfile(IMAGES_FOLDER, 'raw', 'MP2RAGEME_total-field_ppm_head_r2s-header.nii.gz');
ModelParams.BrainMask_file = fullfile(IMAGES_FOLDER, 'Masks', 'brain_mask.nii.gz');

% file that defines regions where R2* values can and can't be trusted
% ModelParams.highGradMask_file = fullfile(IMAGES_FOLDER, 'masks', 'highgrad.nii.gz');

% File with the various parameters to create a susceptbility map
% ModelParams.ChiModulation_file = fullfile(IMAGES_FOLDER, 'data', 'chimodel', 'parametersTensor.mat');
% ModelParams.ChiModulation_file = fullfile(IMAGES_FOLDER, 'parameters_tensor.xlsx');
ModelParams.ChiModulation_file = fullfile(IMAGES_FOLDER, 'parameters_tensor_exp_data.xlsx');

% ModelParams.ChiModulation_file = 'data/chimodel/paramatersNoCalcModel.mat';
ModelParams.OutputLambda = fullfile(OUT_FOLDER, 'eig_phantom_sti.nii.gz');
% ModelParams.OutputChiModel_file = {fullfile(OUT_FOLDER, 'lambda3_sorted.nii.gz'), ...
%     fullfile(OUT_FOLDER, 'lambda2_test_sorted.nii.gz'), ...
%     fullfile(OUT_FOLDER, 'lambda1_test_sorted.nii.gz')};


disp('Generating eigenvalues')
ModelParams.eigenValues = generate_eigenvalues_sti(ModelParams);
% chi_lambda = create_eigen_values(ModelParams);
% CreateOwnRealisticPhantom(ModelParams)

%% Multiply eigenvectors with eigenvalues to generate the tensor

% ModelParams.OutputTensor = [IMAGES_FOLDER, '\fa_tissued_phantom\eig.nii.gz'];
% ModelParams.R1map_file = fullfile(IMAGES_FOLDER, 'maps', 'R1_segmented.nii.gz');
% ModelParams.R2starmap_file = fullfile(IMAGES_FOLDER, 'maps', 'R2s_segmented.nii.gz');
% ModelParams.M0map_file = fullfile(IMAGES_FOLDER, 'maps', 'M0_segmented.nii.gz');

for n = 1:2
    if n == 1  % Tensor with DTI files
        disp('-------------------')
        disp('Simulating tensor with DTI eigenvectors')
        ModelParams.V1_eigenvector = fullfile(DTI_FOLDER, 'dti_reg_V1.nii.gz');
        ModelParams.V2_eigenvector = fullfile(DTI_FOLDER, 'dti_reg_V2.nii.gz');
        ModelParams.V3_eigenvector = fullfile(DTI_FOLDER, 'dti_reg_V3.nii.gz');
        
        ModelParams.Out_V1_filtered = fullfile(DTI_FOLDER, 'V1_filtered.nii.gz');
        ModelParams.Out_V2_filtered = fullfile(DTI_FOLDER, 'V2_filtered.nii.gz');
        ModelParams.Out_V3_filtered = fullfile(DTI_FOLDER, 'V3_filtered.nii.gz');

        ModelParams.OutChi = fullfile(OUT_FOLDER, 'chi_dti.nii.gz');
        ModelParams.Out_filt_chi = fullfile(OUT_FOLDER, 'chi_dti_filt.nii.gz');

    else
        disp('-------------------')
        disp('Simulating tensor with STI eigenvectors')
        ModelParams.V1_eigenvector = fullfile(STI_FOLDER, 'eig_v3_sti.nii.gz');
        ModelParams.V2_eigenvector = fullfile(STI_FOLDER, 'eig_v2_sti.nii.gz');
        ModelParams.V3_eigenvector = fullfile(STI_FOLDER, 'eig_v1_sti.nii.gz');
        
        ModelParams.Out_V1_filtered = fullfile(STI_FOLDER, 'V1_filtered.nii.gz');
        ModelParams.Out_V2_filtered = fullfile(STI_FOLDER, 'V2_filtered.nii.gz');
        ModelParams.Out_V3_filtered = fullfile(STI_FOLDER, 'V3_filtered.nii.gz');

        ModelParams.OutChi = fullfile(OUT_FOLDER, 'chi_sti.nii.gz');
        ModelParams.Out_filt_chi = fullfile(OUT_FOLDER, 'chi_sti_filt.nii.gz');
        
    end

    ModelParams.chi = generate_tensor(ModelParams);
    chi2 = zeros(size(ModelParams.chi.img));
    FWHM1 = 1.7; sigma = FWHM1/2.3548;

    for nn = 1:6
        chi2(:, :, :, nn) = imgaussfilt3(ModelParams.chi.img(:, :, :, nn), sigma);
    end

    ModelParams.chi.img = chi2;
    save_untouch_nii(ModelParams.chi, ModelParams.Out_filt_chi)
    disp([ModelParams.Out_filt_chi, ' saved'])
end

diary off



