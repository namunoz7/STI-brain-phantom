clearvars; close all; clc
cd ../../
principalFolder = pwd;
IMAGES_FOLDER =  '../../Imagenes/Phantom_real_data/';
GROUND_TRUTH_ROOT = fullfile(IMAGES_FOLDER, 'Phantom_tensor');
SIMULATED_ROOT = fullfile(IMAGES_FOLDER, 'Simulated_data');
NOISED_ROOT = fullfile(IMAGES_FOLDER, 'Noised_data');
MASK_FILE = fullfile(IMAGES_FOLDER, 'Masks', 'brain_mask.nii.gz');
nii_mask_1 = load_untouch_nii(MASK_FILE);
mask = logical(repmat(nii_mask_1.img, [1, 1, 1, 6]));
RECOSNTRUCTION_FOLDERS = {'COSMOS_STI', 'STI_pytorch', 'STI_suite'};

PHANTOM_MODELS_ROOT = {'Diffusion_sti', 'Susceptibility_sti'};
ORIENTATION_ROOT = {'6_orientations', '12_orientations'};
TENSOR_NAMES = {'chi_dti_filt.nii.gz', 'chi_sti_filt.nii.gz'};

RECONSTRUCTED_NAMES = {'chi_6_orientations.nii.gz', 'chi_12_orientations.nii.gz'};

error_reconstructions = zeros(3, 2, 2);
error_file = fullfile(NOISED_ROOT, 'error_reconstructions.mat');

for n_reconstruction = 1:length(RECOSNTRUCTION_FOLDERS)
    actual_reconstruction = fullfile(NOISED_ROOT, RECOSNTRUCTION_FOLDERS{n_reconstruction});
    disp('-----------------------------------')
    disp(RECOSNTRUCTION_FOLDERS{n_reconstruction})
    for n_model = 1:2
        actual_model_root = fullfile(actual_reconstruction, PHANTOM_MODELS_ROOT{n_model});
        actual_gt_file = fullfile(GROUND_TRUTH_ROOT, TENSOR_NAMES{n_model});
        nii_gt_chi = load_untouch_nii(actual_gt_file);
        gt_chi = nii_gt_chi.img;
        disp(['... ', PHANTOM_MODELS_ROOT{n_model}])
        for n_orientation = 1:2
            actual_orientation_file = fullfile(actual_model_root, ORIENTATION_ROOT{n_orientation}, RECONSTRUCTED_NAMES{n_orientation});
            nii_chi = load_untouch_nii(actual_orientation_file);
            chi_img = nii_chi.img;
            disp(['... ', ORIENTATION_ROOT{n_orientation}])
            error_chi = sqrt(mean((gt_chi(mask) - chi_img(mask)).^2))/sqrt(mean(gt_chi(mask).^2));
            disp(error_chi)
            error_reconstructions(n_reconstruction, n_model, n_orientation) = error_chi;
        end
    end
end
save(error_file, "error_reconstructions")