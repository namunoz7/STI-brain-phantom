clearvars; clc; close all
cd ../
addpath(genpath('functions/'))
principalFolder = pwd;
f = filesep;
% addpath(genpath('../mritools_Linux_3.4.3/'))
% addpath(genpath('../STISuite_V3.0/'))
diary Simulate_data.txt

% IMAGES_FOLDER = fullfile('..', '..', '..', 'Imagenes', 'Phantom_real_data', 'Noised_data');
IMAGES_FOLDER = fullfile('..', '..', 'Imagenes', 'Phantom_real_data', 'Noised_data');
RECONSTRUCTION_FOLDERS = {'STI_suite', 'COSMOS_STI'};
MODEL_FOLDER = {'Susceptibility_sti', 'Diffusion_sti'};
ORIENTATIONS_FOLDER = {'6_orientations', '12_orientations'};
ANGLES_ORIENTATIONS = {'angles_1.mat', 'angles_2.mat'};
PHASE_NAMES = {'phi_6_orientations.nii.gz', 'phi_12_orientations.nii.gz'};
CHI_RECONSTRUCTIONS = {'chi_6_orientations.nii.gz', 'chi_12_orientations.nii.gz'};

GT_FOLDER = fullfile(IMAGES_FOLDER, '..', 'Phantom_tensor');
TENSOR_NAMES = {'chi_sti_filt.nii.gz', 'chi_dti_filt.nii.gz'};

MASK_FILE = fullfile(IMAGES_FOLDER, '..', 'Masks', 'brain_mask.nii.gz');
nii_mask = load_untouch_nii(MASK_FILE);
mask = logical(repmat(nii_mask.img, [1, 1, 1, 6]));

% STI suite params
parallel_flag = 'off';
B0 = 3;  % T
TE = 12e-3;  % ms
gamma = 42.58;
phase_scale = 2*pi*gamma*B0;
std_noise = 0;

tmp_file = fullfile(GT_FOLDER, 'chi_sti_filt.nii.gz');
nii_tmp = load_untouch_nii(tmp_file);
error_reconstructions = zeros(2, 2, 2);
error_file = fullfile(IMAGES_FOLDER, 'error_reconstructions.mat');

for n_reconstruction = 1:2
    disp('-------------------------')
    disp(RECONSTRUCTION_FOLDERS{n_reconstruction})
    actual_reconstruction = fullfile(IMAGES_FOLDER, RECONSTRUCTION_FOLDERS{n_reconstruction});
    if ~isfolder(actual_reconstruction)
        mkdir(actual_reconstruction)
        disp([actual_reconstruction, ' created'])
    end
    for n_model = 1:2
        disp(['... ', MODEL_FOLDER{n_model}])
        actual_model = fullfile(actual_reconstruction, MODEL_FOLDER{n_model});
        gt_model_folder = fullfile(IMAGES_FOLDER, 'STI_pytorch', MODEL_FOLDER{n_model});
        if ~isfolder(actual_model)
            mkdir(actual_model)
            disp([actual_model, ' created'])
        end
        gt_chi_file = fullfile(GT_FOLDER, TENSOR_NAMES{n_model});
        nii_chi = load_untouch_nii(gt_chi_file);
        gt_chi = nii_chi.img;
        for n_orientation = 1:2
            disp(['... ... ', ORIENTATIONS_FOLDER{n_orientation}])
            actual_orientation = fullfile(actual_model, ORIENTATIONS_FOLDER{n_orientation});
            if ~isfolder(actual_orientation)
                mkdir(actual_orientation)
                disp([actual_orientation, ' created'])
            end
            gt_orientation_folder = fullfile(gt_model_folder, ORIENTATIONS_FOLDER{n_orientation});
            actual_phase = fullfile(gt_orientation_folder, PHASE_NAMES{n_orientation});
            direction_file = fullfile(gt_orientation_folder, ANGLES_ORIENTATIONS{n_orientation});
            chi_rec_file = fullfile(actual_orientation, CHI_RECONSTRUCTIONS{n_orientation});

            nii_phi = load_untouch_nii(actual_phase);
            nii_phi.img = nii_phi.img + std_noise*randn(size(nii_phi.img));
            H_matrix = load(direction_file);
            direction_field = get_direction_components(H_matrix);
            disp(direction_field)
            if n_reconstruction == 1  % STI_suite
                phase_use = nii_phi.img;
                chi_rec_2 = STI_Parfor(phase_use, direction_field, 1, 1, parallel_flag);
                chi_rec = real(chi_rec_2(:, :, 49:176, :))/B0;
            else  % Cosmos_STI
                phase_use = nii_phi.img;
                chi_rec = sti_lsqr(phase_use, direction_field, mask(:, :, :, 1));
            end
            nii_tmp.img = chi_rec;
            save_untouch_nii(nii_tmp, chi_rec_file)
            disp([chi_rec_file, ' saved.'])
        end
    end
end
save(error_file, "error_reconstructions")

function direction_field = get_direction_components(H_matrix)
vec_psi = H_matrix.vec_psi;
vec_theta = H_matrix.vec_theta;
n_orientations = length(vec_psi);
direction_field = zeros(n_orientations, 3);
direction_field(:, 1) = sin(vec_theta) .* cos(vec_psi);  % Hx
direction_field(:,2) = sin(vec_theta) .* sin(vec_psi);  % Hy
direction_field(:, 3) = cos(vec_theta);  % Hz
end

