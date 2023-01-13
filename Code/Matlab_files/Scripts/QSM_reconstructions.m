clearvars; clc; close all
cd ../
addpath(genpath('functions/'))
principalFolder = pwd;
f = filesep;
addpath(genpath('../mritools_Linux_3.4.3/'))
diary Simulate_data.txt

% cd ../../Imagenes/Phantom/
% IMAGES_FOLDER = [pwd, f];
IMAGES_FOLDER = fullfile('..', '..', 'Imagenes', 'Phantom_real_data', 'Orientations');
MODE_FOLDER = {'Susceptibility_sti', 'Diffusion_sti'};
OUT_FOLDER = 'QSM_reconstructions';
MAGN_FOLDER = fullfile(IMAGES_FOLDER, '..', 'Maps');
MAGN_FILE = fullfile(MAGN_FOLDER, 'm0.nii.gz');

B0 = 3; %T
gamma = 42.14; % MHz/T
phase_scale = 2*pi*gamma*B0;

PHI_FILE = 'phi_15_orientations.nii.gz';
N_ORIENTATIONS = 15;

% FANSI options
options = [];
options.tgv = false;
options.nonlinear = true;
options.iterations = 1000;
options.update = 0.1;
options.gradient_mode = 1;
B0_dir = [0, 0, 1];

% Regularization parameters
alpha = 1e-6;
mu = 100*alpha;


%% Define reconstruction parameters

nii_magn = load_untouch_nii(MAGN_FILE);
magn = single(nii_magn.img);
magn_use = (magn - min(magn(:)))/(max(magn(:) - min(magn(:)))); % Normalize magn between [0, 1]
disp(['Magnitude image size: ', num2str(size(magn_use))])

voxel_size = nii_magn.hdr.dime.pixdim(2:4);

for n_tensor = 1:length(PHI_FILE)
    disp('=================================')
    disp(['Simulating ', MODE_FOLDER{n_tensor}, ' mode....'])
    sim_folder = fullfile(IMAGES_FOLDER, MODE_FOLDER{n_tensor});
    phi_file = fullfile(sim_folder, PHI_FILE);
    nii_phi = load_untouch_nii(phi_file);
    n_orientation = 1;
    %     for n_orientation = 1:N_ORIENTATIONS
    qsm_name = ['QSM_', num2str(alpha), '.nii.gz'];
    qsm_file = fullfile(sim_folder, qsm_name);
    
    phase_use = single(nii_phi.img(:, :, :, n_orientation) / phase_scale);    
    disp(['Phase image size: ', num2str(size(phase_use))])
    
    % QSM Reconstruction
    outf = FANSI(phase_use, magn_use, voxel_size, alpha, mu, 0.0054, options, B0_dir);
    nii_magn.img = outf.x;
    save_untouch_nii(nii_magn, qsm_file);
    disp([qsm_file, ' saved'])
%     end
end
diary off