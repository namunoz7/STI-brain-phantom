clearvars; clc; close all
cd ../
addpath(genpath('functions/'))
cd ../
principalFolder = pwd;
f = filesep;
addpath(genpath('../mritools_Linux_3.4.3/'))
diary Simulate_data.txt

% cd ../../Imagenes/Phantom/
% IMAGES_FOLDER = [pwd, f];
IMAGES_FOLDER = fullfile('..', '..', 'Imagenes', 'Phantom_real_data');
SIMULATION_FOLDER = 'Simulated_data';
MODE_FOLDER = {'Susceptibility_sti', 'Diffusion_sti'};
STI_FILE = {'chi_sti_filt.nii.gz', 'chi_dti_filt.nii.gz'};

OUT_FOLDER = {'Phantom_no_microstructure/', 'Phantom_microstructure/'};
B0 = 7;
gyro = 42.58;
TE = 12e-3;

%% Define simulation parameters

SeqParams.n_orietations = 12;
SeqParams.TR = 50e-3;                          % Repetition time in secs
SeqParams.TE = [4, 12, 20, 28] .*1e-3;     % Echo time in seconds
SeqParams.FlipAngle = 15;                      % flip angle in degrees above ernst angle (13) to improve gre contrast between Grey and white matter

SimParams.phase_scale = 2*pi*gyro*TE*B0;
SimParams.B0 = 7; % T
SimParams.PhaseOffset = true;
SimParams.gyro = 42.58; % MHz/T
SimParams.number_orientations = 12;
SimParams.Res = [1 1 1];
SimParams.Shimm = 1;
simulate_microstructure = [false, true];
for n_tensor = 1:length(STI_FILE)
    disp('=================================')
    disp(['Simulating ', MODE_FOLDER{n_tensor}, ' mode....'])
    for n_mu = 1:length(OUT_FOLDER)
        SimParams.microstructure = simulate_microstructure(n_mu);
        SimParams.sim_folder = fullfile(IMAGES_FOLDER, SIMULATION_FOLDER, MODE_FOLDER{n_tensor});

        ModelParams.out_folder = OUT_FOLDER{n_mu};        
        ModelParams.chi_tensor = fullfile(IMAGES_FOLDER, 'Phantom_tensor', STI_FILE{n_tensor});        
        ModelParams.R1map_file = fullfile(IMAGES_FOLDER, 'Maps', 'r1.nii.gz');
        ModelParams.R2starmap_file = fullfile(IMAGES_FOLDER, 'Maps', 'r2s.nii.gz');
        ModelParams.M0map_file = fullfile(IMAGES_FOLDER, 'Maps', 'm0.nii.gz');
        ModelParams.Segmentation_file = fullfile(IMAGES_FOLDER, 'Masks', 'atlas_mask.nii.gz');
        ModelParams.BrainMask_file = fullfile(IMAGES_FOLDER, 'Masks', 'brain_mask.nii.gz');
        
        ModelParams.FA_file = fullfile(IMAGES_FOLDER, 'Diffusion_data', 'DTI_reg_files', 'dti_reg_FA.nii.gz');
        if n_tensor == 1
            ModelParams.V1_file = fullfile(IMAGES_FOLDER, 'Diffusion_data', 'DTI_reg_files', 'V1_filtered.nii.gz');
        else
            ModelParams.V1_file = fullfile(IMAGES_FOLDER, 'Susceptibility_data', 'V1_filtered.nii.gz');
        end
        ModelParams.PhaseOffset = fullfile(SimParams.sim_folder, 'phase_offset.nii.gz');
        
        ModelParams.phi_12_root = fullfile(SimParams.sim_folder, '12_orientations');
        ModelParams.phi_12 = fullfile(ModelParams.phi_12_root, 'phi_12_orientations.nii.gz');
        
        ModelParams.num_orientations = '12_orientations';
        ModelParams.angles = 'angles_2.mat';
        ModelParams.local_phase_name = 'bulk_phase.nii.gz';
        ModelParams.out_magn_name = 'magn.nii.gz';
        ModelParams.out_phase_name = 'phase_unwrapp.nii.gz';
        ModelParams.in_phase = 'phase.nii.gz';
        
        %% Data Simulation        
        simulate_data(ModelParams, SeqParams, SimParams)
        pre_process_tensor(ModelParams, SimParams, SeqParams);
        merge_phase(ModelParams);
    end
end
diary off
