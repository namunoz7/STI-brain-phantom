function pre_process_tensor(ModelParams, SimParams, SeqParams)
% clearvars; close all
% cd ../
% addpath(genpath('functions/'))
% % cd ../
% principalFolder = pwd;
% addpath(genpath('../SEGUE_28012021/'))
mask_hr = load_untouch_nii(ModelParams.BrainMask_file);  % load_nii('../../Imagenes/Phantom/masks/mask_tensor3.nii.gz');
mask_hr = mask_hr.img;
% cd ../../Imagenes/Phantom/Simulated_images/
mu_folder = ModelParams.out_folder;
folder_images = SimParams.sim_folder;
num_orientations = ModelParams.num_orientations;  % '12_orientations';
out_magn_name = ModelParams.out_magn_name;  % 'magn.nii.gz';
out_phase_name = ModelParams.out_phase_name;  % 'phase_unwrapp.nii.gz';
in_phase = ModelParams.in_phase;  % 'phase.nii.gz';
local_phase_name = ModelParams.local_phase_name;  % 'bulk_phase.nii.gz';
n_orientation = SeqParams.n_orietations;  % 12;
TE = SeqParams.TE;  % [4, 12, 20 ,28];

% Parameters ROMEO
params = [];
params.TE = TE; % required for multi-echo
params.mask = mask_hr;
params.calculate_B0 = false;
params.phase_offset_correction = 'bipolar';
params.additional_flags = '-vgiQ';
params.TE = TE;
params.calculate_B0 = false;
params.phase_offset_correction = 'on';
params.additional_flags = '-i';

disp([num2str(n_orientation), ' orientations'])
for n = 1:n_orientation
    disp(['...... Orientation ', num2str(n)])

    actual_orientation = fullfile(folder_images, num_orientations, mu_folder, ...
        ['Orientation_', num2str(n)]);
    params.output_dir = fullfile(actual_orientation, 'Tmp_phase'); % if not set pwd() is used
    if ~isfolder(params.output_dir)
        mkdir(params.output_dir);
        disp([params.output_dir, ' created'])
    end

    actual_magn = fullfile(actual_orientation, out_magn_name);
    actual_phase = fullfile(actual_orientation, in_phase);
    out_phase = fullfile(actual_orientation, out_phase_name);
    out_local_phase = fullfile(actual_orientation, local_phase_name);    
    magn = load_untouch_nii(actual_magn);
    voxel_size = magn.hdr.dime.pixdim(2:4);
    magn = magn.img;
    magn = (magn - min(magn(:)))./max(magn(:));

    params.mag = magn;
    params.voxel_size = voxel_size;

    disp('...... Unwrapping phase images')
    phase = load_untouch_nii(actual_phase);
    phase.img = ROMEO(phase.img, params);
    disp('Saving unwrapped phase')
    save_untouch_nii(phase, out_phase)
    disp([out_phase, ' saved'])

    disp('... Fitting phase echoes')
    [~, deltaB, ~] = ilinearLSQR(phase.img, TE, magn);
    phase.img = deltaB;
    phase.hdr.dime.dim(1) = 3;
    phase.hdr.dime.dim(5) = 1;
    save_untouch_nii(phase, out_local_phase)
    disp([out_local_phase, ' saved'])
end
end

