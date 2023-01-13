function merge_phase(ModelParams)
% clearvars; close all
% cd ../
% addpath(genpath('functions/'))
% % cd ../
% principalFolder = pwd;
% cd ../../Imagenes/Phantom/Simulated_images/
% imagesFolder = SimParams.out_folder;
% addpath(genpath([imagesFolder, '/func/']))

f = filesep;
%%
folder_images = ModelParams.num_orientations;  % '12_orientations';
bulk_phase = ModelParams.local_phase_name;  % 'bulk_phase.nii.gz';
actual_model = ModelParams.phi_12_root;
microstructure_folder = ModelParams.out_folder;

disp(folder_images)
n_orientations = length(dir([microstructure_folder, folder_images, '/Orientation_*']));
phase_file = '';
for n = 1:n_orientations
    actual_folder = fullfile(actual_model, microstructure_folder, ['Orientation_', num2str(n)]);
    phase_file = [phase_file, actual_folder, bulk_phase, ' '];
end
phase_file = phase_file(1:end-1);
out_file = [actual_model, f, bulk_phase, ' '];
% command = ['fsl5.0-fslmerge -t ', out_file, phase_file];
command = ['fslmerge -t ', out_file, phase_file];
disp(command)
system(command)

end
