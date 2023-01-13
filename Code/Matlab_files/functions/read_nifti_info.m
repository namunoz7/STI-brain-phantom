function [TE, voxel_size, img_size, B0, gyro] = read_nifti_info(folder_base_name, nifti_base_name)
f = filesep;
json_files = dir([folder_base_name, nifti_base_name, '.json']);
nii_files = dir([folder_base_name, nifti_base_name, '.nii.gz']);
% READ FIRST FILE TO GET NUMBER OF ECHOES
str = fileread([json_files(1).folder, f, json_files(1).name]);
first_json_file = jsondecode(str);
info = load_untouch_nii([nii_files(1).folder, f, nii_files(1).name]);
B0 = first_json_file.MagneticFieldStrength;
img_size = info.hdr.dime.dim(2:4);
voxel_size = info.hdr.dime.pixdim(2:4);
n_echoes = first_json_file.EchoTrainLength;
TE = zeros(1, n_echoes);
for m = 1:n_echoes
    str = fileread([json_files(m).folder, f, json_files(m).name]);
    actual_json_file = jsondecode(str);
    TE(m) = actual_json_file.EchoTime * 1000; % ms
end
gyro = 42.58;
end