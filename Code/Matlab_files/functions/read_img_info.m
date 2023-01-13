function [TE, voxel_size, img_size, B0, gyro] = read_img_info(folder_base_name, dicom_base_name)
f = filesep;
dicom_files = dir([folder_base_name, dicom_base_name]);
% READ FIRST FILE TO GET NUMBER OF ECHOES
first_dicom_file = dicominfo([dicom_files(1).folder, f, dicom_files(1).name]);
n_echoes = first_dicom_file.EchoTrainLength;
n_slices = length(dicom_files);
TE = zeros(1, n_echoes);
for m = 1:n_echoes
    new_vol = 1 + (m-1)*n_slices/n_echoes;
    actual_dicom = dicominfo([dicom_files(new_vol).folder, f, dicom_files(new_vol).name]);
    TE(m) = actual_dicom.EchoTime;
end
voxel_size = single([actual_dicom.PixelSpacing; actual_dicom.SliceThickness]');
img_size = single([actual_dicom.Columns, actual_dicom.Rows, n_slices/n_echoes]);
B0 = actual_dicom.MagneticFieldStrength;
gyro = 42.58;
end