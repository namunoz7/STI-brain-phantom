function nii_mask = rotate_mask(ref_magn_file, in_magn_file, in_mask_file, out_mask_file, mat_rot)
% Rotate the mask that is in line with the main magnetic field
% param: ref_magn_file, in_magn_file, in_mask_file, out_mask_file, mat_rot
% out: nii_mask
options_fsl = 'fsl5.0-flirt -cost mutualinfo -searchcost mutualinfo -interp spline';
find_mat_rot = [options_fsl, ' -in ',in_magn_file, ' -ref ', ref_magn_file, ' -omat ', mat_rot];
disp(find_mat_rot)
system(find_mat_rot)
command_rotate = ['fsl5.0-flirt -interp nearestneighbour -cost mutualinfo -searchcost mutualinfo -applyxfm -init ', mat_rot,...
    ' -in ', in_mask_file, ' -ref ', ref_magn_file, ' -out ', out_mask_file];
disp(command_rotate)
system(command_rotate)
nii_mask = load_untouch_nii(out_mask_file);
end