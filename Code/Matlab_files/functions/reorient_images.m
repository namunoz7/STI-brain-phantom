function reorient_images(ref_magn_file, in_magn_file, in_phase_file, out_phase_file, mat_rot)
options_flirt = '-cost mutualinfo -searchcost mutualinfo -interp spline';
find_mat_rot = ['flirt ', options_flirt, ' -in ', in_magn_file, ' -ref ', ref_magn_file, ' -omat ', mat_rot];
disp(find_mat_rot)
system(find_mat_rot)
command_rotate = ['flirt ', options_flirt, ' -applyxfm -init ', mat_rot, ' -in ', in_phase_file, ' -ref ', ref_magn_file, ' -out ', out_phase_file];
disp(command_rotate)
system(command_rotate)
disp([out_phase_file, ' created']);
end