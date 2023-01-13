function new_dti = filter_dti(params)

nii_dti = load_nii(params.dti_file);
new_dti = nii_dti;

disp('Wavelet filtering')
disp('Filtering dti tensor')
for nn = 1:6
    disp(['Component ', num2str(nn)])
    actual_img = nii_dti.img(:, :, :, nn);
    tmp=Wavedec3Denoising(actual_img,2.6101e-05,8,'db2','soft',actual_img);
    new_dti.img(:,:,:,nn) = abs(real(tmp)); % same as for segmentation purposes
end
% dti_filtered = new_dti.img;
save_nii(new_dti, params.out_dti)

% new_dti.img = new_dti.img - nii_dti.img;
% save_nii(new_dti, params.difference_file)

end