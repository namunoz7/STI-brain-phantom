function process_dti_files(params)
% Process diffusion tensor images. Filers the tensor, generates the
% eigen-decomposition, i.e. eigen-values and eigen-vectors, and finally
% generates the fractional anisotropy to add texture to the phantom

dti_tensor = load_nii(params.tensor_file);
nii_v1 = load_nii(params.v1_file);
nii = load_nii(params.s0_file);
mask = nii.img>0;
disp('... Filtering diffusion tensor')
for nn = 1:6
    disp(nn)
    tmp=Wavedec3Denoising(dti_tensor.img(:, :, :, nn),0.0014/87,8,'db2','soft',double(mask));
    dti_tensor.img(:, :, :, nn) = abs(real(tmp)); % same as for segmentation purposes
end

if isfile(params.out_eig)
    disp('... Loading eigen decomposition')
    nii_eig = load_nii(params.out_eig);
    eig = nii_eig.img;
    clear nii_eig
    disp(['... ', params.out_eig, ' loaded'])
else
    disp('...Generating eigen-decomposition')
    [eig_vec, eig] = eig_dti(dti_tensor.img, mask);
    
    disp('...Saving eigen decomposition images')
    nii_v1.img = eig;
    save_nii(nii_v1, params.out_eig)
    disp(['... ', params.out_eig, ' saved'])

    nii_v1.img = eig_vec(:, :, :, :, 1);
    save_nii(nii_v1, params.out_v1)
    disp(['... ', params.out_v1, ' saved'])    
    nii_v1.img = eig_vec(:, :, :, :, 2);
    save_nii(nii_v1, params.out_v2)
    disp(['... ', params.out_v2, ' saved'])    
    nii_v1.img = eig_vec(:, :, :, :, 3);
    save_nii(nii_v1, params.out_v3)
    disp(['... ', params.out_v3, ' saved'])
end
MD = mean(eig, 4);
FA = sqrt(1.5).*(vecnorm(MD-eig,2,4)./vecnorm(eig,2,4));
nii.img = FA;
save_nii(nii, params.out_fa)
disp(['... ', params.out_fa, ' saved'])

end