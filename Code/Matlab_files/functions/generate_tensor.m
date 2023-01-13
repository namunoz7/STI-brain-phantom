function nii_sti = generate_tensor(params)
% 
nii_eig = params.eigenValues;

nii_v1 = load_untouch_nii(params.V1_eigenvector);
nii_v2 = load_untouch_nii(params.V2_eigenvector);
nii_v3 = load_untouch_nii(params.V3_eigenvector);

N = nii_v1.hdr.dime.dim(2:4);
% voxel_size = nii_eig.hdr.dime.pixdim(2:4);

V = zeros([N, 3, 3]);
V(:, :, :, :, 1) = nii_v1.img;
V(:, :, :, :, 2) = nii_v2.img;
V(:, :, :, :, 3) = nii_v3.img;

disp('Filtering eigenvectors')
V = filt_eigenvectors(V);
disp('Saving filtered eigenvectors')

nii_v1.img = V(:, :, :, :, 1);
save_untouch_nii(nii_v1, params.Out_V1_filtered);
disp([params.Out_V1_filtered, ' saved'])

nii_v2.img = V(:, :, :, :, 2);
save_untouch_nii(nii_v2, params.Out_V2_filtered);
disp([params.Out_V2_filtered, ' saved'])

nii_v3.img = V(:, :, :, :, 3);
save_untouch_nii(nii_v3, params.Out_V3_filtered);
disp([params.Out_V3_filtered, ' saved'])

L = zeros([N, 3, 3]);
L(:, :, :, 1, 1) = nii_eig.img(:, :, :, 3);
L(:, :, :, 2, 2) = nii_eig.img(:, :, :, 2);
L(:, :, :, 3, 3) = nii_eig.img(:, :, :, 1);

disp('Constructing tensor')
tensor = mult_tensors(N, mult_tensors(N, V, L), permute(V, [1, 2, 3, 5, 4]));
disp('Saving tensor')
out_sti = zeros([N, 6]);
out_sti(:, :, :, 1) = tensor(:, :, :, 1, 1);
out_sti(:, :, :, 2) = tensor(:, :, :, 1, 2);
out_sti(:, :, :, 3) = tensor(:, :, :, 1, 3);
out_sti(:, :, :, 4) = tensor(:, :, :, 2, 2);
out_sti(:, :, :, 5) = tensor(:, :, :, 2, 3);
out_sti(:, :, :, 6) = tensor(:, :, :, 3, 3);
clear tensor

disp('... Creating nifti file');
nii_eig.img = out_sti;
nii_eig.hdr.dime.dim(5) = 6;
save_untouch_nii(nii_eig, params.OutChi)
disp([params.OutChi, ' saved'])

nii_sti = nii_eig;
end