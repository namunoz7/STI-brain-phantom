function [nii_v1, nii_v2, nii_v3, nii_eig] = eig_dti(params)
nii_dti = params.dti_filtered;
voxel_size = nii_dti.hdr.dime.pixdim(2:4);

dti = nii_dti.img;
N = size(dti, 1:3);

% Create mask
disp('Creating mask')
nii_S0 = load_nii(params.S0_file);
mask = nii_S0.img>0;
% mask = dti(:,:,:,1)>4e-4;
% for nn = 1:N(3)
%     mask(:, :, nn) = imfill(mask(:,:,nn), 8, 'holes');
% end

Dti_tensor = zeros([N, 3, 3]);
Dti_tensor(:,:,:,1,1) = dti(:,:,:,1);
Dti_tensor(:,:,:,1,2) = dti(:,:,:,2);
Dti_tensor(:,:,:,2,1) = dti(:,:,:,2);
Dti_tensor(:,:,:,1,3) = dti(:,:,:,3);
Dti_tensor(:,:,:,3,1) = dti(:,:,:,3);
Dti_tensor(:,:,:,2,2) = dti(:,:,:,4);
Dti_tensor(:,:,:,2,3) = dti(:,:,:,5);
Dti_tensor(:,:,:,3,2) = dti(:,:,:,5);
Dti_tensor(:,:,:,3,3) = dti(:,:,:,6);

mask_tensor = mask(:);

Dti_tensor = permute(Dti_tensor, [4,5,1,2,3]);
dti_tensor = reshape(Dti_tensor, [3,3, numel(mask_tensor)]);

vec_dti = zeros(prod(N), 3, 3);
eig_dti = zeros(prod(N),3);

disp('Generating eig decomposition')
for v = 1:length(mask_tensor)
    if mask_tensor(v) ~= 0
        [V,D] = eig(dti_tensor(:,:,v));
        [~,ind] = sort(diag(D), 'descend');
        Ds = D(ind, ind);
        Vs = V(:,ind);
        eig_dti(v,:) = diag(Ds)';
        vec_dti(v,:,:) = Vs;
    end
end

vec_dti = reshape(vec_dti, [N, 3, 3]);
eig_dti = reshape(eig_dti, [N, 3]);

V1 = vec_dti(:,:,:,:,1);
V2 = vec_dti(:,:,:,:,2);
V3 = vec_dti(:,:,:,:,3);

disp('Saving nifti files')
nii_v1 = make_nii(V1, voxel_size);
save_nii(nii_v1, params.out_v1)

nii_v2 = make_nii(V2, voxel_size);
save_nii(nii_v2, params.out_v2)

nii_v3 = make_nii(V3, voxel_size);
save_nii(nii_v3, params.out_v3)

nii_eig = make_nii(eig_dti, voxel_size);
save_nii(nii_eig, params.out_eig)

end