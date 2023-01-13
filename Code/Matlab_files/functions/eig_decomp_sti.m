function [Chi_eig, Chi_vec] = eig_decomp_sti(chi_res, mask)
N = size(mask);
Chi_res = zeros([N, 3, 3]);
Chi_res(:,:,:,1,1) = chi_res(:,:,:,1);
Chi_res(:,:,:,1,2) = chi_res(:,:,:,2);
Chi_res(:,:,:,2,1) = chi_res(:,:,:,2);
Chi_res(:,:,:,1,3) = chi_res(:,:,:,3);
Chi_res(:,:,:,3,1) = chi_res(:,:,:,3);
Chi_res(:,:,:,2,2) = chi_res(:,:,:,4);
Chi_res(:,:,:,2,3) = chi_res(:,:,:,5);
Chi_res(:,:,:,3,2) = chi_res(:,:,:,5);
Chi_res(:,:,:,3,3) = chi_res(:,:,:,6);

mask_tensor = mask(:);

Chi_res = permute(Chi_res, [4,5,1,2,3]);
chi_tensor = reshape(Chi_res, [3,3, numel(mask_tensor)]);

Chi_eig = zeros(numel(mask_tensor),3);
Chi_vec = zeros(numel(mask_tensor), 3, 3);

tic
for v = 1:length(mask_tensor)
    if mask_tensor(v) ~= 0
        [V,D] = eig(chi_tensor(:,:,v));
        Chi_eig(v,:) = diag(D)';
        Chi_vec(v,:,:) = V;
    end
end
toc

Chi_eig = reshape(Chi_eig, [N, 3]);
Chi_vec = reshape(Chi_vec, [N, 3, 3]);
end