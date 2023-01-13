function chi_res = sti_lsqr(phase_tissue, direction_field, mask)
N = size(mask);
[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1,-N(1)/2:N(1)/2-1,-N(3)/2:N(3)/2-1);      % k-space grid
kx = fftshift(kx);      ky = fftshift(ky);      kz = fftshift(kz);

N_direction = size(phase_tissue, 4);        % no of directions

Phase_tissue = zeros(size(phase_tissue));
for ind = 1:N_direction
    Phase_tissue(:,:,:,ind) = fftn(phase_tissue(:,:,:,ind));
end

param = [];
param.SS = N;
param.N_direction = N_direction;
param.H_Matrix = direction_field;
param.kx = kx;      param.ky = ky;      param.kz = kz;
param.k2 = kx.^2 + ky.^2 + kz.^2;

lsqr_tol = 5e-3;                            % LSQR tolerance
lsqr_iter = 30;                             % no of LSQR iterations    

tic
    [res, flag, relres, iter] = lsqr(@apply_STI, Phase_tissue(:), lsqr_tol, lsqr_iter, [], [], [], param);  
    disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
toc

Fchi_res = reshape(res, [N,6]);             % susceptibility tensor in k-space

chi_res = zeros(size(Fchi_res));
for ind = 1:6
    chi_res(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
end