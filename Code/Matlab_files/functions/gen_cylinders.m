function gen_cylinders(s_cylinders,n_cylinders,s_theta,s_psi,center,FOV,N)


kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);

delta_kx = FOV(1)/N(1);
delta_ky = FOV(2)/N(2);
delta_kz = FOV(3)/N(3);

[kxx,kyy,kzz] = meshgrid(kx,ky,kz);

Rx_theta  = [1 0 0 0;
    0 cos(s_theta) -sin(s_theta) 0;
    0 sin(s_theta) cos(s_theta) 0;
    0 0 0 1];
% Rotation with psi

Rz_psi  = [cos(s_psi) -sin(s_psi) 0 0;
    0 1 0 0;
    -sin(s_psi) 0 cos(s_psi) 0;
    0 0 0 1];

T_2 = [1 0 0 -Nx;
    0 1 0 -Ny;
    0 0 1 -Nz;
    0 0 0 1];

T_1 = [1 0 0 Nx;
    0 1 0 Ny;
    0 0 1 Nz;
    0 0 0 1];

kxx = reshape(kxx,1,numel(xx));
kyy = reshape(kyy,1,numel(yy));
kzz = reshape(kzz,1,numel(zz));
unos = ones(1,numel(zz));

original_positions = [kxx;
    kyy;
    kzz;
    unos];

new_positions = T_2*Rz_psi*Rx_theta*T_1*original_positions;

kxx = reshape(new_positions(1,:),matrix_size);
kyy = reshape(new_positions(2,:),matrix_size);

kxx = single(kxx)*delta_kx;
kyy = single(kyy)*delta_ky;


end