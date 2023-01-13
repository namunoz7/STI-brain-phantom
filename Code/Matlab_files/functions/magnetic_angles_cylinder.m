function [v_thetaMagn,v_psiMagn] = magnetic_angles_cylinder(n_orientations)
% Generate the different angles of the orientation of the magnetic field
% Inputs:
%   n_orientation : number of different orientations of the magnetic field
% Outputs:
%   v_thetaMagn : azimuth angles with respect to the axis of the magnetic 
%   field
%   v_psiMagn : rotation angles with respect to the axis of the magnetic
%   field

big = 90;
medium = 45;
small = 0;

theta_max = 90;
psi_max = 180;

tmp_theta = zeros(6,1);
tmp_psi = zeros(6,1);

tmp_theta(1) = big*pi/180;
tmp_psi(1) = small*pi/180;

tmp_theta(2) = big*pi/180;
tmp_psi(2) = medium*pi/180;

tmp_theta(3) = medium*pi/180;
tmp_psi(3) = medium*pi/180;

tmp_theta(4) = big*pi/180;
tmp_psi(4) = big*pi/180;

tmp_theta(5) = medium*pi/180;
tmp_psi(5) = big*pi/180;

tmp_theta(6) = small*pi/180;
tmp_psi(6) = small*pi/180;

v_theta = 2*theta_max*rand(n_orientations,1) - theta_max;
v_theta = v_theta*pi/180;

v_psi = 2*psi_max*rand(n_orientations,1) - psi_max;
v_psi = v_psi*pi/180;

v_thetaMagn = [tmp_theta;v_theta];
v_psiMagn = [tmp_psi;v_psi];

end

