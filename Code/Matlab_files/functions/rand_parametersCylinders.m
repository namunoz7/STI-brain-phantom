function [r_cylinders,n_cylinders,v_theta,v_psi,centers] = rand_parametersCylinders(n_maxCylinders,N)
% Generates random parameters of the cylinder
% Inputs:
%   n_maxCylinders : Maximum number of cylinders in the image
%   N : matrix size of the image
% Outputs:
%   r_cylinders : vector containing the radius of the cylinders
%   n_cylinders : number of cylinders in the image
%   v_theta,v_psi : angular values of the cylinder
%   centers : center of the cilinders.


n_maxTheta = 180.0;
n_maxPsi = 180;
n_cylinders = randi(n_maxCylinders);

max_rCylinders = min(N)/20;
max_rCylinders = round(max_rCylinders);
r_cylinders = randi([5,max_rCylinders],[1,n_cylinders]);

v_theta = n_maxTheta*rand([n_cylinders,1]) - 0.5*n_maxTheta;
v_theta = v_theta*pi/180;

v_psi = n_maxPsi*rand([n_cylinders,1]) - 0.5*n_maxPsi;
v_psi = v_psi*pi/180;

tmp_1 = min(N)/4;
tmp_2 = max(N) - tmp_1;

centers = randi([tmp_1,tmp_2],[3,n_cylinders]);

end