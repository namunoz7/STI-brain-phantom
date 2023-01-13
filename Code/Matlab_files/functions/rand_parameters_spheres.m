function [r_spheres,n_spheres,center,Xin,Xout] = rand_parameters_spheres(n_maxSpheres,N)
% Generate random radius and number of spheres 
% Parameters:
%   n_maxSpheres : maximum number of spheres
%   N : size of the matrix{
% Outputs:
%   r_spheres : vector with the of all the spheres in the image.
%   n_spheres : number of spheres in the image
%   center : center of the spheres in the image
%   Xi : internal susceptibility
%   Xout : external susceptibility

min_rSpheres = 3;
n_minSpheres = 3;
Xin = 3e-6;
Xout = 1e-6;

max_rSpheres = min(N)/10;
max_rSpheres = round(max_rSpheres);
n_spheres= randi([n_minSpheres,n_maxSpheres]);
r_spheres = randi([min_rSpheres,max_rSpheres],[1,n_spheres]);

tmp_min_radius = max(r_spheres);
tmp_max_radius = min(N) - tmp_min_radius;

center = tmp_min_radius + (tmp_max_radius - tmp_min_radius)...
    *rand([3,n_spheres]);
 
end