clearvars; clc; close all
cd ../
addpath(genpath('functions/'))
principalFolder = pwd;
f = filesep;
addpath(genpath('../mritools_Linux_3.4.3/'))
diary Simulate_data.txt

IMAGES_FOLDER = fullfile('..', '..', '..', 'Imagenes', 'Phantom_real_data', 'Orientations');
MODE_FOLDER = {'Susceptibility_sti', 'Diffusion_sti'};
OUT_FOLDER = 'QSM_reconstructions';
PEV_NAME = 'V1_filtered.nii.gz';
SUSCEPTIBILITY_FOLDER = fullfile(IMAGES_FOLDER, '..', 'Susceptibility_data');
DIFFUSION_FOLDER = fullfile(IMAGES_FOLDER, '..', 'Diffusion_data', 'DTI_reg_files');

N_DATA = 245; % Number of data points in the roi
mean_qsm = zeros(15, 2);
std_qsm = zeros(15, 2);
x_angles = linspace(0, 90, 15);

B0_dir = single([0, 0, 1])';

for n_mode = 1:length(MODE_FOLDER)
    disp('=================================')
    disp(['Simulating ', MODE_FOLDER{n_mode}, ' mode....'])
    sim_folder = fullfile(IMAGES_FOLDER, MODE_FOLDER{n_mode}, OUT_FOLDER);
    if n_mode == 1
        pev_file = fullfile(SUSCEPTIBILITY_FOLDER, PEV_NAME);
    else
        pev_file = fullfile(DIFFUSION_FOLDER, PEV_NAME);
    end
    nii_pev = load_untouch_nii(pev_file);
    fiber_dir = squeeze(nii_pev.img(75, 75, 47, :));
    angle_dephase = abs(90 - acosd(dot(B0_dir, fiber_dir)));
    disp(['Dephase angle: ', num2str(angle_dephase)])
    num_orientations = length(dir(fullfile(sim_folder, '*.nii.gz')));
    data_roi = zeros(N_DATA, num_orientations);
    for n_orientation = 1:num_orientations
        qsm_name = ['QSM_', num2str(n_orientation), '_orientation.nii.gz'];
        qsm_file = fullfile(sim_folder, qsm_name);
        nii_qsm = load_untouch_nii(qsm_file);
        qsm_roi = nii_qsm.img(72:78, 72:78, 45:49);
        data_roi(:, n_orientation) = qsm_roi(:);
        mean_qsm(n_orientation, n_mode) = mean(qsm_roi(:));
        std_qsm(n_orientation, n_mode) = std(qsm_roi(:));
    end
    sin_eq_2 = ['a*sind(x+', num2str(angle_dephase), ').^2+b'];
    [f2, gof] = fit(x_angles', mean_qsm(:,n_mode), sin_eq_2);
    y2 = f2.a*sind(x_angles + angle_dephase).^2+f2.b;
    fig = figure('Position', [500, 500, 320, 270], 'Units', 'centimeters'); 
    plot(x_angles, 100*mean_qsm(:,n_mode), 'r*')
    hold on
    plot(x_angles, 100*y2, 'k')
end