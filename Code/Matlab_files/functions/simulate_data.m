function simulate_data(ModelParams, SeqParams, SimParams)
% ModelParams.R1map_file
% ModelParams.R2starmap_file
% ModelParams.M0map_file
% ModelParams.Segmentation_file
% ModelParams.Chimap_file
% ModelParams.BrainMask_file

% SeqParams.TR          Repetition time in secs
% SeqParams.TE          Echo time in secs
% SeqParams.FlipAngle   flip angle in degrees

% SimParams.PhaseOffset multiplier term of a quadratic phase over the brain
%                       0 no phase offset; 2pi phase difference inside
%                       brain mask
% SimParams.Shimm       boolean 0 if no additional shimms are applied 1 if
%                       2nd order shimms are applied inside brain mask
% SimParams.BrainExtracted  boolean 0 makes the simulation with whole head
%                       susceptibility, 1 makes only head
% SimParams.Res         Resolution of output
% SimParams.Output_dir  Output Directory

% modeln refers to which susceptibility model is being used
% SimulateData is wheather the simulation part is bing done or not

%% Data needed to be loaded for simulation
disp('Loading images ...')
M0 = load_untouch_nii(ModelParams.M0map_file);
R1 = load_untouch_nii(ModelParams.R1map_file);
R2star = load_untouch_nii(ModelParams.R2starmap_file);
v1_dti = load_untouch_nii(ModelParams.V1_file);
FA_nii = load_untouch_nii(ModelParams.FA_file);
WM_mask = load_untouch_nii(ModelParams.Segmentation_file);
WM_mask = WM_mask.img == 3;
nii_tmp = R2star;
nii_tmp.hdr.dime.dim(1) = 4;
nii_tmp.hdr.dime.dim(5) = length(SeqParams.TE);
nii_tmp.hdr.dime.pixdim(5) = SeqParams.TR;

% Params
voxel_size = round(M0.hdr.dime.pixdim(2:4)*100)/100;
Brain = load_untouch_nii(ModelParams.BrainMask_file);
mask = logical(Brain.img);
TissueParam.M0 = (M0.img).*mask;                    %Water Concentration (arbitrry Units);
TissueParam.R1 = (R1.img).*mask;                    %longitudinal relaxation rate (1/s);
TissueParam.R2star = single(R2star.img).*mask;      %apparent transverse relaxation rate (1/s);
TissueParam.res = voxel_size;                     %resolution in mm;
TissueParam.mask = mask;
SequenceParam.TR = SeqParams.TR;                     %repetition time (in seconds);
SequenceParam.theta = SeqParams.FlipAngle;                     %flip angle in degrees;

if length(SimParams.Res)==1
    SequenceParam.res=[1 1 1]*SimParams.Res;             %resolution in mm;
elseif length(SimParams.Res)==3
    SequenceParam.res=SimParams.Res;             %resolution in mm;
else
    disp('you have not entered the resolution in the correct format, we [ 1 1 1 ] will be assumed')
    SequenceParam.res=[1 1 1 ];                             %resolution in mm;
end

%% Images in High Resolution
disp('Wavelet filtering')
disp('Filtering R1 map')
tmp=Wavedec3Denoising(TissueParam.R1,25/1000,8,'db2','verysoft',double(R1.img~=0));
TissueParam.R1 = abs(real(tmp)); % same as for segmentation purposes

disp('Filtering R2* map')
tmp=Wavedec3Denoising(TissueParam.R2star,45/10,8,'db2','soft',double(R2star.img~=0)); % same as for segmentation purposes
TissueParam.R2star = abs(real(tmp));

TissueParam.M0 = abs(TissueParam.M0);

dims = size(M0.img);
% phase_scale = SimParams.phase_scale;

disp('Phase offset ...') %PhaseOffset at TE=0;
if SimParams.PhaseOffset
    TissueParam.PhaseOffset = get_phase_offset(M0.img, dims, mask);
    disp('Saving phase offset ...')
    nii_tmp.img = TissueParam.PhaseOffset;
    save_untouch_nii(nii_tmp, ModelParams.PhaseOffset)
    disp([ModelParams.PhaseOffset, ' saved'])
    TissueParam.PhaseOffset = TissueParam.PhaseOffset.*mask;
else
    disp('No phase offset is calculated')
    TissueParam.PhaseOffset = 0;
end

%%
disp(['Simulating ', num2str(SimParams.number_orientations), ' orientations...'])
disp('... Loading local field')
field = load_untouch_nii(ModelParams.phi_12);
out_dir = fullfile(ModelParams.phi_12_root, ModelParams.out_folder);

field = single(field.img);

for nn_orient = 1:SimParams.number_orientations
    disp(['... Simulating ', num2str(nn_orient), ' orientation.'])
    TissueParam.field = field(:,:,:,nn_orient);   %field should be (Rads/s);
    Folder = fullfile(out_dir, ['Orientation_', num2str(nn_orient)]);
    if ~isfolder(Folder)
        mkdir(Folder)
        disp([Folder, ' created'])
    end

    angles = load(fullfile(ModelParams.phi_12_root, ModelParams.angles));
    if SimParams.microstructure
        new_v1 = get_v1_rotated(v1_dti.img, angles.vec_theta(nn_orient));

        microstructure = add_microstructure(new_v1, FA_nii.img/0.59, WM_mask);
        nii_tmp.img = microstructure;
        save_untouch_nii(nii_tmp, fullfile(Folder,'microstructure.nii.gz'))
        disp([fullfile(Folder,'microstructure.nii.gz') ' saved'])

        TissueParam.field = TissueParam.field + (microstructure*2*pi);
        nii_tmp.img = TissueParam.field;
        save_untouch_nii(nii_tmp, fullfile(Folder,'freq_shift.nii.gz'))
        disp([fullfile(Folder,'freq_shift.nii.gz') ' saved'])
    end

    out_magnHR = zeros([size(TissueParam.field), length(SeqParams.TE)]);
    out_phaseHR = zeros([size(TissueParam.field), length(SeqParams.TE)]);

    n_tmp = 1;
    for TE = SeqParams.TE
        disp(['...... ', num2str(n_tmp), ' echo'])

        % config is a structure with all sort of parameters
        SequenceParam.TE=TE; %echo time (in seconds);
        sigHR= Tensor_DataSimulation(SequenceParam,TissueParam);

        out_magnHR(:, :, :, n_tmp) = abs(sigHR).*mask;
        out_phaseHR(:, :, :, n_tmp) = angle(sigHR).*mask;
        n_tmp = n_tmp + 1;
    end
    disp('... Saving images')
    nii_tmp.img = out_magnHR;
    save_untouch_nii(nii_tmp, fullfile(Folder,'magn.nii.gz'))
    disp([fullfile(Folder,'magn.nii.gz'),' saved'])
    
    nii_tmp.img = out_phaseHR;
    save_untouch_nii(nii_tmp, fullfile(Folder,'phase.nii.gz'))
    disp([fullfile(Folder,'phase.nii.gz'), ' saved'])
end
end


% function lr_img = get_lr_image(hr_image, dims_end)
% lr_img = real(ifftn(ifftshift(crop(fftshift(fftn(hr_image)),dims_end))));
% end

function phi0 = get_phase_offset(M0, dims, mask)
[c , w ] = centerofmass(M0);
[y,x,z] = meshgrid((1:dims(2))-c(2), (1:dims(1))-c(1), (1:dims(3))-c(3));
temp = (x/w(1)).^2 + (y/w(2)).^2 + (z/w(3)).^2 ;
phi0 = - temp/(max(temp(mask==1))-min(temp(mask==1)))*pi;
end

function new_v1 = get_v1_rotated(v1, theta_n)
new_v1 = zeros(size(v1));
new_v1(:,:,:,1) = v1(:,:,:,1)*cos(theta_n) + v1(:,:,:,3)*sin(theta_n);
new_v1(:,:,:,2) = v1(:,:,:,2);
new_v1(:,:,:,3) = -v1(:,:,:,1)*sin(theta_n) + v1(:,:,:,3)*cos(theta_n);
end

function new_freq = add_microstructure(new_v1, FA_norm, WM_mask)
sin2_WM = (sum(new_v1(:,:,:,1:2).^2,4))./(sum(new_v1.^2,4));
new_freq = (-5*(sin2_WM - 2/3).*FA_norm - 3).*WM_mask;
new_freq(isnan(new_freq)) = 0;
end

