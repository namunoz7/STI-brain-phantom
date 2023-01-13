function [ sig, vox_out,sigHR] = DataSimulation(SequenceParam,TissueParam)
% sig = DataSimulation(SequenceParam,TissueParam);
% data s
% config is a structure with all sort of parameters
%     SequenceParam.TR;                     repetition time (in seconds);
%     SequenceParam.TE;                     echo time (in seconds);
%     SequenceParam.theta;                  flip angle in degrees;
%     SequenceParam.res;                    resolution in mm;
%     TissueParam.M0;                      Water Concentration (arbitrry Units); 
%     TissueParam.R1;                      longitudinal relaxation rate (1/s); 
%     TissueParam.R2star;                  apparent transverse relaxation rate (1/s);
%     TissueParam.field;                   field (1/s); - we could change
%     this to have the info of the 
%     TissueParam.PhaseOffset;             PhaseOfsett at TE=0;
%     TissueParam.res        ;             resolution in mm;

%% reading sequence parameters

if isfield(SequenceParam,'TR')
    TR=SequenceParam.TR;
else
    TR=1;
end

if isfield(SequenceParam,'TE')
    TE=SequenceParam.TE;
else
    TE=30e-3;
end

if isfield(SequenceParam,'theta')
    theta=SequenceParam.theta;
else
    theta=90;
end

if isfield(SequenceParam,'res')
    vox=SequenceParam.res;
else
    vox=[1 1 1];
end

%% reading tissue model parameters

if isfield(TissueParam,'R1')
    R1=TissueParam.R1;
else
    R1=1;
end

if isfield(TissueParam,'R2star')
    R2star=TissueParam.R2star;
else
    R2star=50;
end

if isfield(TissueParam,'M0')
    M0=TissueParam.M0;
else
    M0=1;
end



if isfield(TissueParam,'field')
    field=TissueParam.field;
else
    field=0;
end

if isfield(TissueParam,'PhaseOffset')
    PhaseOffset=TissueParam.PhaseOffset;
else
    PhaseOffset=0;
end

if isfield(TissueParam,'res')
    res=TissueParam.res;
else
    res=[1 1 1];
end

% keyboard

%% calculating signal at full resolution
sigHR = M0.*exp(1i.* (2 * pi * field * TE + PhaseOffset)) .*exp(-TE.*R2star)...
    .*(1-exp(-TR.*R1)).*sind(theta)./(1-cosd(theta).*exp(-TR.*R1));

sigHR(isnan(sigHR))=0;
% resampling, keeping the same field of view 
dims_start = size(R1);
dims_end = round(dims_start.*res./vox);
dims_end(dims_end>dims_start) = dims_start(dims_end>dims_start);
vox_out = res.*dims_start./dims_end;
 sig = ifftn(ifftshift(crop(fftshift(fftn(sigHR)),dims_end)));
% sig =
% ifftshift(ifftn(ifftshift(crop(fftshift(fftn(fftshift(sigHR))),dims_end))));%del
% sig =
% ifftshift(ifftn(fftshift(crop(fftshift(fftn(ifftshift(sigHR))),dims_end))));%this
% is formaly correct but messes up the phase when downsampling, 

ShowResults=0;
if ShowResults==1
    figureJ();
    subplot1(2,2,'Gap',[0.02 0.02]);
    subplot1(1);
    Orthoview(abs(sigHR));
    title('High Res absolute image');
    subplot1(2);
    Orthoview(angle(sigHR));
    title('High Res phase image');
    subplot1(3);
    Orthoview(abs(sig));
    title('Low Res absolute image');
    subplot1(4);
    Orthoview(angle(sig));
    title('Low Res phase image');
end


