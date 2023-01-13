function sigHR = Tensor_DataSimulation(SequenceParam,TissueParam)
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

%% reading tissue model parameters

if isfield(TissueParam,'R1')
    R1=single(TissueParam.R1);
else
    disp('No R1 image loaded!!!')
    R1=1;
end

if isfield(TissueParam,'R2star')
    R2star=single(TissueParam.R2star);
else
    disp('No R2s image loaded!!!')
    R2star=50;
end

if isfield(TissueParam,'M0')
    M0=TissueParam.M0;
else
    disp('No M0 image loaded!!!')
    M0=1;
end

if isfield(TissueParam,'field')
    field=TissueParam.field;
else
    disp('No field loaded!!!')
    field=0;
end

if isfield(TissueParam,'PhaseOffset')
    PhaseOffset=TissueParam.PhaseOffset;
else
    disp('No phase offset loaded!!!')
    PhaseOffset=0;
end

% keyboard

%% calculating signal at full resolution
T1 = exp(-TR.*R1);
T2 = exp(-TE.*R2star);
magn = M0.*(1 - T1)*sind(theta)./(1 - cosd(theta).*T1).*T2;
phase = field*TE + PhaseOffset;
sigHR = magn.*exp(1i*(phase));
% sigHR = M0.*(1-exp(-TR.*R1)).*sind(theta)./(1-cosd(theta).*exp(-TR.*R1))...
%     .*exp(-TE.*R2star).*exp(1i*(field*TE + PhaseOffset));

sigHR(isnan(sigHR))=0;
