function [vox_out, dims_end] = calc_tensor_dims(SequenceParam,TissueParam)
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
if isfield(SequenceParam,'res')
    vox=SequenceParam.res;
else
    vox=[1 1 1];
end

if isfield(TissueParam,'res')
    res=TissueParam.res;
else
    res=[1 1 1];
end

if isfield(TissueParam, 'R1')
    R1 = TissueParam.R1;
else
    R1=1;
end

dims_start = size(R1);
dims_end = round(dims_start.*res./vox);
dims_end(dims_end>dims_start) = dims_start(dims_end>dims_start);
vox_out = res.*dims_start./dims_end;

end

