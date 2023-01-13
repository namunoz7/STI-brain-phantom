function off_resonance_field = tensor_phase_processing(phase_use, brain_mask, TE, N, voxel_size, B0, params)
Rvector = 30:-1:1;
% delta_TE = (TE(2) - TE(1))/1000; % Time in seconds
% phase_scale = 2*pi*gamma*B0*delta_TE;
phase_scale = 2*pi*B0;
% unwL = zeros([N, length(TE)]);
disp('=============================')
[unw_phase, ~] = ROMEO(phase_use, params);
unw_phase = (TE(1)*unw_phase(:, :, :, 1) - TE(2)*unw_phase(:, :, :, 2))./(TE(2) - TE(1));
disp('Removing background');
disp('LBV ...')
field_LBV = LBV(unw_phase,brain_mask,N,voxel_size);
disp('vSharp ...');
[~, RDF, ~]=  m1_bb_vSHARP(double(field_LBV),double(brain_mask),double(Rvector));
off_resonance_field = RDF/phase_scale;
disp('=============================')
% field = magn_use.*exp(1i*single(RDF));
% disp('Fitting image echoes')
% off_resonance_field = Fit_ppm_complex_TE(field , TE);
% off_resonance_field = off_resonance_field/phase_scale;
end