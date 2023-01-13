function WeightingBetweenModels = apply_gradient(HighGradMask, RealFreqMap, FinalSegm)
% Apply gradient to smooth R2* map in high intensity regions
% Params: HighGradMask, RealFreqMap, FinalSegm
% Return: WeightingBetweenModels
disp('Applying 3D gradient field')
[fieldgradient(:,:,:,1), fieldgradient(:,:,:,2), fieldgradient(:,:,:,3)]=gradient_3D(RealFreqMap.img,[],0);
fieldgradient=(sos(fieldgradient,4));
fieldgradient(FinalSegm>=11)=0;

fieldgradient=fieldgradient.*single(HighGradMask.img);
%     %define acceptable threshold ???
% thresholdgrad=prctile(fieldgradient(fieldgradient~=0),[10 90]); % is this a good threshold, check on the fv of the fiedlgradient
% looking at the figure 0.05 seems to be the gradient close to the edges
% looking at the figure 0.2 seems to be a too large gradient

fieldgradient(and(fieldgradient~=0,fieldgradient>0.2))=0.2;
fieldgradient(and(fieldgradient~=0,fieldgradient<0.05))=0.05;
fieldgradient(fieldgradient~=0)=(fieldgradient(fieldgradient~=0)-0.05)/(0.2-0.05);
WeightingBetweenModels=cos(fieldgradient);
end