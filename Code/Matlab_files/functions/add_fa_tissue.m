% function Chimap3 = add_fa_tissue(Mask_smooth, LG_mask_index, indexes, FA, FAtoChi, chiref, Chimap3, smoothThresh, smoothThresh2, NSTD)
function Chimap3 = add_fa_tissue(Mask_smooth, indexes, FA, FAtoChi, chiref, Chimap3, smoothThresh, smoothThresh2, NSTD)
tissue=...
    FAtoChi*(FA-mean(FA(indexes)));

% trying to find a way to avoid creating outlyers because of errors on the R2* and R1 maps...
modulation_std = std(tissue(Mask_smooth>smoothThresh));
modulation_mean = mean(tissue(Mask_smooth>smoothThresh));
% makes all the outliers have the average value of the segmented tissue
tissue(and(Mask_smooth>smoothThresh2,abs(tissue-modulation_mean)>NSTD*modulation_std))=modulation_mean;

Chimap3=Chimap3+Mask_smooth.*(chiref+1*tissue);

end