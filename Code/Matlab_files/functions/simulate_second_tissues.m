function [Chimap3] = simulate_second_tissues(Mask, indexes, R2s, R1, chiref, R2star2chi, R1only2chi, R12chi, Chimap3, Chimap4, smoothThresh, smoothThresh2, NSTD, ProbabilityAccumulated)
chitemp=...
    R2star2chi*(R2s(indexes)-mean(R2s(indexes)))+...
    R12chi*(R1(indexes)-mean(R1(indexes)));

modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;

Chimap3(indexes)=Chimap3(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes))...
    .*(chiref+chitemp);

% only based on R1
% chitemp=R1only2chi*(R1(indexes)-mean(R1(indexes)));
% 
% modulation_std = std(chitemp(Mask(indexes)>smoothThresh));
% modulation_mean = mean(chitemp(Mask(indexes)>smoothThresh));
% chitemp(and(Mask(indexes)>smoothThresh2,chitemp>modulation_mean+NSTD*modulation_std))=modulation_mean;
% Chimap4(indexes)=Chimap4(indexes)+(Mask(indexes)-ProbabilityAccumulated(indexes)).*(chiref+chitemp);
end