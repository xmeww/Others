time = (1:length(maskSND))/Sd.freqBase;

me_MaskSnd = [];
for iFrequency = 1:26;
    nrPlayed = 0;
    e_MaskSnd = [];
    for iTrial = 1:S.nrTrials
        load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/Sub1/S1/maskSND_T' num2str(iTrial) '.mat']);
        if ~o.timing.audMaskTime(iTrial, iFrequency) == 0
            nrPlayed  = nrPlayed + 1;
            time2trig = o.timing.triggerOnsets(1, iTrial) - o.timing.audMaskTime(iTrial, iFrequency);
            samp2trig = find((time > time2trig - 0.5) & (time < time2trig + 1.5));
            e_MaskSnd(nrPlayed, :) = maskSND(iFrequency, samp2trig);
        end
    end
    me_MaskSnd(iFrequency,:) = mean(e_MaskSnd,1);
end