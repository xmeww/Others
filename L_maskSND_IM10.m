%% Generates 3 seconds masking sounds. 1 for each trial.

% subjects ID
SubID = 'Sub14';
% duration of mask in seconds
durationMask = 3.5;
% duration of beep in seconds
durationBeep = 0.3;
nTrials   = 250;
nSessions = 4;
% for all sessions
for iSession=1:nSessions
    % for all trials
    for iTrial = 1:nTrials;
        Sd.durMask  = durationMask;
        Sd.durBeep  = durationBeep;
        Sd.freqBase = 48000;
        Sd.AMmask   = 1;
        Sd.meanSOA  = 1050;
        % generate Mask
        [PossfreqEX IM SOAt] = IMsound_IM10(Sd.durMask, Sd.durBeep, Sd.freqBase, Sd.AMmask, Sd.meanSOA);
        maskSND = IM;
        % save
        eval(['save /Users/anette/Documents/MATLAB/Experiment2/Tones/' SubID '/S' num2str(iSession) '/maskSND_T' num2str(iTrial) ' maskSND'])
    end
end
Sd.PossFreq = PossfreqEX;


%% Generate Target tones
% AM modulation (= 1), no AM modulation (= 0)
AMmod = 1;

% target frequencies 
freqs  = [1222, 1583, 2049, 2654, 3437];
nFreqs = length(freqs);

% auditory sampling frequency
freqBase = Sd.freqBase;


for iFreq = 1:nFreqs
    % TARGET GENERATION
    samplesBeep = freqBase * durationBeep;
    time        = (0:samplesBeep-1) / freqBase; % time vector
    fc          = freqs(iFreq); % base frequency
    fcAM        = 40; % 40 Hz AM modulated
    beep        = sin(2*pi*fc*time);
    AM          = (1-cos(2*pi*fcAM*time))/2;
    if AMmod == 1
        beep = beep.*AM;
    end
    Sd.beep(iFreq,:) = beep;
end

Sd.tarFreqs = freqs;



% save sound information
eval(['save /Users/anette/Documents/MATLAB/Experiment2/Tones/' SubID '/SNDinfo Sd'])




% % timing
% % for timing test
% for iSession=1:1
%     for iTrial = 1:250;
%         maskSND = zeros(26,168000);
%         eval(['save /Users/anette/Documents/MATLAB/Experiment2/Tones/TimingTest/S' num2str(iSession) '/maskSND_T' num2str(iTrial) ' maskSND'])
%     end
% end