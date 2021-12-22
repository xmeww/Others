function [PossfreqEX IM SOAt] = IMsound_IM10(durMask, durBeep, fs, AMmask, meanSOA)
% function to create an informational sound mask, such as in the paper of
% Gutschalk et al., 2008 (PLOS Biology). !!!Important: Just as in the
% manual the mask is suited for a 1000 Hz target. (Since the surrounding
% frequencies are excluded. 
% durMask = duration of the overall mask sound (in seconds)
% durBeep = duration of one individual beep within the mask (in seconds)
% fs = sampling rate

% Gutschalk at al. used the following frequencies: [239 286 342 409 489 585
% 699 836 1000 1196 1430 1710 2045 2445 2924 3497 4181 5000]. Here, we will
% take the PossfreqEX frequencies. Hence, the 1000 Hz and 2 sideband
% frequencies have been excluded to minimize energetic masking
% PossfreqEX = [239 286 342 409 489 585 699 836 1000 1196 1430 1710 2045 2445 2924 3497 4181 5000];
generateFreq;
PossfreqEX = f;
nrMaskFreq = length(PossfreqEX);
% time axis of one beep
tTar = linspace(0,durBeep,fs*durBeep);

% calculate a cosine shaped ramp for each of the masks
durRamp = 0.01;
fcRamp  = 1/(durRamp*4);
tRamp   = linspace(0,durRamp,fs*durRamp);
RampON  = sin(2*pi*tRamp*fcRamp);
RampOFF = cos(2*pi*tRamp*fcRamp);
Ramp = [RampON ones(1,length(tTar)-2*length(tRamp)) RampOFF];

% calculate the masker sound
% SOA by which the masker beeps will be seperated
if meanSOA == 750
    SOAmasker = [0.1:0.01:0.65 0.85:0.01:1.4] - durBeep; % mean SOA of 750
elseif meanSOA == 1050
    SOAmasker = [0.55:0.01:0.85 1.25:0.01:1.55] - durBeep; % mean SOA = 1000 ms
end
tMask  = linspace(0,durMask,fs*durMask); % time axis
IM     = zeros(nrMaskFreq,length(tMask));
SOAt   = [];
nBeeps = 0;

for iMaskerFreq = 1:nrMaskFreq
    MaskSound = [];
    nSamp = length(MaskSound) - length(tMask);
    FreqIdx = iMaskerFreq;
    while nSamp < 0;      
        ERB = 24.7*(4.37*PossfreqEX(FreqIdx)/1000+1); % estimated equivavlent rectancular bandwidth
        fc  = RandSample((PossfreqEX(FreqIdx)-ERB):1:(PossfreqEX(FreqIdx)+ERB), [1,1]);
        MaskBeep = sin(2*pi*tTar*fc); % one beep
        % AM modulation of the sound. Only 20 Hz can really fit into a 250
        % ms interval, but 24 - 48 can approximately. (With the additional
        % ramping no clicks should be perceivable).
        fcAMTar  = RandSample([32 36 44 48], [1,1]); % 20 24 28 
        AMTar    = (1-cos(2*pi*tTar*fcAMTar))/2;
        if AMmask == 1
            MaskBeep = MaskBeep.*AMTar;
        end
        MaskBeep = MaskBeep.*Ramp; % ramp sound, to avoid clicks
        % SOA between sounds
        if isempty(MaskSound)
            SOA = zeros(1,round(fs*RandSample(0:0.01:02, [1,1])));
        else
            SOA = zeros(1,round(fs*RandSample(SOAmasker, [1,1])));
        end
        nSamp = length([MaskSound SOA MaskBeep ]) - length(tMask);
        if nSamp > 0
            break;
        end
        nBeeps = nBeeps + 1;
        SOAt.SOA(iMaskerFreq,nBeeps) = length(SOA)/fs;
        SOAt.frq = PossfreqEX(FreqIdx);
        MaskSound = [MaskSound SOA MaskBeep ];
    end
    nSamp = abs(length(MaskSound) - length(tMask));
    MaskSound = [MaskSound  zeros(1,nSamp)]; %#ok<*AGROW>
    IM(iMaskerFreq, :) =  MaskSound;
end


