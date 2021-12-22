function dB = AudiometryOldSystem(S)
% This function has been written to do an audiometry at the MEG Center in
% Tuebingen. Short 500 Hz tones will be presented to a subject a variable
% ISI.


%% ------------------------------------------------------------------------
% INPUT STRUCT
%
% if no input is given to the system, S will be created and all default
% settings will be used
if ~exist('S', 'var'),           S = [];               end


% beep length
if ~isfield(S, 'durStim'),       S.durStim = 4;        end

if ~isfield(S, 'fc'),            S.fc = 2000;          end

if ~isfield(S, 'attenuation'),   S.attenuation = 1/30; end

if ~isfield(S, 'dBdiff'),        S.dBdiff = 2;         end

if ~isfield(S, 'offset'),        S.offest = -50;       end



%%
% always needs to be put if using KbName!!!!!
KbName('UnifyKeyNames');
% define names of keys
quit   = KbName('q');
softer = KbName('s');
louder = KbName('l');



%%                  PREPARE SOUNDS
% =====================================================================
% Open sound driver for high timing precision, mono playback with 48
% kHZ sampling frequency:
InitializePsychSound(1);
freqBase = 48000;

% pa = audioport handle if the sound should be presented in stereo, one
% can simply change the channels from 1 to 2 and fill buffers for the
% different channels. Creating one audioport that plays the noise and
% one that plays the beeps.
%
% pahandle = PsychPortAudio('Open' [, deviceid][, mode][,
% reqlatencyclass][, freq][, channels][, buffersize][,
% suggestedLatency][, selectchannels][, specialFlags=0]);
paBeep    = PsychPortAudio('Open', [], [], 2, freqBase, 1);

% siehe auch BasicSoundScheduleDemo schaltet Treiber in Betriebsmode
% --> kuerzere soundonset latencies:
% In 'runMode' 1, the audio hardware and processing don't shut down at
% the end of audio playback. Instead, everything remains active in a
% ''hot standby'' state. This allows to very quickly (with low latency)
% restart sound playback via the 'RescheduleStart' function. The
% downside is a permanent use of system ressources even if no sound is
% playing. Future runMode settings may provide more interesting
% options, stay tuned...
%
% oldRunMode = PsychPortAudio('RunMode', pahandle [,runMode]);
PsychPortAudio('RunMode', paBeep, 1);


%----------------------------------------------------------------------
% SOUND MODULATION
fc   = S.fc; % carrier frequency
fam  = 40;
samplesBeep = freqBase * S.durStim;
time        = (0:samplesBeep-1) / freqBase; % time vector
beep        = sin(2 * pi * time * fc);
AM          = (1-cos(2 * pi * fam * time))/2;
beep        = beep.*AM;

% ramping the first 10ms of the sound
ramp = ones(1,length(time));
ramp(1:sum(time<0.01))         = linspace(0,1,sum(time<0.01));
ramp(end-sum(time<0.01)+1:end) = linspace(1,0,sum(time<0.01));
beep        = ramp.*beep;
mySoundBeep = beep*S.attenuation;


%----------------------------------------------------------------------
% PLAY SOUNDS
%
% this is important in order to not miss the very first deadline in the
% loop. The very first time the AudioPort is used, it needs some time.
%
% Load a sound that only contains zeros (i.e. no sound) into buffer
PsychPortAudio('FillBuffer', paBeep, zeros(1,100));

% startTime = PsychPortAudio('Start', pahandle [, repetitions=1] [,
% when=0] [, waitForStart=0] [, stopTime=inf] [, resume=0]);
PsychPortAudio('Start', paBeep, 0, 0, 1);
PsychPortAudio('Stop', paBeep);

ListenChar(2);
dB = S.offset;

% % start 40 dB less
% while dB > -70
%     mySoundBeep = mySoundBeep/(10^(S.dBdiff/20));
%     dB = dB - S.dBdiff;
%     KbWait([], 1);
% end
% fprintf('\n\n %i dB less - that are in total: %i dB \n\n', S.dBdiff, dB)

while 1
    % check if sound is not too loud
    if max(max(mySoundBeep)) > 0.4
        error('!!! check loudness !!! mask loudness: %f', max(max(mySoundBeep)))
    end
    
    % fill buffers with the actual sound
    PsychPortAudio('FillBuffer', paBeep, mySoundBeep);
    
    
    audiostatus = PsychPortAudio('GetStatus', paBeep);
    tDeadlineExact = (audiostatus.StartTime + RandSample(1:0.1:2));
    PsychPortAudio('Start', paBeep, 1, tDeadlineExact);
    
    %%                  START AUDIOMETRY
    while 1
        audiostatus = PsychPortAudio('GetStatus', paBeep);
        % only if the sound is not active
        if (audiostatus.Active == 0) && (audiostatus.RequestedStartTime <= audiostatus.StartTime)
            % give the new daedline for the sound, somewhen between
            % RandSample(1:0.5:3) + S.durStim
            tDeadlineExact = (audiostatus.StartTime + RandSample(1:0.5:3) + S.durStim);
            PsychPortAudio('Start', paBeep, 1, tDeadlineExact);
        end
        
        % -----------------------------------------------------------------
        % QUIT
        [~,~,keyCode] = KbCheck;
        % stopping the loop if user presses 'q'
        if keyCode(quit)
            break;
        elseif keyCode(softer) || keyCode(louder)
            break;
        end
    end
    
    
    if keyCode(softer)        
        mySoundBeep = mySoundBeep/(10^(S.dBdiff/20));
        dB = dB - S.dBdiff;
        KbWait([], 1);
        fprintf('\n %i dB less - that are in total: %i dB \n\n', S.dBdiff, dB)
    end
    
    
    if keyCode(louder)        
        mySoundBeep = mySoundBeep*(10^(S.dBdiff/20));
        dB = dB + S.dBdiff;
        KbWait([], 1);
        fprintf('\n\n %i dB more - that are in total: %i dB \n\n', S.dBdiff, dB)
    end
    
    PsychPortAudio('Stop', paBeep);
    
    % QUIT
    [~,~,keyCode] = KbCheck;
    % stopping the loop if user presses 'q'
    if keyCode(quit)
        break;
    end
end

ListenChar;


