function dB = AudiometryComparison(S)
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
if ~isfield(S, 'durStim'),       S.durStim = 0.3;      end

if ~isfield(S, 'fc'),            S.fc = 1122;          end

if ~isfield(S, 'attenuation'),   S.attenuation = 1/30; end

if ~isfield(S, 'dBdiff'),        S.dBdiff = 2;         end

if ~isfield(S, 'offset'),        S.offset = 0;         end

if ~isfield(S, 'MEG'),           S.MEG = 1;            end



%%
% always needs to be put if using KbName!!!!!
KbName('UnifyKeyNames');
% define names of keys
quit     = KbName('q');
refSound = KbName('r');
softer   = KbName('s');
louder   = KbName('l');



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
paBeepRef = PsychPortAudio('Open', [], [], 2, freqBase, 1);

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
PsychPortAudio('RunMode', paBeepRef, 1);


%----------------------------------------------------------------------
% SOUND MODULATION
fc   = S.fc; % carrier frequency
fcRef= 2000;
fam  = 40;
samplesBeep = freqBase * S.durStim;
time        = (0:samplesBeep-1) / freqBase; % time vector
beep        = sin(2 * pi * time * fc);
beepRef     = sin(2 * pi * time * fcRef);
AM          = (1-cos(2 * pi * fam * time))/2;
beep        = beep.*AM;
beepRef     = beepRef.*AM;

% ramping the first 10ms of the sound
ramp = ones(1,length(time));
ramp(1:sum(time<0.01))         = linspace(0,1,sum(time<0.01));
ramp(end-sum(time<0.01)+1:end) = linspace(1,0,sum(time<0.01));
beep           = ramp.*beep;
beepRef        = ramp.*beepRef;
mySoundBeep    = beep*S.attenuation;
mySoundBeepRef = beepRef*S.attenuation;


%----------------------------------------------------------------------
%             BUTTON RESPONSES
% Open bitwhacker box:
% defines values of Button Presses. 1 - normal, 0 - inverted
if S.MEG == 1
    % MEG system:
    bitwhacker = CMUBox('Open', 'bitwhacker', FindSerialPort('usbserial-00001004', 1), 'norelease', [], [0 0 1 1 1 1 1 1 1]);
    BuYes = 1;
    BuNo  = 2;
elseif S.MEG == 0
    % Test mit buttons:
    bitwhacker = CMUBox('Open', 'bitwhacker', FindSerialPort('usbserial-00001004', 1), 'norelease', [], [1 1 1 1 1 1 1 1 1]);
    BuYes = 9; % button ganz aussen in der Ecke
    BuNo  = 8; % links neben Button 9
end

% Clear button queue:
while ~isempty(CMUBox('GetEvent', bitwhacker, 0)); end;




%----------------------------------------------------------------------
% PLAY SOUNDS

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
    PsychPortAudio('FillBuffer', paBeepRef, mySoundBeepRef);
    
    % -----------------------------------------------------------------
    % check for responses
    Bu = 0;
    evt = CMUBox('GetEvent', bitwhacker, 0);
    if ~isempty(evt)
        Bu = evt.state;
    end
    
    [~,~,keyCode] = KbCheck;
    % stopping the loop if user presses 'q'
    if keyCode(quit)
        break;
    end
   
    
    if keyCode(refSound) 
        PsychPortAudio('Start', paBeepRef);
        PsychPortAudio('Start', paBeep, 1, GetSecs + 0.5);
    end
    
    
    if keyCode(softer) || (Bu == BuYes)       
        mySoundBeep = mySoundBeep/(10^(S.dBdiff/20));
        dB = dB - S.dBdiff;
        KbWait([], 1);
        fprintf('\n %i dB less - that are in total: %i dB \n\n', S.dBdiff, dB)
    end
    
    
    if keyCode(louder) || (Bu == BuNo)       
        mySoundBeep = mySoundBeep*(10^(S.dBdiff/20));
        dB = dB + S.dBdiff;
        KbWait([], 1);
        fprintf('\n\n %i dB more - that are in total: %i dB \n\n', S.dBdiff, dB)
    end
    
    
    % QUIT
    [~,~,keyCode] = KbCheck;
    % stopping the loop if user presses 'q'
    if keyCode(quit)
        break;
    end
end

CMUBox('Close', bitwhacker);
ListenChar;


