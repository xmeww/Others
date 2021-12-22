function I = Exp2_IM_10(S)
% this function presents 2 beeps (AM modulated) within and informational
% sound mask, which is generated according to Gutschalk et al. The 2 beeps 
% seperated by 800 ms. According to the modality (A, V, AV paradigm) a
% flash my (additionally) occur. The flash occurs at the time of the second
% beep with a certain offset (see Debener et al.). After each trial the
% subject has to indicate whether he/she heard the tone.
%
% The function depends on:
% SoundAllFreqComb.mat (Generated with IMsound.m & GenerateMask.m)



% check how many default variables are used
defaults = 0;

% -------------------------------------------------------------------------
% INPUT STRUCT
%
% if no input is given to the system, S will be created and all default
% settings will be used
if ~exist('S', 'var'),           S = [];             end              


% -------------------------------------------------------------------------
% SUBJECT INFO
if ~isfield(S, 'SubID'),         S.SubID = 'Sub4'; defaults = defaults + 1;  end

sdir = '/Users/anette/Documents/MATLAB/Experiment2/Subdata/Test';
if ~isfield(S, 'saveDir'),       S.saveDir = sdir; defaults = defaults + 1;  end

if ~isfield(S, 'Session'),       S.Session = 1; defaults = defaults + 1;     end

if ~isfield(S, 'debug'),         S.debug = 0; defaults = defaults + 1;       end


% -------------------------------------------------------------------------
% AUDITORY
% percent aud signal strength
if ~isfield(S, 'SNRaud'),        S.SNRaud = 1/30; defaults = defaults + 1;      end

% percent aud noise strength
if ~isfield(S, 'SNRaudNoise'),   S.SNRaudNoise = 1/30; defaults = defaults + 1; end


% -------------------------------------------------------------------------
% VISUAL
% percent vis signal strength
if ~isfield(S, 'SNRvis'),        S.SNRvis = 0.07; defaults = defaults + 1;   end

% 1 = sine-wave grating, 2 = gabor grating
if ~isfield(S, 'visStim'),       S.visStim = 1; defaults = defaults + 1;     end


% -------------------------------------------------------------------------
% TRIAL & TARGETS
% nr of Tragets
if ~isfield(S, 'nrTrials'),     S.nrTrials = 240; defaults = defaults + 1;   end

% duration of one Aud stim in seconds
if ~isfield(S, 'durTarAud'),    S.durTarAud = 0.3; defaults = defaults + 1;  end

% duration of one Vis stim in seconds
if ~isfield(S, 'durTarVis'),    S.durTarVis = 0.2; defaults = defaults + 1;  end

% play masking sound or not
if ~isfield(S, 'mask'),         S.mask = 1; defaults = defaults + 1;         end

% play masking sound or not
if ~isfield(S, 'training'),     S.training = 0; defaults = defaults + 1;     end

% How many target frequencies
if ~isfield(S, 'nrFreqs'),      S.nrFreqs = 1:5; defaults = defaults + 1;    end

% SOA
if ~isfield(S, 'SOA'),          S.SOA = 1.05; defaults = defaults + 1;          end

% defines whether the visual stimulus wil be presented with the first
% (=1), second (=2) auditory target or both (= 3)
if ~isfield(S, 'timeAV'),       S.timeAV = 3; defaults = defaults + 1;       end


% -------------------------------------------------------------------------
% define loudness of sounds
if ~isfield(S, 'audScaling')
    dBdiff   = 2;
    audTimes = [30 30 30 30 30 30 30 30 30 30];
    for iTimes = 1: 10
       S.audScaling(iTimes) = (1/(10^(dBdiff/20)^audTimes(iTimes)));
    end    
    defaults = defaults + 1;
end

if ~isfield(S, 'audScalingTar')
    dBdiff   = 2;
    audTimes = [30 30 30 30 30];
    for iTimes = 1:5
       S.audScalingTar(iTimes) = (1/(10^(dBdiff/20)^audTimes(iTimes)));
    end
    defaults = defaults + 1;
end

if ~isfield(S, 'loudDB')
    dBdiff   = 2;
    audTimes = 25;
    S.loudDB = ((10^(dBdiff/20))^audTimes);
    defaults = defaults + 1;
end

if ~isfield(S, 'loudDBTar')
    dBdiff      = 2;
    audTimes    = 25;
    S.loudDBTar = ((10^(dBdiff/20))^audTimes);
    defaults = defaults + 1;
end


% -------------------------------------------------------------------------
% MEG
% MEG recording = 1
if ~isfield(S, 'MEG'),   S.MEG = 1; defaults = defaults + 1; end


% -------------------------------------------------------------------------
% load mask sound if necessary
load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/SNDinfo.mat']);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              START                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    
    %%                 START PASYCHTOOLBOX
    % =====================================================================
    % Make sure we're running on PTB-3
    AssertOpenGL;
    
    % always needs to be put if using KbName!!!!!
    KbName('UnifyKeyNames');   
    % define names of keys
    quit = KbName('q');
    
    
    % ---------------------------------------------------------------------
    % DEBUG
    % check for breakpoints
    brkp = dbstatus;
    if (S.debug == 1) || ~isempty(brkp)
        % simplyfies debugging by making windows transparent
        PsychDebugWindowConfiguration;
    else
        % only these keys will be checked
        % RestrictKeysForKbCheck([6 20]);
        % to prevent keyboard inputs during stimulstion inside the script
        % or matlab window. Important needs to be reset.
        ListenChar(2);
    end
     
    
    
    
    
    %%                 OPEN & PREPARE SCREENS
    % =====================================================================
    
    % Get the list of Screens and choose the one with the highest screen
    % number. Screen 0 is, by definition, the display with the menu bar.
    % Often when two monitors are connected the one without the menu bar is
    % used as the stimulus display.  Chosing the display with the highest
    % dislay number is a best guess about where you want the stimulus
    % displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    
    % Open double-buffered window: Optionally enable stereo output if
    % stereo > 0. 1 = Frame-sequential, 4 = Stereoscope (dual-display), 8 =
    % Anaglyph Red-Blue.
    PsychImaging('PrepareConfiguration');
    
    % Flip the image (only necessary if subject is in suspine position)
    % PsychImaging('AddTask', 'AllViews', 'FlipVertical');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'CheckOnly');
    
    
    backgroundvalue = 0.5;
    % Prepare setup of imaging pipeline for onscreen window. This is the
    % first step in the sequence of configuration steps
    [w, winRect]=PsychImaging('OpenWindow', screenNumber, 255*backgroundvalue); 
    Screen('ColorRange',w, 1);
    
    % Using this blendfunction the information which is presented to the
    % screen is non-transperent, i.e. info coming later will overwrite info
    % from before. If you want an additive screen function you can do the
    % following:
    % Screen(w,'BlendFunction',GL_ONE, GL_ONE);
    % BlendFunctions can rapidly be changed all over the script.
    Screen(w,'BlendFunction',GL_ONE, GL_ONE);
    % PsychColorCorrection('SetColorClampingRange', w, 0.01, 0.99)
    
    % define window height
    % win_w = (winRect(3)-winRect(1));
    win_h = (winRect(4)-winRect(2));
    
    % define text and layout of text that is used during test trials
    Screen('TextFont', w, 'Arial');
    Screen('TextSize', w, 23);
    Screen('TextStyle', w, 1);
    

    
    
    
    %%                MEASURING THE REFRESH RATE
    % =====================================================================
    
    % Finding out the RefreshRate of the monitor Query nominal framerate as
    % returned by Operating system: If OS returns 0, then we assume that we
    % run on a flat-panel with fixed 60 Hz refresh interval.
    framerate=Screen('NominalFramerate', w);
    if (framerate==0)
        framerate=60;
    end;
    
    % gives ms of one time frame eg. 1/60 = 16.6667 ms
    ifinominal=1 / framerate;
    fprintf('\n \n The refresh interval reported by the operating system is %2.5f ms.\n', ifinominal*1000);
    fprintf('\n \n The measured framerate is %2.5f ms.\n', framerate);
    
    % Perform a calibration loop to determine the "real" interframe
    % interval for the given gfx-card + monitor combination:
    %
    %
    % Measure monitor refresh interval again, just for fun... This will
    % trigger a calibration loop of minimum 100 valid samples and return
    % the estimated ifi in 'ifi': We require an accuracy of 0.05 ms ==
    % 0.00005 secs. If this level of accuracy can't be reached, we time out
    % after 20 seconds...
    [ifi nvalid stddev ]= Screen('GetFlipInterval', w, 100, 0.00005, 5);
    fprintf('\n \n Measured refresh interval, as reported by "GetFlipInterval" is %2.5f ms. (nsamples = %i, stddev = %2.5f ms)\n \n', ifi*1000, nvalid, stddev*1000);
    
   
    
    
    
        
    %%             BUTTON RESPONSES
    % Open bitwhacker box:
    % defines values of Button Presses. 1 - normal, 0 - inverted
    if S.MEG == 1
        % MEG system:
        bitwhacker = CMUBox('Open', 'bitwhacker', FindSerialPort('usbmodem', 1), 'norelease', [], [0 0 1 1 1 1 1 1 1]); 
        BuYes = 1; 
        BuNo  = 2;
    elseif S.MEG == 0
        % Test mit buttons:
        bitwhacker = CMUBox('Open', 'bitwhacker', FindSerialPort('usbmodem', 1), 'norelease', [], [1 1 1 1 1 1 1 1 1]); 
        BuYes = 9; % button ganz aussen in der Ecke
        BuNo  = 8; % links neben Button 9 
    end
    
    % Clear button queue:
    while ~isempty(CMUBox('GetEvent', bitwhacker, 0)); end;
    
    
    
    
    
    %%                      TARGETS
    % =====================================================================
    % number of targets. Note that AV & A conditions occur more often since
    % they will be split into detected and undetected. For training all
    % trials will occur evenly often.
    nTrials = S.nrTrials; 
    if (S.training == 0) && (S.mask == 1)
        nA = 3;
        nV = 1;
    elseif (S.training == 1) || (S.mask == 0)
        nA = 2;
        nV = 2;
    end
    nAtar  = nA/8*nTrials; % 3/8*240 =90
    nAVtar = nA/8*nTrials; % 90
    nVtar  = nV/8*nTrials; % 1/8*240 =30
    nNOtar = nV/8*nTrials; % 30  
    
    % Important the jitter MUST be a multiple of the refresh rate!!!!!!!!!!
    % minimum SOA befor the first target stimulus may occur.
    minSOA  = ifi;
    % maximum SOA until the first target stimulus has to occur.
    maxSOA  = 0.5;
    % mean SOA
    meanSOA = 0.25;
    % there may be as many possible SOAs as there ar V/noStim trials.
    % However the flip may occur only for SOAs which are a multiple of the
    % SOA. If there are less possible SOAs than targets, we fill the rest
    % with the mean SOA.
    % create a vector of mean SOA values that is as long as the shortest
    % number of targets (which are V or noStim targets)
    possSOAs = ones(1,nVtar)*meanSOA; % [1 30] 0.25
    % create a vector of all possible SOAs that are multiples of the flip.
    % Fill the rest with mean SOAs
    SOA = minSOA:ifi:maxSOA;
    possSOAs(1:length(SOA)) = SOA; 
    
    % if the SOA vector is longer than the number of targets, cut the rest.
    possSOAs = possSOAs(randperm(length(possSOAs)));
    if length(possSOAs)>nVtar
        possSOAs(nVtar+1:end) = [];
    end
    
    % create a vector for A and AV targets that might need to be longer,
    % since there are more A and V stimuli.
    possSOAs2 = possSOAs;
    while length(possSOAs2) < nAtar
        possSOAs2 = [possSOAs2 possSOAs];
    end       
    
    
    % Target modality vectors
    %               noStim            AV targets            V targets                        A targets                     
    Targets      = [zeros(2,nNOtar),  ones(2,nAVtar),      [zeros(1,nVtar); ones(1,nVtar)], [ones(1,nAtar); zeros(1,nAtar)]];
    Targets(3,:) = [possSOAs,         possSOAs2,            possSOAs,                        possSOAs2];
    Targets(4,:) = [ones(1,nNOtar)*4, ones(1,nAVtar),       ones(1,nVtar)*3,                 ones(1,nAtar)*2]; 
    
    % randomize onsets
    Targets = Targets(:,randperm(nTrials));
    
    % add one AV trial, since the timing of the first target may be bad
    Targets = [[1; 1; 0.25; 1], Targets];
    nTrials = nTrials + 1;
    
    % get out events
    AudTarget   = Targets(1,:); % 0 = no A target, 1 = A target
    VisTarget   = Targets(2,:); % 0 = no V target, 1 = V target
    firstTarSOA = Targets(3,:); % onset of the first target sound
    Conditions  = Targets(4,:); % 1 = AV, 2 = A, 3 = V, 4 = NoStim;
    
    

    
    
    %%                  PREPARE SOUNDS
    % =====================================================================
    % Open sound driver for high timing precision, mono playback with 48
    % kHZ sampling frequency:
    freqBase = Sd.freqBase;
    freqs    = Sd.tarFreqs;
    InitializePsychSound(1);
    
    nMaskFreq = length(Sd.PossFreq);
    for iMaskFreq = 1:nMaskFreq
        % pa = audioport handle if the sound should be presented in stereo, one
        % can simply change the channels from 1 to 2 and fill buffers for the
        % different channels. Creating one audioport that plays the noise and
        % one that plays the beeps.
        %
        % pahandle = PsychPortAudio('Open' [, deviceid][, mode][,
        % reqlatencyclass][, freq][, channels][, buffersize][,
        % suggestedLatency][, selectchannels][, specialFlags=0]);
        eval(['paMask' num2str(iMaskFreq) ' = PsychPortAudio(''Open'', [], [], 2, freqBase, 1);']);
           
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
        eval(['PsychPortAudio(''RunMode'', paMask' num2str(iMaskFreq) ', 1);']);
    end
    paBeep = PsychPortAudio('Open', [], [], 2, freqBase, 1);
    PsychPortAudio('RunMode', paBeep, 1);

    
    
    
    
    %%                 CREATING VISUAL STIMULUS
    % =====================================================================
    % generate visual stimulus using gaborPatch.mat
    [grating gabor] = gaborPatch(75, 50, 40); 
    rotationAngle   = 0;
    
    % sine wave grating
    if S.visStim == 1
        visStim = grating;
    % square wave grating
    elseif S.visStim == 2
        visStim = gabor;
    % use white square    
    elseif S.visStim == 3
        visStim = ones(size(grating));
    end   
    
    % make Texture
    gaborTex = Screen('MakeTexture', w, visStim, [], [], 2);
    
   
    
   
 
    %%                INITIALIZE VARIABLES
    % =====================================================================
    % onset of the trial
    tStart   = zeros(1,nTrials);
    wakeup   = zeros(6,nTrials);
    % flips at targets
    timeFlip = zeros(2,nTrials);
    % aud target sounds
    audTarTime = zeros(2,nTrials);
    
    % onset of the trigger
    triggerOnset = zeros(2,nTrials);
    
    BuTooEarly = zeros(1,nTrials);
    
    % buttonResponse
    BuResp = zeros(1,nTrials);
    hit    = zeros(1,nTrials);
    tResp  = zeros(1,nTrials);
    
    % find locations of Target frequencies
    for iT = S.nrFreqs; 
        fr(iT) = find(Sd.PossFreq == Sd.tarFreqs(iT)); 
    end
    
    
    % audiovisual offset, i.e. the time the auditory stimulus is delayed
    % with respect to the visual stimulus
    AVoffset = 0.05;
    
    
    
    
    
    %% SOUND 1st trial
    %----------------------------------------------------------------------
    % offset of the sound in seconds that accounts for any audiovisual
    % delays.
    if S.MEG == 1
        offsetAudio = 0.012;
    elseif S.MEG == 0
        offsetAudio = 0;
        fprintf('\n \n !!!WARNING!!!! \n \n OFFSET = 0!!!! \n \n')
    end
    
    % scale the auditory target stimulus
    for iTimes = 1:5
        Sd.beep(iTimes,:) = Sd.beep(iTimes,:) * S.audScalingTar(iTimes) * S.loudDBTar;
    end
    
    if S.mask == 1
        % load mask for the first trial
        load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/S' num2str(S.Session) '/maskSND_T1.mat']);
        
        
        % scale mask according to hearing threshold of the subject
        maskSND = maskSND * S.SNRaudNoise; %#ok<*NODEF>
        
        maskSND(1:5,:)    = maskSND(1:5,:)   * S.audScaling(1);
        maskSND(6:11,:)   = maskSND(6:11,:)  * S.audScaling(2);
        maskSND(12:14,:)  = maskSND(12:14,:) * S.audScaling(3);
        maskSND(15:16,:)  = maskSND(15:16,:) * S.audScaling(4);
        maskSND(17,:)     = maskSND(17,:)    * S.audScaling(5);
        maskSND(18,:)     = maskSND(18,:)    * S.audScaling(6);
        maskSND(19:20,:)  = maskSND(19:20,:) * S.audScaling(7);
        maskSND(21:22,:)  = maskSND(21:22,:) * S.audScaling(8);
        maskSND(23,:)     = maskSND(23,:)    * S.audScaling(9);
        maskSND(24:end,:) = maskSND(24:end,:)* S.audScaling(10);
        
        maskSND = maskSND * S.loudDB;
        % check if mask is not too loud
        if sum(max(maskSND,[],2)) >= 1.1 
            error('!!! check loudness !!! mask loudness: %f', max(max(maskSND)))
        end
    end
    tarFreqIdx = RandSample(S.nrFreqs, [1,1]);
    tarFreq    = fr(tarFreqIdx);
    tarFreqTrial(1,1) = freqs(tarFreqIdx);
    
    % creat masking sounds
    % check which masks will be used in the next trial
    maskOnNext = ones(1,nMaskFreq);
    maskOnNext(tarFreq-3:tarFreq+3) = 0;
    
    if (S.mask == 1)
        for iMaskFreq = 1:nMaskFreq
            % stop mask sounds
            eval(['PsychPortAudio(''Stop'', paMask' num2str(iMaskFreq) ');']);
            % refill the buffer
            eval(['PsychPortAudio(''FillBuffer'', paMask' num2str(iMaskFreq) ', maskSND(iMaskFreq,:));']);
            % start only relevant streams
            if (S.mask == 1) && (maskOnNext(iMaskFreq) == 1)
                eval(['PsychPortAudio(''Start'', paMask' num2str(iMaskFreq) ', 0);']);
            end
        end
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                          START                                    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iTrial = 1:nTrials
        %% --------------------- SOUND ------------------------------------

        beep(1,:) = Sd.beep(tarFreqIdx,:);
        % cosine shaped ramp (10 ms on and offset)
        samplesBeep = freqBase * S.durTarAud;
        time        = (0:samplesBeep-1) / freqBase; % time vector
        durRamp = 0.01;
        fcRamp  = 1 /(durRamp*4);
        tRamp   = linspace(0,durRamp,freqBase*durRamp);
        RampON  = sin(2*pi*tRamp*fcRamp);
        RampOFF = cos(2*pi*tRamp*fcRamp);
        ramp    = [RampON ones(1,length(time)-2*length(tRamp)) RampOFF];
        beep    = ramp .* beep;
        
        % IMPORTANT The two streams will be added. While being added they
        % are not allowed to exceed values between -1 and 1. Otherwise they
        % will be clipped.
        
        % sound vector containing the modulated tones. Vector beep ranges
        % from 0 to 1. Multiplying by scaling factor will make sure that
        % beep + noise does not exceed 1.
        mySoundBeep = beep * S.SNRaud;
        % check if sound is not too loud
        if max(max(mySoundBeep)) > 0.2
            error('!!! check loudness !!! mask loudness: %f', max(max(mySoundBeep)))
        end
        
        % fill buffers with the target sound
        PsychPortAudio('FillBuffer', paBeep, mySoundBeep);
        
        % get the start time of the mask
        for iMaskFreq = 1:nMaskFreq
            eval(['audioMaskStatus = PsychPortAudio(''GetStatus'', paMask' num2str(iMaskFreq) ');'])
            if audioMaskStatus.Active == 1
                audMaskTime(iMaskFreq, iTrial) = audioMaskStatus.StartTime;        
            end
        end
   
        
        
        
        
        %%                  SYNCHRONIZING
        % =================================================================
        
        % Perform some initial Flip to get us in sync with retrace: tvbl is
        % the timestamp (system time in seconds) when the retrace started.
        % We need it as a reference value for our display timing:
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        tvbl         = Screen('Flip', w, GetSecs + 0.5*ifi);
        tDeadlineVis = tvbl + ifi - 0.5 * ifi;
        
        
        
        
        
        %%                     1st INTERVAL   
        % =================================================================
        %
        % ---------------------- 750 ms -----------------------------------
        % 1. FIXATION (before trial starts)
        % Draw fixation cross in one color
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('DrawingFinished', w, 0);
        tStart(1,iTrial) = Screen('Flip', w, tDeadlineVis);
        
        % load mask as we have 750 ms to do so, we will not have any
        % further delays 
        if S.mask ==1
            load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/S' num2str(S.Session) '/maskSND_T' num2str(iTrial + 1) '.mat']);
        end
        
        % ---------------- firstTarSOA ------------------------------------
        % 1. FIXATION (for firstTarSOA seconds)
        % Draw fixation cross in one color
        tDeadlineVis = tStart(1,iTrial) + 0.75 - 0.5 * ifi;
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('DrawingFinished', w, 0);
        wakeup(1, iTrial) = Screen('Flip', w, tDeadlineVis);
        
        
        
        % ------------------ S.durTarVis ms -------------------------------
        % 2. Targets (for S.durTar Seconds)
        % calculate deadline for target presentation
        tDeadlineVis = wakeup(1, iTrial) + firstTarSOA(iTrial) - 0.5 * ifi;
        
        % Visual Targets
        if (VisTarget(iTrial) == 1) && ((S.timeAV == 1) || (S.timeAV == 3))
            Screen(w,'BlendFunction',GL_ONE, GL_ONE);
            Screen('DrawTexture', w, gaborTex, [], [], rotationAngle, [], [], S.SNRvis);
            Screen(w,'BlendFunction',GL_ONE, GL_ZERO);
        end
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('DrawingFinished', w, 0);
        
        % Flip
        timeFlip(1, iTrial) = Screen('Flip', w, tDeadlineVis); 
        
        % Send Trigger
        try
            % USBTTL('Open'); should be done in the script to start this
            % up, since it needs some time.
            USBTTL('Set',1);
        catch
            % fprintf('\n \n !!! WARNING !!!  No triggers have been sent! \n \n');
        end
        triggerOnset(1, iTrial) = GetSecs;
        
        % Auditory Targets
        if AudTarget(iTrial) == 1
            PsychPortAudio('Start', paBeep, 1, triggerOnset(1, iTrial) + offsetAudio + AVoffset);
        end
        
        % save flip time
        wakeup(2, iTrial) = timeFlip(1, iTrial);
              
        
        
        % ---------------------- SOA --------------------------------------
        % 3. FIXATION (time between 2 targets)
        % calculate deadline for presentation
        
        tDeadlineVis = wakeup(2,iTrial) + S.durTarVis - 0.5 * ifi;
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('DrawingFinished', w, 0);
        wakeup(3,iTrial) = Screen('Flip', w, tDeadlineVis); 
        % Close Trigger
        try
            USBTTL('Set',0);
        catch
        end
        
        % check exact timing of the auditory stimulus
        if AudTarget(iTrial) == 1
            audiostatus = PsychPortAudio('GetStatus', paBeep);
            audTarTime(1, iTrial) = audiostatus.StartTime;
        end
        
        
        % --------------------- 250 ms ------------------------------------
        % 4. Second Target
        tDeadlineVis = wakeup(3, iTrial) + (S.SOA - S.durTarVis) - 0.5 * ifi;
        
        % Visual Targets
        if (VisTarget(iTrial) == 1) && ((S.timeAV == 2) || (S.timeAV == 3))
            Screen(w,'BlendFunction',GL_ONE, GL_ONE);
            Screen('DrawTexture', w, gaborTex, [], [], rotationAngle, [], [], S.SNRvis);
            Screen(w,'BlendFunction',GL_ONE, GL_ZERO);
        end
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('DrawingFinished', w, 0);
        
        % Flip
        timeFlip(2, iTrial) = Screen('Flip', w, tDeadlineVis); 
        
        % Send Trigger
        try
            % USBTTL('Open'); should be done in the script to start this
            % up, since it needs some time.
            USBTTL('Set',1);
        catch
            % fprintf('\n \n !!! WARNING !!!  No triggers have been sent! \n \n');
        end
        triggerOnset(2, iTrial) = GetSecs;
        wakeup(4, iTrial)       = timeFlip(2, iTrial);
        
        % Auditory Targets
        if AudTarget(iTrial) == 1
            PsychPortAudio('Start', paBeep, 1, triggerOnset(2, iTrial) + offsetAudio + AVoffset);
        end
        
        
        
        % -------------------- END ----------------------------------------
        % 4. Fixation
                
        tDeadlineVis = wakeup(4, iTrial) + S.durTarVis - 0.5 * ifi;
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        wakeup(5,iTrial) = Screen('Flip', w, tDeadlineVis); 
        
        % Close Trigger
        try
            USBTTL('Set',0);
        catch
        end
        
        % check exact timing of the auditory stimulus
        if AudTarget(iTrial) == 1
            audiostatus = PsychPortAudio('GetStatus', paBeep);
            audTarTime(2, iTrial) = audiostatus.StartTime; 
        end
        
        
        
        % scale mask as we have 750 ms to do so, we will not have any
        % further delays 
        if S.mask == 1
            % scale mask according to hearing threshold of the subject
            maskSND = maskSND*S.SNRaudNoise;
            
            maskSND(1:5,:)    = maskSND(1:5,:)   * S.audScaling(1);
            maskSND(6:11,:)   = maskSND(6:11,:)  * S.audScaling(2);
            maskSND(12:14,:)  = maskSND(12:14,:) * S.audScaling(3);
            maskSND(15:16,:)  = maskSND(15:16,:) * S.audScaling(4);
            maskSND(17,:)     = maskSND(17,:)    * S.audScaling(5);
            maskSND(18,:)     = maskSND(18,:)    * S.audScaling(6);
            maskSND(19:20,:)  = maskSND(19:20,:) * S.audScaling(7);
            maskSND(21:22,:)  = maskSND(21:22,:) * S.audScaling(8);
            maskSND(23,:)     = maskSND(23,:)    * S.audScaling(9);
            maskSND(24:end,:) = maskSND(24:end,:)* S.audScaling(10);
            
            maskSND = maskSND*S.loudDB;
            % check if mask is not too loud
            if sum(max(maskSND,[],2)) >= 1.1
                error('!!! check loudness !!! mask loudness: %f', max(max(maskSND)))
            end
        end
        
        
        
        
        %%                     RESPONSE
        % =================================================================       
        % if you are sure about early button responses present questions
        % and wait for the actual or wanted button response
        tDeadlineVis = wakeup(5, iTrial) + (0.75 + AVoffset + (S.durTarAud - S.durTarVis)) - 0.5 * ifi; 
        % print question to be asked
        DrawFormattedText( w, 'Gab es einen Zielton?', 'center',  win_h/2 - 60, 0.4);
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        tvbl = Screen('Flip', w, tDeadlineVis); 
        wakeup(6,iTrial) = tvbl;
        
        
        % empty CMU box and check if the subject has answered before he/she
        % should have answered.
        % A problem still is that you cannot distinguish between too early
        % (within this trial) and too late or double responses within the
        % last trial.
        Early = 1;
        % is CMU box empty
        evtEarly  = CMUBox('GetEvent', bitwhacker, 0);
        while ~isempty(evtEarly); % if not empty (i.e. early button press)
            respEarly(Early, iTrial) = evtEarly;            
            BuTooEarly(1, iTrial)    = 1;
            % check if there was another button press
            evtEarly = CMUBox('GetEvent', bitwhacker, 0);
            Early    = Early + 1;
        end
        
        
        %% meanwhile start the new mask to not loose too much time
        % -----------------------------------------------------------------        
        tarFreqIdx = RandSample(S.nrFreqs, [1,1]);
        tarFreq    = fr(tarFreqIdx);
        tarFreqTrial(1,iTrial+1) = freqs(tarFreqIdx);
        
        % creat masking sounds
        % check which masks will be used in the next trial
        maskOnNext = ones(1,nMaskFreq);
        maskOnNext(tarFreq-3:tarFreq+3) = 0;
        
        if (S.mask == 1)
            for iMaskFreq = 1:nMaskFreq
                % stop mask sounds
                eval(['PsychPortAudio(''Stop'', paMask' num2str(iMaskFreq) ');']);
                % refill the buffer
                eval(['PsychPortAudio(''FillBuffer'', paMask' num2str(iMaskFreq) ', maskSND(iMaskFreq,:));']);
                % start only relevant streams
                if (S.mask == 1) && (maskOnNext(iMaskFreq) == 1)
                    eval(['PsychPortAudio(''Start'', paMask' num2str(iMaskFreq) ', 0);']);
                end
            end
        end
        
        
        
        %% participants get maximally 2 seconds to respond
        % -----------------------------------------------------------------
        maxResponseTime = 2;
        evt  = [];
        while isempty(evt)
            evt   = CMUBox('GetEvent', bitwhacker, 0);
            tResp(1,iTrial) = GetSecs - tvbl;
            if tResp(1,iTrial) > maxResponseTime
                break;
            end
            
        end
        
        % QUIT
        [~,~,keyCode] = KbCheck;
        % stopping the loop if user presses 'q'
        if keyCode(quit)
            break;
        end
        
        % save response
        if ~isempty(evt)
            resp(iTrial) = evt; 
            if resp(iTrial).state == BuYes
                BuResp(1,iTrial) = 1;
                if AudTarget(1,iTrial) == 1
                    hit(1,iTrial) = 1; % Hit
                else
                    hit(1,iTrial) = 2; % FA
                end
            elseif resp(iTrial).state == BuNo
                BuResp(1,iTrial) = -1;
                if AudTarget(1,iTrial) == keyCode
                    hit(1,iTrial) = -2; % Miss 
                else
                    hit(1,iTrial) = -1; % correct rejection
                end
            end
        else
            resp(iTrial).time       = []; %#ok<*AGROW>
            resp(iTrial).streamTime = [];
            resp(iTrial).state      = [];
            resp(iTrial).trouble    = [];
            resp(iTrial).deltaScan  = [];
            BuResp(1,iTrial)  = 0; % no response
            hit(1,iTrial)     = 0;
        end
        
    if S.training == 1
        if abs(hit(1,iTrial)) == 1
            txt = 'RICHTIG';
        elseif abs(hit(1,iTrial)) == 2
            txt = 'FALSCH';
        else
            txt = 'zu langsam!';
        end
        DrawFormattedText( w, txt, 'center',  win_h/2 - 80, 1);
        DrawFormattedText( w, '+', 'center',  'center', [0.2 0.5 0.2]);
        Screen('Flip', w, GetSecs + 0.5*ifi);
        WaitSecs(1);
    end
        
       
        
    end % Draw next frame
    
    % print results:
    fprintf('\n \n RESULTS \n ');
    fprintf('\n nr of trials: %i \n nr of hits: %i \n nr of correct rejections: %i \n', iTrial, sum(hit==1), sum(hit == -1));
    fprintf('\n nr of misses: %i \n nr of FA: %i \n too late: %i \n', sum(hit==-2), sum(hit==2), sum(hit == 0));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          END LOOP                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%                         SAVE
    % =====================================================================
    
    I.Subinfo.subID   = S.SubID;
    I.Subinfo.session = S.Session;
    
    I.targets.possSOA       = possSOAs;
    I.targets.minSOA        = minSOA;
    I.targets.maxSOA        = maxSOA;
    I.targets.targets       = Targets;
    I.targets.audTargets    = AudTarget;
    I.targets.visTargets    = VisTarget;
    I.targets.SOA           = firstTarSOA;
    I.targets.tarFreqs      = tarFreqTrial;
      
    I.conditions = Conditions;
    
    I.visStim = visStim;
    
    I.timing.AVoffset      = AVoffset;
    I.timing.offsetAudio   = offsetAudio;
    I.timing.wakeup        = wakeup;
    I.timing.startTrial    = tStart;
    I.timing.flipTarOnset  = timeFlip; 
    I.timing.audTarOnset   = audTarTime;
    if S.mask == 1
        I.timing.audMaskTime   = audMaskTime;
    end
    I.timing.triggerOnsets = triggerOnset;
    
    I.sound.frequencies = Sd.PossFreq;
    I.sound.targetFreqs = Sd.tarFreqs;
    
    I.response.BW             = resp;
    I.response.button         = BuResp;
    I.response.hit            = hit;
    I.response.time           = tResp;
    I.response.tooEarly.trial = BuTooEarly;
    if exist('respEarly', 'var')
        I.response.tooEarly.evt   = respEarly;
    end
    
    
    savedfname=[S.saveDir, '/Subject_', S.SubID, '_', num2str(S.Session), '_Training', num2str(S.training), '.mat'];
    save(savedfname, 'I')
    
    
    
    
    
    %%                         CLOSE
    % =====================================================================
    % ---------------------------------------------------------------------
    % Keypress
    % makeing keypresses being 'heard' again
    ListenChar;
    
    % ---------------------------------------------------------------------
    % AUDIO
    % Stop soundplayback:
    PsychPortAudio('Stop', paBeep);
    PsychPortAudio('Stop', paMask1);    
    PsychPortAudio('Stop', paMask2);    
    % Shutdown sound driver:
    PsychPortAudio('Close');
    
    % ---------------------------------------------------------------------
    % SCREENS
    %
    % Close display: If we skipped/missed any presentation deadline during
    % Flip, Psychtoolbox will automatically display some warning message on
    % the Matlab console:
    Screen('CloseAll');
    
    % ---------------------------------------------------------------------
    % Shutdown realtime scheduling:
    Priority(0);
    
  
    % Close bitwhacker:
    CMUBox('Close', bitwhacker);
  
    % check number of default variables used
    fprintf('\n \n Number of defaults used: %i \n \n', defaults)
    
    
     
    
    
    
catch  %#ok<*CTCH>
    % This "catch" section executes in case of an error in the "try"
    % section above. Importantly, it closes the onscreen window if its open
    % and shuts down realtime-scheduling of Matlab
    
    % ---------------------------------------------------------------------
    % Keypress
    % makeing keypresses being 'heard' again
    ListenChar;
    
    try %#ok<TRYNC>
        % ---------------------------------------------------------------------
        % Close bitwhacker:
        CMUBox('Close', bitwhacker);
        
        % ---------------------------------------------------------------------
        % AUDIO
        % Stop soundplayback:
        PsychPortAudio('Stop', paBeep);
        PsychPortAudio('Stop', paMask1);
        PsychPortAudio('Stop', paMask2);
        % Shutdown sound driver:
        PsychPortAudio('Close');
    end
    
    % ---------------------------------------------------------------------
    % SCREENS
    %
    % Close display: If we skipped/missed any presentation deadline during
    % Flip, Psychtoolbox will automatically display some warning message on
    % the Matlab console:
    Screen('CloseAll');    
    
    % ---------------------------------------------------------------------
    % Shutdown realtime scheduling:
    Priority(0);  
    
    % ---------------------------------------------------------------------
    % output error that occured
    psychrethrow(psychlasterror); 
end

% Done!!!! :)