function info = DEC_2_Present(S) 
%
%

%%                         PARSING INPUT PARAMETERS
% =========================================================================
p = inputParser;

fieldsReqS = {'audioParam','auditoryLocationLevels','auditoryStimuli',...
    'condition','diameterSphere','expMode','expName','expStage','feedback',...
    'hand','iDay','iSessionInDay','iSessionOverall','locationAuditory',...
    'locationVisual','maxVisAngle','nTrials','reliabilityVisual',...
    'runScriptName','screenRes','sessionOrder','setupSpec','subID',...
    'syncMode','task','testMode','timing','user','viewDist',...
    'visualLocationLevels','visualStimuli'};
fieldsReqTiming = {'t_endResponse','t_feedbackDuration','t_respCueDuration',...
    't_stimDuration','t_toResponse','t_toStimOnset'};
fieldsReqAudioPara m = {'details','mode','nChannels','sr'};
fieldsReqSetupSpec = {'audioDeviceID','corrGamma','ID','info','latBias','room'};

checkS = @(x) all(isfield(x,fieldsReqS)) && ...
              all(isfield(x.timing,fieldsReqTiming)) && ...
              all(isfield(x.audioParam,fieldsReqAudioParam)) && ...
              all(isfield(x.setupSpec,fieldsReqSetupSpec));

addRequired(p,'S',checkS);

parse(p,S);

S = p.Results.S;

%%                                 START
%%=========================================================================
presentFunctName = mfilename;

try
    % Make sure we're running on PTB-3
    AssertOpenGL;
    
    % Always needs to be put if using KbName! 
    KbName('UnifyKeyNames');
    
    % to exit the program
    quitKey   = KbName('ESCAPE');
    returnKey = KbName('RETURN');
    % These keys are necessary for controlling the eyelink form the
    % presentation computer
    aKey      = KbName('a');
    cKey      = KbName('c');
    vKey      = KbName('v');
    leftArrow = KbName('LeftArrow');
    rightArrow = KbName('RightArrow');
    % defining keys for input 
    % right hand keys
    loc1KeyRight = KbName('4');    % corresponding to far left
    loc2KeyRight = KbName('5');    % corresponding to close left
    loc3KeyRight = KbName('6');    % corresponding to close right
    loc4KeyRight = KbName('+');    % corresponding to far right
    % left hand keys
    loc1KeyLeft = KbName('7');    % corresponding to far left
    loc2KeyLeft = KbName('8');    % corresponding to close left
    loc3KeyLeft = KbName('9');    % corresponding to close right
    loc4KeyLeft = KbName('-');    % corresponding to far right
    
    respKeysRight = [loc1KeyRight,loc2KeyRight,loc3KeyRight,loc4KeyRight];
    respKeysLeft = [loc1KeyLeft,loc2KeyLeft,loc3KeyLeft,loc4KeyLeft];
    RestrictKeysForKbCheck([quitKey,returnKey,aKey,cKey,vKey,leftArrow,...
        rightArrow,respKeysRight,respKeysLeft]); 
    
    % Hand to use for responding
    if all(S.hand == 1) % for left hand, hand == 1
        hand = 1;
    elseif all(S.hand == 2) % for right hand, hand == 2
        hand = 2;
    else
        error('Inconsistent response hand definition');
    end
    
    % Defining correct and incorrect response keys
    if hand == 1
        respKeysCorrect = respKeysLeft;
        respKeysInCorrect = respKeysRight;
    elseif hand == 2
        respKeysCorrect = respKeysRight;
        respKeysInCorrect = respKeysLeft;
    end
    
    % sync mode
    syncMode = S.syncMode;
    
    % test mode
    testMode = S.testMode;
    
    % DebugMode
    debug = false;
    if debug
        PsychDebugWindowConfiguration;
        dbstop if error
    end
    
    %%                     OPEN & PREPARE SCREENS
    % =====================================================================
    % Get the list of Screens and choose the one with the highest screen
    % number. Screen 0 is, by definition, the display with the menu bar.
    % Often when two monitors are connected the one without the menu bar is
    % used as the stimulus display.  Chosing the display with the highest
    % dislay number is a best guess about where you want the stimulus
    % displayed.
    screens = Screen('Screens');
    screenNumber = max(screens);
    
    % Load gamma table for Flat panel gamma correction:
    % We load identical gamma tables for red, green and blue channels.
    % Load gamma table for Flat panel gamma correction:
    if ~testMode
        if ~isempty(S.setupSpec.corrGamma)
            oldLut = Screen('LoadNormalizedGammaTable',screenNumber,S.setupSpec.corrGamma);        
        end
    end
    
    % Prepare setup of imaging pipeline for onscreen window.
    % This is the first step in the sequence of configuration steps
    PsychImaging('PrepareConfiguration');
    % Disable psychtoolbox start screen and present a blank black screen
    % instead. 
    Screen('Preference','VisualDebuglevel',3);
    
    % Colors
    [white,black] = deal(WhiteIndex(screenNumber),BlackIndex(screenNumber));
    grey = black + round((white - black)*0.33);
    blue = [0,0,255];
    stimColor = white;
    % black + round((white - black)*0.5);
    
    % Background luminance value for the screen. 
    if syncMode
        backgroundColor = black;
    else
        backgroundColor = grey;
    end
    
    % Color of the feedback cue frame. 
    feedbCueFrameColor = white;
    
    % Open window. 
    [w,winRect] = PsychImaging('OpenWindow',screenNumber,backgroundColor);
        
    HideCursor;
    % prevent that keyboard inputs are inserted somewhere inside the script
    % or the Matlab Window. IMPORTANT: needs to be reset at the end of the
    % script.
    ListenChar(2);
    
    % Define window center. 
    [win_center_x,win_center_y] = RectCenter(winRect);
    
    %%                      FILL VISUAL BUFFERS
    % =====================================================================
    % determines the text size etc.
    Screen('TextFont',w,'Arial');
    Screen('TextSize',w,30);
    Screen('TextStyle',w,1);
    
    % Clear screen to black background color
    Screen('Flip',w,[],0);
    % ...to have a black Screen all the time
    
    % Switch to realtime-priority to reduce timing jitter and interruptions
    % caused by other applications and the operating system itself:
    Priority(MaxPriority(w));
    
    Screen('Flip',w,[],0);
    
    % Measure monitor refresh. This will trigger a calibration loop of 
    % minimum 100 valid samples and return the estimated ifi in 
    % 'ifi': We require an accuracy of 0.05 ms == 0.00005 secs. 
    % If this level of accuracy can't be reached, we time out after 20 
    % seconds. 
    % ifi: interflip interval in seconds.
    if debug
        ifi = 0.016666666667;
    else
        [ifi,nvalid,stddev] = Screen('GetFlipInterval',w,100,0.00005,5);
        fprintf(['Measured refresh interval, as reported by "GetFlipInterval" ' ...
            'is %2.5f ms. (nsamples = %i, stddev = %2.5f ms)\n'],...
            ifi*1000,nvalid,stddev*1000);
    end

    %                           Fixation cross
    %----------------------------------------------------------------------
    % Parameters in azimuth
    [fixWidth,lineWidth] = deal(round(posdeg2pix(1,S.viewDist,S.screenRes)),...
                                round(posdeg2pix(0.1,S.viewDist,S.screenRes))); 
    
    [fixImage,fixRect] = makefixationcross(fixWidth,lineWidth,backgroundColor,white);
    fixCrossWhite.texture = Screen('MakeTexture',w,fixImage);
    fixCrossWhite.rect = CenterRectOnPoint(fixRect,win_center_x,win_center_y);
    
    [fixImage,fixRect] = makefixationcross(fixWidth,lineWidth,backgroundColor,blue);
    fixCrossBlue.texture = Screen('MakeTexture',w,fixImage);
    fixCrossBlue.rect = CenterRectOnPoint(fixRect,win_center_x,win_center_y);
        
    %                           Visual stimuli
    %----------------------------------------------------------------------
    % Defining one single visual stimulus, a sphere
    sphereRect = [0 0 repmat(round(posdeg2pix(S.diameterSphere,S.viewDist,S.screenRes)),1,2)];
    
    % Pre-allocating an array for storing the spheres' coordinates. 
    % nSpheresPerStimulus x 4 x nTrials
    spheres = NaN(size(S.visualStimuli,3),4,size(S.visualStimuli,1));
    
    % Converting visual stimulus locations to pixel coordinates. 
    % Positions are given in azimuth relative to the center of the screen. 
    % In order to convert them to the screen's coordinate system, their 
    % values must be converted to pixel value using posdeg2pix and the 
    % center coordinates of the screen must be added to them.
    for i = 1:size(spheres,3)
        spheres(:,:,i) = CenterRectOnPoint(sphereRect,...
            win_center_x+round(posdeg2pix(squeeze(S.visualStimuli(i,1,:)),S.viewDist,S.screenRes)),...
            win_center_y+round(posdeg2pix(squeeze(S.visualStimuli(i,2,:)),S.viewDist,S.screenRes)));
    end
    
    %                          Feedback cues
    %----------------------------------------------------------------------
    % Width and height of the response cues in azimuth. 
    feedbCueSizeAzim = 2;
    feedbCueRect = [0 0 repmat(round(posdeg2pix(feedbCueSizeAzim,S.viewDist,S.screenRes)),1,2)];
    
    % Width of the frame of the response cues. 
    feedbCueFrameWidthPix = round(posdeg2pix(0.1,S.viewDist,S.screenRes));
    
    feedbCueRects = NaN(S.nTrials,4);
    
    % All possible locations of the response cues converted to pixel
    % coordinates. 
    for i = 1:S.nTrials
        
        if ~isnan(S.locationAuditory(i))
            feedbCueRects(i,:) = CenterRectOnPoint(feedbCueRect,...
                win_center_x+round(posdeg2pix(S.locationAuditory(i),S.viewDist,S.screenRes)),...
                win_center_y);
        end
        
    end
    
    
    %%                       PREPARING AUDIO DEVICE
    %%=====================================================================
    
    % Open sound driver for high timing precision, stereo playback with
    % 44.1 kHZ sampling frequency: 
    % This has to be called with input argument 1 in order to achieve
    % required AV synchrony! For test mode we don't really need low
    % latency. 
    if testMode
        reallyNeedLowLatency = 0;
    else
        reallyNeedLowLatency = 1;
    end
    InitializePsychSound(reallyNeedLowLatency);
    
    % Use sugLat seconds minimum latency on Windows to not overload the
    % system:
    if IsWin
        % According to the PsychPortAudioTimingTest I specify 0.015 s for
        % suggested latency. 
        sugLat = 0.015;
    else
        sugLat = [];
    end
    
    % The measured visuo-audio delay (visual onset minus auditory onset 
    % mesured with PsychPortAudioTimingTest) of the system depends on the 
    % setup we use. This is specific to each setup, hence must be measured
    % for each individual setup we use. 
    latBias = S.setupSpec.latBias; % -0.0067
    
    % reqLatencyClass level 2 means: Take full control over the audio 
    % device, even if this causes other sound applications to fail or
    % shutdown.
    if testMode
        reqLatencyClass = 0;
    else
        reqLatencyClass = 2;
    end
    audioDeviceID = S.setupSpec.audioDeviceID;
    
    % The final volume for the audio devices
    volume = 0.05;
    
    % Open audio device. 
    if strcmp(S.audioParam.mode,'speaker')
        paMaster = PsychPortAudio('Open',audioDeviceID,1+8,reqLatencyClass,S.audioParam.sr,...
            S.audioParam.nChannels,[],sugLat);
        % The visuo-auditory delay compensation takes place here. 
        PsychPortAudio('LatencyBias',paMaster,latBias);
        % Start master immediately:
        PsychPortAudio('Start',paMaster,0,0,1);
        
        paSlave = NaN(1,S.audioParam.nChannels);
        for i = 1:S.audioParam.nChannels
            paSlave(i) = PsychPortAudio('OpenSlave',paMaster,1,1,i);
            PsychPortAudio('Volume',paSlave(i),1.0,volume*S.audioParam.details.channelComp(i));
        end
        
        respCueChannels = find(ismember(S.audioParam.details.channelMap,...
                            S.audioParam.details.auditoryCueLocations));
        
        % Generate some beep sound 1000 Hz, 0.1 secs, 50% amplitude and fill it
        % in the buffer for preheting playback.
        mynoise = 0.5 * MakeBeep(1000,0.1,S.audioParam.sr);
        
        % Preheat: run  audio device once silently, with volume set to zero.
        PsychPortAudio('Volume',paMaster,0);
        for i = 1:numel(paSlave)
            PsychPortAudio('FillBuffer',paSlave(i),mynoise);
            PsychPortAudio('Start',paSlave(i),1,0,1);
            PsychPortAudio('Stop',paSlave(i),1);
        end
        PsychPortAudio('Volume',paMaster,1);
        
    else % itd rec
        paHandle = PsychPortAudio('Open',audioDeviceID,[],reqLatencyClass,S.audioParam.sr,...
            S.audioParam.nChannels,[],sugLat);
        % The visuo-auditory delay compensation takes place here. 
        PsychPortAudio('LatencyBias',paHandle,latBias);
        
        % Generate some beep sound 1000 Hz, 0.1 secs, 50% amplitude and fill it
        % in the buffer for preheting playback.
        mynoise = 0.5 * MakeBeep(1000,0.1,S.audioParam.sr);
        mynoise = repmat(mynoise,S.audioParam.nChannels,1);
        PsychPortAudio('FillBuffer',paHandle,mynoise);
        % Preheat: run  audio device once silently, with volume set to zero.
        PsychPortAudio('Volume',paHandle,0);
        PsychPortAudio('Start',paHandle,1,0);
        PsychPortAudio('Stop',paHandle,1);
        
        PsychPortAudio('Volume',paHandle,volume);
    end
    
    %%                         PREPARING LABJACK 
    %%=====================================================================
    if ~exist('lj','class')
        %open a U3 with verbose logging on 
        lj = labJack('deviceID',3,'verbose',true);
        
    end
    
    if strfind(lj.version,'FAILED')
        % error('Opening labjack failed! ');
        warning('Opening labjack failed! Unable to send triggers! ');
        ljPresent = false;
    else
        ljPresent = true;
    end
    
    %%                      PREPARING THE EYETRACKER
    %%=====================================================================
    S.disp.grey = grey;
    S.disp.white = white;
    S.disp.win = w;
    
    S.eyelink.dispRect = CenterRectOnPoint(round(posdeg2pix([0,0,32,25],S.viewDist,S.screenRes)),win_center_x,win_center_y);
    S.eyelink.edfFileName = ...
        fullfile(DEC_2_setupdir(S.expStage,'data_behav_sub',S.subID),...
        [S.subID,'_',S.expMode,'_',num2str(S.iDay),'_',num2str(S.iSessionInDay),...
        '_','eyelink','_',datestr(now,'ddmmyyyy_HHMM'),'.edf']);
    
    S = runeyelink('setup',S);
    
    %%                       DEFINING TRIGGER CODES
    %%=====================================================================
    % 1 - 76  = corresponding conditions
    % Trigger code for indicating the start of the run. 
    if hand == 1
        trStartSession = 101;
    elseif hand == 2
        trStartSession = 102;
    end
    % Trigger code for indicating the start of the trial. 
    trStartTrial = 100;
    
    trRespCue = 103;
    % Vector of trigger codes for indicating the answers. 
    trResponse = 110+(1:length(respKeysRight))-1;
    
    
    %%                PREALLOCATE DATA COLLECTION ARRAYS
    % =====================================================================
    % Pre-allocate an array for storing iformation from all presented 
    % frames in the whole run. 
    % The estimated duration of the session in seconds.
    if S.feedback
        sessionDuration = S.nTrials*(max(S.timing.t_toStimOnset) + S.timing.t_stimDuration + ...
            S.timing.t_toResponse + S.timing.t_endResponse + S.timing.t_feedbackDuration);
    else
        sessionDuration = S.nTrials*(max(S.timing.t_toStimOnset) + S.timing.t_stimDuration + ...
            S.timing.t_toResponse + S.timing.t_endResponse);
    end
    % The estimated number of frames in the session with an extra 10 s to
    % be on the safe side: duration of the session x refresh rate (Hz). 
    estimatedNumOfFrames = ceil((sessionDuration + 10)/ifi);
    % ts: vbl timestamps
    % beampos: beam positions
    % missest: misses
    % flipfin: the time when flip is returned to MatLab
    % td: time drawing is finished
    % so: estimated stimulus onset time aka end of VBL
    [ts,beampos,missest,flipfin,td,fo] = deal(NaN(estimatedNumOfFrames,1));
    
    % Arrays for collecting timing details from trials. 
    % t_startTrial: The time of the first flip to which every other flips are 
    %   timed. (Not identical to the start of the trial as it is the 
    %   appearence of the fixation dot).
    % t_supposedStartVis: Supposed start time of visual stimulus. 
    % t_actualStartVis: The time when the visual stimulus actually started. 
    % t_actualStartAud: The time of the actual auditory stimulus onset. 
    [t_startTrial,t_supposedStimulusOnset,t_actualStartVis,...
        t_actualStartAud,t_actualStartAudCue] = deal(NaN(S.nTrials,1));
    
    % Arrays for collecting answers. 
    [locResp,locRespTimes] = deal(NaN(S.nTrials,1));
    wrongHand = false(S.nTrials,1);
    
    %%                          START SESSION
    %%=====================================================================
    
    % Number of monitor refresh intervals before flipping front- and 
    % backbuffer. 
    numifis = 1;
    % We have to let the system know a couple of refresh cycles earlier
    % before auditory onset that we're going to start the audio device
    % at a particular time to maximize the possibility of the system
    % hitting the onset deadline. This is the amount of time (in
    % frames) before stimulus onset, we tell the system to schedule
    % the presentation of the auditory stimulus. This is a fairly
    % conservative value, possibly less frames wolud do as well, but I
    % haven't tested. 
    audioOnsetBuffer = ceil((2*sugLat)/ifi)+1;
    
    % Running index counting the presented frames from the start of the
    % run. 
    iFramesInSession = 0;
    
    % Instructions. 
    if ~syncMode
        if S.task == 1
            task = 'sound';
        else
            task = 'center of the cloud';
        end
        
        if hand == 1
            handStr = 'LEFT';
        elseif hand == 2
            handStr = 'RIGHT';
        end
        
        line1 = sprintf('Please report the location of the %s! ',task);
        line2 = sprintf('Use your %s hand to answer! ',handStr);
        DrawFormattedText(w,sprintf('%s\n\n%s\n',line1,line2),'center','center',white);
        Screen(w,'Flip');
        
        WaitSecs(5);
    end
    
    % Start recroding eyetracking data
    runeyelink('record',S);
    
    % Drawing the fixation cross. 
    Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
    Screen(w,'Flip');
    
    WaitSecs(2);
    
    for iTrial=1:S.nTrials
        %                  Initializing variables
        %------------------------------------------------------------------
        % Converting timing parameters from seconds to frames. 
        f_toStimOnset = round(S.timing.t_toStimOnset(iTrial)/ifi);
        f_stimDuration = round(S.timing.t_stimDuration/ifi);
        f_toResponse = round(S.timing.t_toResponse/ifi);
        f_respCueDuration = round(S.timing.t_respCueDuration/ifi);
        f_endResponse = round(S.timing.t_endResponse/ifi);
        f_feedbackDuration = round(S.timing.t_feedbackDuration/ifi);
        
        % Flag indicating the task being done. 
        responseDone = false;
        % Flag indicating if feedback is done. 
        feedbackDone = false;
        
        % Number of frames presented in the trial. 
        iFramesInTrial = 0;
        % Number of frames presented for the stimulus. 
        f_presentedStim = 0;
        % Number of frames presented for the response. 
        f_presentedResp = 0;
        % Number of frames presented for the response. 
        f_presentedFeedb = 0;
        
        %                         Fill sound buffer. 
        %------------------------------------------------------------------
        if ~isnan(S.locationAuditory(iTrial))
            if strcmp(S.audioParam.mode,'speaker')
                PsychPortAudio('FillBuffer',...
                    paSlave(S.audioParam.details.channelMap == S.locationAuditory(iTrial)),...
                    squeeze(S.auditoryStimuli(iTrial,:,:))');
            else
                PsychPortAudio('FillBuffer',paHandle,squeeze(S.auditoryStimuli(iTrial,:,:)));
            end
        end
        
        %                           Synchronizing
        %------------------------------------------------------------------
        % Drawing the fixation cross. 
        Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
        % Perform some initial Flip to get us in sync with retrace:
        % tvbl is the timestamp (system time in seconds) when the retrace
        % started. We need it as a reference value for our display timing:
        tvbl = Screen('Flip',w);
                
        % SEND TRIGGER HERE
        if ~syncMode
            if iTrial == 1
                if ljPresent
                    % Prepare a trigger output with value given in trig.
                    lj.prepareStrobe(trStartSession);
                    % Send the strobed word, which is 11 bit max length on EIO
                    % and CIO via the DB15.
                    lj.strobeWord;
                    lj.prepareStrobe(trStartTrial);
                    lj.strobeWord;
                end
                if S.eyelink.online
                    Eyelink('Message','%ld',trStartSession);
                    Eyelink('Message','%ld',trStartTrial);
                end
            else
                if ljPresent
                    lj.prepareStrobe(trStartTrial);
                    lj.strobeWord;
                end
                if S.eyelink.online
                    Eyelink('Message','%ld',trStartTrial);
                end
            end
        end
        
        t_startTrial(iTrial) = (tvbl+ifi);
        
        % Check input
        [~,~,keyCode] = KbCheck;
        
        
        %%                           START TRIAL                      
        %%=================================================================
        while 1
            % Keep track of the number of frames presented in the run. 
            iFramesInSession = iFramesInSession + 1;
            
            % Keep track of the number of frames presented in the trial. 
            iFramesInTrial = iFramesInTrial + 1;
            
            % The deadline for the next stimulus presentation. We substract
            % 0.5*ifi so the screen gets the flip command a bit earlier 
            % and flips on the next possibility.
            t_deadline = tvbl + numifis * ifi - 0.5 * ifi;
            
            % If user supplied numifis=0, we force t_deadline_visual=0,
            % so Flip will actually ignore the deadline and just Flip
            % at the next possible retrace...
            if numifis == 0
                t_deadline = 0;
            end
            
            % Reset condition trigger to zero
            trig = 0;
            
            
            %                    Pre-stimulus fixation
            %--------------------------------------------------------------
            if iFramesInTrial < f_toStimOnset
                Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
                
                % To provide the system enough time to prepare the sound
                % device for starting, I send out the start command a
                % couple of video refresh cycles earlier (audioOnsetBuffer).
                % The visual stimuli will be presented at around half the
                % height of the monitor, so I delay sound onsets relative to
                % the visual stimuli with half an ifi.
                if iFramesInTrial == f_toStimOnset - audioOnsetBuffer
                    % Start the sound device if applicable.
                    if ~isnan(S.locationAuditory(iTrial))
                        if strcmp(S.audioParam.mode,'speaker')
                            PsychPortAudio('Start',...
                                paSlave(S.audioParam.details.channelMap == S.locationAuditory(iTrial)),...
                                [],(visOnsetAtAudioOnsetBufferBeforeDeadline + (audioOnsetBuffer + 0.5) * ifi));
                        else
                            PsychPortAudio('Start',paHandle,[],...
                                (visOnsetAtAudioOnsetBufferBeforeDeadline + (audioO] h,.nsetBuffer + 0.5) * ifi));
                        end
                    end
                    % Saving timestamp for supposed AudioVisual stimulus onset.
                    t_supposedStimulusOnset(iTrial) = ...
                        (visOnsetAtAudioOnsetBufferBeforeDeadline + (audioOnsetBuffer + 0.5) * ifi);
                end
                
            %                     Stimulus presentation
            % -------------------------------------------------------------
            elseif  iFramesInTrial >= f_toStimOnset && f_presentedStim <= f_stimDuration
                
                f_presentedStim = f_presentedStim + 1;
                
                % Setting the trigger code. 
                if iFramesInTrial == f_toStimOnset
                    trig = S.condition(iTrial);
                end
                
                % Fixation cross
                Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
                
                % Visual or audiovisual tirals. 
                if ~isnan(S.locationVisual(iTrial))
                    
                    if syncMode
                        % Drawing blank white screen
                        Screen('FillRect',w,white,winRect);
                    else
                        % Drawing spheres.
                        Screen('FillOval',w,stimColor,spheres(:,:,iTrial)');
                    end
                end
                % For auditory only trials the fixation cross alone is
                % presented. 
                
            %                    Post-stimulus fixation
            % -------------------------------------------------------------
            elseif iFramesInTrial < f_toStimOnset + f_stimDuration + f_toResponse
                
                % Resetting screen to black if in sync mode
                if syncMode
                    Screen('FillRect',w,black,winRect);
                end
                
                % Fixation cross
                Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
                
                if iFramesInTrial == f_toStimOnset + f_stimDuration + f_toResponse - audioOnsetBuffer
                    
                    % Log audio stimulus onset time and stop audio device
                    % if applicable
                    if ~isnan(S.locationAuditory(iTrial))
                        if strcmp(S.audioParam.mode,'speaker')
                            % Log sound onset time:
                            audiostatus = PsychPortAudio('GetStatus',...
                                paSlave(S.audioParam.details.channelMap == S.locationAuditory(iTrial)));
                            t_actualStartAud(iTrial) = audiostatus.StartTime;
                            % Stop soundplayback:
                            PsychPortAudio('Stop',paSlave(S.audioParam.details.channelMap == S.locationAuditory(iTrial)));
                        else
                            % Log sound onset time:
                            audiostatus = PsychPortAudio('GetStatus',paHandle);
                            t_actualStartAud(iTrial) = audiostatus.StartTime;
                            % Stop soundplayback:
                            PsychPortAudio('Stop',paHandle);
                        end
                    end
                    
                    % Present response cue if applicable
                    if f_respCueDuration > 0
                        % Fill audio device with some noise
                        if strcmp(S.audioParam.mode,'speaker')
                            for i = 1:numel(respCueChannels)
                                PsychPortAudio('FillBuffer',...
                                    paSlave(respCueChannels(i)),mynoise);
                            end
                        else
                            PsychPortAudio('FillBuffer',paHandle,mynoise);
                        end
                        
                        
                        % Start the sound device if applicable.
                        if strcmp(S.audioParam.mode,'speaker')
                            for i = 1:numel(respCueChannels)
                                PsychPortAudio('Start',...
                                    paSlave(respCueChannels(i)),...
                                    [],(visOnsetAtAudioOnsetBufferBeforeDeadline + (audioOnsetBuffer + 0.5) * ifi));
                            end
                        else
                            PsychPortAudio('Start',paHandle,[],...
                                (visOnsetAtAudioOnsetBufferBeforeDeadline + (audioOnsetBuffer + 0.5) * ifi));
                        end
                    end
                    
                end
                
                
            %                        Response period
            %--------------------------------------------------------------
            elseif iFramesInTrial >= f_toStimOnset + f_stimDuration + f_toResponse
                
                f_presentedResp = f_presentedResp + 1;
                
                % Fixation cross
                if f_presentedResp <= f_respCueDuration
                    Screen('DrawTexture',w,fixCrossBlue.texture,[],fixCrossBlue.rect);
                else
                    Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
                end
                
                if iFramesInTrial == f_toStimOnset + f_stimDuration + f_toResponse
                    if ~syncMode, trig = trRespCue; end
                end
                
                if ~S.feedback 
                    
                    if f_presentedResp <= f_endResponse
                        
                        % Check the keyboard to see if a button has been pressed
                        [keyIsDown,secs,keyCode] = KbCheck;
                        keyIDx = find(keyCode,1,'first');
                        
                        if ~responseDone
                            
                            if keyIsDown && any(respKeysCorrect == keyIDx)
                                % Send response trigger immediately
                                if ljPresent
                                    lj.prepareStrobe(trResponse(respKeysCorrect == keyIDx));
                                    lj.strobeWord;
                                end
                                if S.eyelink.online
                                    Eyelink('Message','%ld',trResponse(respKeysCorrect == keyIDx));
                                end
                                locResp(iTrial) = S.auditoryLocationLevels(respKeysCorrect == keyIDx);
                                locRespTimes(iTrial) = secs;
                                responseDone = true;
                            elseif keyIsDown && any(respKeysInCorrect == keyIDx)
                                % Send response trigger immediately
                                if ljPresent
                                    lj.prepareStrobe(trResponse(respKeysInCorrect == keyIDx));
                                    lj.strobeWord;
                                end
                                if S.eyelink.online
                                    Eyelink('Message','%ld',trResponse(respKeysInCorrect == keyIDx));
                                end
                                locResp(iTrial) = S.auditoryLocationLevels(respKeysInCorrect == keyIDx);
                                locRespTimes(iTrial) = secs;
                                wrongHand(iTrial) = true;
                                responseDone = true;
                            end
                            
                        end
                        
                    end
                
                % Feedback period if applicable
                else
                    
                    if ~responseDone && f_presentedResp <= f_endResponse
                        
                        % Check the keyboard to see if a button has been pressed
                        [keyIsDown,secs,keyCode] = KbCheck;
                        keyIDx = find(keyCode,1,'first');
                        
                        if ~responseDone
                            
                            if keyIsDown && any(respKeysCorrect == keyIDx)
                                % Send response trigger immediately
                                if ljPresent
                                    lj.prepareStrobe(trResponse(respKeysCorrect == keyIDx));
                                    lj.strobeWord;
                                end
                                if S.eyelink.online
                                    Eyelink('Message','%ld',trResponse(respKeysCorrect == keyIDx));
                                end
                                locResp(iTrial) = S.auditoryLocationLevels(respKeysCorrect == keyIDx);
                                locRespTimes(iTrial) = secs;
                                responseDone = true;
                            elseif keyIsDown && any(respKeysInCorrect == keyIDx)
                                % Send response trigger immediately
                                if ljPresent
                                    lj.prepareStrobe(trResponse(respKeysInCorrect == keyIDx));
                                    lj.strobeWord;
                                end
                                if S.eyelink.online
                                    Eyelink('Message','%ld',trResponse(respKeysInCorrect == keyIDx));
                                end
                                locResp(iTrial) = S.auditoryLocationLevels(respKeysInCorrect == keyIDx);
                                locRespTimes(iTrial) = secs;
                                wrongHand(iTrial) = true;
                                responseDone = true;
                            end
                            
                        end
                        
                    else
                        
                        f_presentedFeedb = f_presentedFeedb + 1;
                        
                        if f_presentedFeedb < f_feedbackDuration
                            if ~isnan(S.locationAuditory(iTrial))
                                Screen('FrameRect',w,feedbCueFrameColor,feedbCueRects(iTrial,:),feedbCueFrameWidthPix);
                            end
                        else
                            feedbackDone = true;
                        end
                        
                    end
                end
                
            end
            
            
            % Tell the computer that no more stimuli will be drawn, in 
            % order to make it faster. 
            td(iFramesInSession) = Screen('DrawingFinished',w,0);
            
            %                             FLIP
            % -------------------------------------------------------------
            [tvbl,visOnset,flipfin(iFramesInSession),...
                missest(iFramesInSession),beampos(iFramesInSession)] = ...
                Screen('Flip',w,t_deadline);
            
            % Send trigger. 
            if trig ~= 0
                if ljPresent
                    lj.prepareStrobe(trig);
                    lj.strobeWord;
                    % Trigger code will be resetted at the beginning of the
                    % loop.
                end
                if S.eyelink.online
                    Eyelink('Message','%ld',trig);
                end
            end
            
            % Record timestamp for vbl and visual onset respectively for
            % later use: 
            ts(iFramesInSession) = tvbl;
            fo(iFramesInSession) = visOnset;
            
            % Saving the vbl timestamp watiFramesAud frames before the
            % the desired stimulus onset for scheduling sound onset. 
            if iFramesInTrial == f_toStimOnset - (audioOnsetBuffer + 1)
                % Quite a crappy name for a variable I know, but I couldn't
                % come up with a better one... 
                visOnsetAtAudioOnsetBufferBeforeDeadline = visOnset;
            elseif iFramesInTrial == f_toStimOnset + f_stimDuration + f_toResponse - (audioOnsetBuffer + 1)
                visOnsetAtAudioOnsetBufferBeforeDeadline = visOnset;
            end
            
            % Saving vbl timestamp for the stimulus onset. 
            if iFramesInTrial == f_toStimOnset
                t_actualStartVis(iTrial) = tvbl;
            end
            
            
            if S.feedback
                % Terminate trial loop if feedback is done.
                if feedbackDone
                    break;
                end
            else
%                 % Terminate trial loop immediately after response. 
%                 if responseDone
%                     break;
%                 elseif iFramesInTrial == f_toStimOnset + f_stimDuration + f_toResponse + f_endResponse
%                     break;
%                 end

                % Also terminate trial loop a certain amount of time after the
                % answering period has started.
                if iFramesInTrial == f_toStimOnset + f_stimDuration + f_toResponse + f_endResponse
                    break;
                end
            end
            
            [~,~,keyCode] = KbCheck;
            % Stopping the loop if user wants to quit. 
            if keyCode(quitKey)
                break;
            end
   
        end % Draw next frame...
        
        %                        END TRIAL LOOP                               
        %------------------------------------------------------------------        
        if strcmp(S.audioParam.mode,'speaker')
            % Log sound onset time:
            for i = 1:numel(respCueChannels)
                audiostatus = PsychPortAudio('GetStatus',paSlave(respCueChannels(i)));
            end
            t_actualStartAudCue(iTrial) = audiostatus.StartTime;
            % Stop soundplayback:
            for i = 1:numel(respCueChannels)
                PsychPortAudio('Stop',paSlave(respCueChannels(i)));
            end
            
        else
            % Log sound onset time:
            audiostatus = PsychPortAudio('GetStatus',paHandle);
            t_actualStartAudCue(iTrial) = audiostatus.StartTime;
            % Stop soundplayback:
            PsychPortAudio('Stop',paHandle);
        end
        
        if keyCode(quitKey)
            break;
        end
        
        % Drawing the fixation cross. 
        Screen('DrawTexture',w,fixCrossWhite.texture,[],fixCrossWhite.rect);
        Screen(w,'Flip');
        
    %                          END SESSION LOOP
    %----------------------------------------------------------------------
    end
    
    
    %%                             TIMING
    %%=====================================================================
    % Number of cycles (in loops) in ts are the timestamps. 
    n = iFramesInSession;
    % Count and output number of missed flip on VBL deadlines:
    numbermisses = 0;
    if numifis > 0
        for iFramesInSession = 2:n
            if (ts(iFramesInSession)-ts(iFramesInSession-1) > ifi*(numifis+0.5))
                numbermisses = numbermisses+1;
            end
        end
    else
        for iFramesInSession = 2:n
            if (ts(iFramesInSession)-ts(iFramesInSession-1) > ifi*1.5)
                numbermisses = numbermisses+1;
            end
        end
    end
    
    % Output some summary and say goodbye...
    fprintf('PTB missed %i out of %i stimulus presentation deadlines.\n',...
            numbermisses,n);
    
    
    %%                      COLLECTING OUTPUT DATA
    %%=====================================================================
    
    % general information
    info.iDay                           = S.iDay;
    info.inputParameter                 = S;
    info.inputParameter.presentFuncName = presentFunctName;
    info.iSessionInDay                  = S.iSessionInDay;
    info.iSessionOverall                = S.iSessionOverall;
    info.subID                          = S.subID;
    
    % timing
    info.timing.actAuditoryStartTime      = t_actualStartAud;
    info.timing.actVisualStartTime        = t_actualStartVis;
    info.timing.beampos                   = beampos(1:iFramesInSession);
    info.timing.drawingFinised            = td(1:iFramesInSession);
    info.timing.flipOnsets                = flipfin(1:iFramesInSession);
    info.timing.ifi                       = ifi;
    info.timing.latBias                   = latBias;
    info.timing.missest                   = missest(1:iFramesInSession);
    info.timing.nFramesPresentedInSession = iFramesInSession;
    info.timing.frameOnsets               = fo(1:iFramesInSession);
    info.timing.suppStimulusOnset         = t_supposedStimulusOnset;
    info.timing.suppTrialStartTime        = t_startTrial;
    info.timing.vblTimeStamps             = ts(1:iFramesInSession);
    
    % responses
    info.resp.locResp      = locResp;
    info.resp.locRespTimes = locRespTimes;
    info.resp.wrongHand    = wrongHand;
    
    if any(wrongHand)
        warning('Participant used the wrong hand for response!');
    end
    
    % Cleaning up the most spacious and duplicat variables before saving 
    % all data. 
    [ts,beampos,missest,flipfin,td,fo] = deal([]);
    S.auditoryStimuli = [];
    
    % Saving eyetracking data and closing eyetracker
    runeyelink('close',S);
    
    % Saving experiment data
    savedfname = fullfile(DEC_2_setupdir(S.expStage,'data_behav_sub',S.subID),...
        [S.subID,'_',S.expMode,'_',num2str(S.iDay),'_',num2str(S.iSessionInDay),...
        '_','all_present','_',datestr(now,'ddmmyyyy_HHMM'),'.mat']);
    
    if exist(savedfname,'file') == 2
        warning('Attempted to overwrite already existing file! ');
        savedfname = fullfile(DEC_2_setupdir(S.expStage,'data_behav_sub',S.subID),...
            [S.subID,'_',S.expMode,'_',num2str(S.iDay),'_',num2str(S.iSessionInDay),...
            '_','all_present','_',datestr(now,'ddmmyyyy_HHMM'),'_warn','.mat']);
    end
    save(savedfname)
    
    
    %%                          END PRESENTATION
    %%=====================================================================
    % Shutdown realtime scheduling:
    Priority(0);
    
    % Shutdown sound driver:
    PsychPortAudio('Close');
    
    % Close display: If we skipped/missed any presentation deadline during
    % Flip, Psychtoolbox will automatically display some warning message on
    % the Matlab console:
    Screen('CloseAll');
    
    % Restore original system gamma table
    if ~testMode
        if ~isempty(S.setupSpec.corrGamma)
            Screen('LoadNormalizedGammaTable',screenNumber,oldLut);        
        end
    end
    
    ShowCursor;
    ListenChar;
    RestrictKeysForKbCheck([]);
    
    
catch %#ok<CTCH>
    % This "catch" section executes in case of an error in the "try" 
    % section above. Importantly, it closes the onscreen window if its 
    % open and shuts down realtime-scheduling of Matlab:
    
    % Disable realtime-priority in case of errors.
    Priority(0);
    
    % Shutdown sound driver 
    PsychPortAudio('Close');
    
    % Close display 
    Screen('CloseAll');
    
    % Restore original system gamma table
    if ~testMode
        if ~isempty(S.setupSpec.corrGamma)
            Screen('LoadNormalizedGammaTable',screenNumber,oldLut);        
        end
    end
    
    ShowCursor;
    ListenChar;
    RestrictKeysForKbCheck([]);
    
    psychrethrow(psychlasterror);
    
end %try..catch..

return