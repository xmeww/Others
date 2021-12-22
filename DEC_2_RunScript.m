clear all; clc;
commandwindow;

% Checking matlab version
if verLessThan('matlab', '8.3.0.532')
    error('Please start MATLAB 2014a in order to run the experiment! ');
end

%%                           Getting input
%%=========================================================================
% Request parameters via dialog box
prompt = {'Subject ID','Experiment mode','Audio mode','Session order','Day','Session to start with','Sessionf for day'};
title = 'Input parameters';
lines = [1 35]; % 1 row, 35 chars in each line
defaults = {'test','exp','speaker','1','1','1',''};
params = inputdlg(prompt,title,lines,defaults);
if ~isempty(params) % ok is pressed
    subID = params{1};
    expMode = params{2};
    audioMode = params{3};
    params = cellfun(@str2num,params(4:end),'UniformOutput',false);
    [actSessionOrder,iDay,sessionToStartWith,sessionsForCurrDay] = params{:};
else % cancel is pressed
    fprintf('Run has been aborted...\n')
    return
end

% Checking input. 
checkExpMode = @(x) any(validatestring(x,{'exp','train'}));
checkAudioMode = @(x) any(validatestring(x,{'rec','itd','speaker'}));
checkPositiveScalar = @(x) isscalar(x) && x > 0;

if ~isempty(sessionsForCurrDay)
    validateattributes(sessionsForCurrDay,{'numeric'},{'row','positive','finite','integer','increasing'});
end

if ~checkExpMode(expMode)
    error('Wrong input for ''Experiment mode''! ');
elseif ~checkAudioMode(audioMode)
    error('Wrong input for ''Audio mode''! ');
elseif ~checkPositiveScalar(sessionToStartWith)
    error('Wrong input for ''Session to start with''! ');
elseif ~checkPositiveScalar(iDay)
    error('Wrong input for ''Day''! ');
end


%%                   Setting up experiment environment. 
%%=========================================================================
% Testing mode
testMode = false;
if testMode
    warning('Experiment runs in test mode! ');
end

% Synchrony measurement mode
syncMode = false;
if syncMode
    warning('Experiment runs in sync mode! ');
end

% Get file name of currently running function
runScriptName = mfilename;

% Experiment stage. 
expStage = 'final';

% Adding MatlabGarbageCollector.jar to the dynamic java path
javaaddpath(fullfile(DEC_2_setupdir(expStage,'pres'),'MatlabGarbageCollector.jar'));

% always needs to be put if using KbName!!!!!
KbName('UnifyKeyNames');

% to exit the program
quit = KbName('ESCAPE');

%                   Loading setup specific information
%--------------------------------------------------------------------------
setupID = getenv('computername');

temp = load('setup_spec');
setupIDlist = [temp.setup_spec.ID];
idx = find(strcmp(setupID,setupIDlist));
if ~isempty(idx)
    setupSpec = temp.setup_spec(idx);
else
    error('Could not find setup specific information! ');
end


%%                      Defining experiment parameters 
%%=========================================================================

%                     Stimulus parameters and conditions
%--------------------------------------------------------------------------
% The locations of the auditory and visual stimuli in horizontal plane. 
% Unit: azimuth (angle in the horizontal plane).
[auditoryLocationLevels,visualLocationLevels,respCursorLocationLevels] = deal([-10 -3.33 3.33 10]);

% The viewing distance in mm. 
viewDist = setupSpec.view_dist;

% Screen widht in mm. 
reqSreenWidth = setupSpec.screen_coords_mm(3);
% Screen height required in mm.
reqSreenHeight = setupSpec.screen_coords_mm(4);

% Screen widht in mm. 
reqScreenWidthPix = setupSpec.screen_coords_pix(3);
% Screen height required in pixels. 
reqScreenHeightPix = setupSpec.screen_coords_pix(4);
% The resolution required for the experiment in pixels/mm.
reqScreenRes = reqScreenWidthPix/reqSreenWidth; 

% The maximum visual angle required for the experiment
% How to calculate (right half): the most extreme right location +
% 2*greatest stimulus SD + right delta value. The left is the same degree
% in the other direction evidently. (azimuth)
maxVisAngle = 68;

% Get the actual screen size in pixels. 
screens = Screen('Screens');
scrInfo = Screen('Resolution',max(screens));

% Checking if the specified screen is suitable for the experiment
if scrInfo.width ~= reqScreenWidthPix || scrInfo.height ~= reqScreenHeightPix
    warning('The specified screen is not suitable for the experiment. ');
end

% Checking if the specified viewing distance is suitable for the experiment
if round(pospix2deg(reqScreenWidthPix/2,viewDist,reqScreenRes)) < maxVisAngle/2
    warning('The specified viewing distance is not suitable for the experiment. ');
end

% The horizontal STDs for gaussian blobs. This determines the reliability
% of the stimulus, hence the use of this word as a synonym. The bigger the 
% STD, the less the reliability. (azimuth)
highVisReliability = 2; 
lowVisReliability = 12;
reliabilityLevels = [highVisReliability lowVisReliability];
% The vertical STD for gaussian blobs. This is fixed for all reliabilities
% (i.e. horizontal STDs)
vertSTDVisual = 2;

% Number of spheres in one visual stimulus (gaussian cloud)
numOfSpheres = 20;

% Diameter of spheres (azimuth)
diameterSphere = 0.43;

% Number of audio channels to use. 
if strcmp(audioMode,'speaker')
    nAudioChannels = 4;
else
    nAudioChannels = 2;
end

% Audio sampling frequency in Hz. 
audioSamplingRate = 44100;

% Stimulus duration in seconds. 
stimulusDuration = 0.05;

% Task levels. 1 - auditory report, 2 - visual report
auditoryReport = 1;
visualReport = 2;
taskLevels = [auditoryReport visualReport];

% Coding of hands
leftHand = 1;
rightHand = 2;

% Number of possible stimulus locations. 
nLocationsStimuli = length(auditoryLocationLevels);

% Number of possible reliability levels. 
nLevelsReliability = length(reliabilityLevels);



%                             Trial numbers
%--------------------------------------------------------------------------

if syncMode
    % Number of days
    nDays = 1;
    % Overriding the iDay
    iDay = 1;
    % Number of sessions per day
    nSessionsPerDay = 1;
    % Total number of trials per session.
    nTrialsPerSession = 100;
    
else
    % Number of days
    nDays = 3;
    
    if strcmp(expMode,'exp')
        % Number of sessions per day
        nSessionsPerDay = 20;
        nSessionsPerDayBisAud = 8;
        nSessionsPerDayBisVis = 8;
        nSessionsPerDayUnisAud = 2;
        nSessionsPerDayUnisVis = 2;
        % Total number of trials per session.
        nTrialsPerSession = 128;
        % Number of trials per condition.
        % In this case this is exactly the number of trials per conditions. 
        nTrialsPerCond = 96;
        
    elseif strcmp(expMode,'train')
        % Number of sessions per day
        nSessionsPerDay = 5;
        nSessionsPerDayBisAud = 1;
        nSessionsPerDayBisVis = 1;
        nSessionsPerDayUnisAud = 2;
        nSessionsPerDayUnisVis = 1;
        % Total number of trials per session.
        nTrialsPerSession = 64;
        % Number of trials per condition.
        % In this case this is the maximum of the number of trials per 
        % condition, just as a rough guide. Not all conditions are going to
        % be presented this many times.
        nTrialsPerCond = 48;
        
    end
    
end

% % Total number of sessions
% nSessions = nDays*nSessionsPerDay;

%      Loading or generating file containing all trials for a session
%--------------------------------------------------------------------------
if ~syncMode
    
    trialFname = fullfile(DEC_2_setupdir(expStage,'data_behav_sub',subID),[subID,'_',expMode,'_allTrials','.mat']);
    
    if exist(trialFname,'file')
        
        load(trialFname);
        
    else
        
        % Reseeding the random number generator
        rng('shuffle');
        
        % Design for bisensory trials.
        setsBis = {visualLocationLevels,auditoryLocationLevels,reliabilityLevels,taskLevels};
        coordsBis = cell(1,4);
        [coordsBis{:}] = ndgrid(setsBis{:});
        
        locationAuditory = coordsBis{2}(:);
        locationVisual = coordsBis{1}(:);
        reliabilityVisual = coordsBis{3}(:);
        task = coordsBis{4}(:);
        condition = (1:size(locationAuditory,1))';
        
        designBis = table(condition,task,locationAuditory,locationVisual,reliabilityVisual);
        
        % Design for unisensory trials.
        % Auditory only trials.
        setsUnisAud = {1,auditoryLocationLevels,reliabilityLevels,taskLevels(1)};
        coordsUnisAud = cell(1,4);
        [coordsUnisAud{:}] = ndgrid(setsUnisAud{:});
        
        locationAuditory = coordsUnisAud{2}(:);
        locationVisual = NaN(size(coordsUnisAud{1}(:)));
        reliabilityVisual = NaN(size(coordsUnisAud{3}(:)));
        task = coordsUnisAud{4}(:);
        condition = ((size(designBis.condition,1) + 1):(size(designBis.condition,1)+size(auditoryLocationLevels,2)))';
        condition = repmat(condition,size(reliabilityLevels,2),1);
        
        designUnisAud = table(condition,task,locationAuditory,locationVisual,reliabilityVisual);
        
        % Visual only trials.
        setsUnisVis = {visualLocationLevels,1,reliabilityLevels,taskLevels(2)};
        coordsUnisVis = cell(1,4);
        [coordsUnisVis{:}] = ndgrid(setsUnisVis{:});
        
        locationAuditory = NaN(size(coordsUnisVis{2}(:)));
        locationVisual = coordsUnisVis{1}(:);
        reliabilityVisual = coordsUnisVis{3}(:);
        task = coordsUnisVis{4}(:);
        condition = ((size(designBis.condition,1)+size(auditoryLocationLevels,2)+1):...
            (size(designBis.condition,1)+size(auditoryLocationLevels,2)+size(locationVisual,1)))';
        
        designUnisVis = table(condition,task,locationAuditory,locationVisual,reliabilityVisual);
        
        designFull = cat(1,designBis,designUnisAud,designUnisVis);
        
        % Generating all the stimuli for the whole experiment
        allTrialsBis = repmat(designBis,nTrialsPerCond,1);
        allTrialsBisAud = allTrialsBis(allTrialsBis.task == auditoryReport,:);
        allTrialsBisVis = allTrialsBis(allTrialsBis.task == visualReport,:);
        clear allTrialsBis;
        allTrialsUnisAud = repmat(designUnisAud,nTrialsPerCond,1);
        allTrialsUnisVis = repmat(designUnisVis,nTrialsPerCond,1);
        
        % Session types coding
        bis_a_L = 1;
        bis_a_R = 2;
        bis_v_L = 3;
        bis_v_R = 4;
        uni_a_L = 5;
        uni_a_R = 6;
        uni_v_L = 7;
        uni_v_R = 8;
        
        % Defining the sessions. This is going to be the actual session order
        % for the training sessions. For the experimental sessions this is
        % going to be randomized.
        if strcmp(expMode,'train')
            sessions = [[uni_a_L,uni_a_R,uni_v_L,bis_a_L,bis_v_R]',...
                        [uni_a_R,uni_a_L,uni_v_R,bis_v_R,bis_a_L]',...
                        [uni_a_L,uni_a_R,uni_v_L,bis_a_R,bis_v_L]'];
        elseif strcmp(expMode,'exp')
            load('session_order.mat');
            sessions = NaN(nSessionsPerDay,nDays);
            % All session orders should look like this when sorted. 
            sortedTemplate = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8];
            if any(sort(sessionOrder(actSessionOrder,:)) ~= sortedTemplate)
                error('Current session order does not match the template!')
            end            
            sessions(:,1) = sessionOrder(actSessionOrder,:)';
            temp = sessions(:,1);
            temp(1:10) = flipud(sessions(11:20,1));
            temp(11:20) = flipud(sessions(1:10,1));
            sessions(:,2) = temp;
            temp(1:5) = flipud(sessions(1:5,1));
            temp(6:10) = flipud(sessions(6:10,1));
            temp(11:15) = flipud(sessions(11:15,1));
            temp(16:20) = flipud(sessions(16:20,1));
            sessions(:,3) = temp;
        end
        
        % Checking if the number of trials per session is divisible by the
        % minimum number of trials covering all possible conditions for each
        % session type.
        if mod(nTrialsPerSession,size(designBis(designBis.task == auditoryReport,:),1))
            error('nTrialsPerSession is not divisible by number of bisAud conditions!');
        elseif mod(nTrialsPerSession,size(designBis(designBis.task == visualReport,:),1))
            error('nTrialsPerSession is not divisible by number of bisVis conditions!');
        elseif mod(nTrialsPerSession,size(designUnisAud,1))
            error('nTrialsPerSession is not divisible by number of unisAud conditions!');
        elseif mod(nTrialsPerSession,size(designUnisVis,1))
            error('nTrialsPerSession is not divisible by number of unisVis conditions!');
        end
        % Packing all the trials into a final table (according to the session
        % order, randomized as well).
        allTrialsForExp = table;
        allTrialsBisAudPack = allTrialsBisAud;
        allTrialsBisVisPack = allTrialsBisVis;
        allTrialsUnisAudPack = allTrialsUnisAud;
        allTrialsUnisVisPack = allTrialsUnisVis;
        
        for i = 1:numel(sessions)
            if any(sessions(i) == [bis_a_L,bis_a_R])
                temp = cat(1,allTrialsBisAudPack(1:nTrialsPerSession,:));
                allTrialsBisAudPack(1:nTrialsPerSession,:) = [];
            elseif any(sessions(i) == [bis_v_L,bis_v_R])
                temp = cat(1,allTrialsBisVisPack(1:nTrialsPerSession,:));
                allTrialsBisVisPack(1:nTrialsPerSession,:) = [];
            elseif any(sessions(i) == [uni_a_L,uni_a_R])
                temp = cat(1,allTrialsUnisAudPack(1:nTrialsPerSession,:));
                allTrialsUnisAudPack(1:nTrialsPerSession,:) = [];
            elseif any(sessions(i) == [uni_v_L,uni_v_R])
                temp = cat(1,allTrialsUnisVisPack(1:nTrialsPerSession,:));
                allTrialsUnisVisPack(1:nTrialsPerSession,:) = [];
            end
            
            temp.session = repmat(i,size(temp.condition));
            temp.iTrialInSession = (randperm(size(temp.condition,1)))';
            if mod(sessions(i),2)
                % hand = 1 if left
                temp.hand = ones(size(temp.condition))*leftHand;
            else
                % hand = 2 if right
                temp.hand = ones(size(temp.condition))*rightHand;
            end
            allTrialsForExp = [allTrialsForExp;temp]; %#ok<AGROW>
        end
        
        allTrialsForExp = allTrialsForExp(:,[6:7,1:5,8]);
        allTrialsForExp = sortrows(allTrialsForExp,{'session','iTrialInSession'});
        
        % Saving the data
        save(trialFname,'allTrialsForExp','allTrialsBisAud','allTrialsBisVis',...
            'allTrialsUnisAud','allTrialsUnisVis','sessions');
    end
    
end

% Loading recorded sounds if applicable
if strcmp(audioMode,'rec')
    
    if strcmp(subID,'101') % for testing, that's me
        startID = 100;
    elseif strcmp(subID,'102') % Tao
        startID = 30;
    else
        error('There are no recordings for the subject ID');
    end
    
    REC.audioFreq = audioSamplingRate;
    REC.cutoff = '12000';
    REC.folderPath = DEC_2_setupdir(expStage,data_sound_sub);
    REC.recIDs = startID + (1:8);
    REC.stimDuration = stimulusDuration;
    REC.stimLocations = auditoryLocationLevels;
    REC.subID = subID;
    
    % Loading sounds form the sound database for the given subject to an
    % m x n cell array, where m = the number of locations, n = number of
    % recordings for the given location.
    REC.recordedSounds = getaudio(REC);
    
end


%%                       Loop for experiment runs
%%=========================================================================
% Setting the sessions for the current day
if strcmp(expMode,'train')
    sessionsForCurrDay = (1:nSessionsPerDay)+((iDay - 1)*nSessionsPerDay);
    nSessionsCurrDay = nSessionsPerDay;
else
    if isempty(sessionsForCurrDay)
        sessionsForCurrDay = (1:nSessionsPerDay)+((iDay - 1)*nSessionsPerDay);
        nSessionsCurrDay = nSessionsPerDay;
    else
        nSessionsCurrDay = size(sessionsForCurrDay,2);
    end
end

if sessionToStartWith > nSessionsCurrDay
   error('The number of the session to start with exceeds the number of sessions for the current day! ');
end

for iSessionInDay = sessionToStartWith:nSessionsCurrDay
    
    % The index number of the current session indexed linearly from the
    % beginning of the experiment. 
    currSession = sessionsForCurrDay(iSessionInDay);
    
    % Extracting the the trials for the current session.
    if syncMode
        locationAuditory = ones(nTrialsPerSession,1)*auditoryLocationLevels(3);
        locationVisual = ones(nTrialsPerSession,1)*auditoryLocationLevels(3);
        reliabilityVisual = ones(nTrialsPerSession,1)*highVisReliability;
        task = ones(nTrialsPerSession,1)*auditoryReport;
        condition = ones(nTrialsPerSession,1)*255;
        hand = ones(nTrialsPerSession,1)*leftHand;
    else
        locationAuditory = allTrialsForExp.locationAuditory(allTrialsForExp.session == currSession);
        locationVisual = allTrialsForExp.locationVisual(allTrialsForExp.session == currSession);
        reliabilityVisual = allTrialsForExp.reliabilityVisual(allTrialsForExp.session == currSession);
        task = allTrialsForExp.task(allTrialsForExp.session == currSession);
        condition = allTrialsForExp.condition(allTrialsForExp.session == currSession);
        hand = allTrialsForExp.hand(allTrialsForExp.session == currSession);
    end
    
    
    %                       Trial timing parameters
    %----------------------------------------------------------------------
    if syncMode
        % The amount of time between the onset of the fixation cross
        % and the stimulus onset (in seconds).
        t_toStimOnset = ones(nTrialsPerSession,1)*0.1;
        % Stimulus duration (in seconds).
        t_stimDuration = stimulusDuration;
        % The amount of time between the offset of the stimulus to the onset
        % of the response (in seconds).
        t_toResponse = 0.8;
        % Response cue duration
        t_respCueDuration = 0;
        % The amount of time available for answering. 
        t_endResponse = 0;
        % Feedback duration - if applicable (in seconds)
        t_feedbackDuration = 0;
    else
        % The amount of time between the onset of the fixation cross
        % and the stimulus onset (in seconds).
        t_toStimOnset = (rand(nTrialsPerSession,1)*0.5)+0.6;
        % Stimulus duration (in seconds).
        t_stimDuration = stimulusDuration;
        % The amount of time between the offset of the stimulus to the onset
        % of the response (in seconds).
        t_toResponse = 0.95;
        % Response cue duration
        t_respCueDuration = 0.1;
        % The amount of time available for answering. 
        t_endResponse = 2;
        % Feedback duration - if applicable (in seconds)
        t_feedbackDuration = 1.5;
    end
    
    %                                 ITD
    %----------------------------------------------------------------------
    if strcmp(audioMode,'itd')
        
        % Radius of the head in m. (average, it can potentially be
        % individualized)
        ITD.headRadius = 0.09;
        % The speed of sound in dry air at 20 ?C in m/s.
        ITD.soundSpeed = 343.2;
        % Length of the ramping up and down in the beginning and the end of
        % the sounds in seconds.
        ITD.ramp = 0.005;
        % Compute all possible stimulus locations.
        ITD.azimuthLevels = unique(locationAuditory(~isnan(locationAuditory)))';% +-10,+-3.3
        % Signed value of timeshift in seconds(ITD delay). Negative, if sounds are 
        % coming from the left.
        if ~isempty(ITD.azimuthLevels)
            ITD.timeShiftLevels = computeitd(ITD.azimuthLevels,'headRadius',...
                ITD.headRadius,'soundSpeed',ITD.soundSpeed);
        else
            ITD.timeShiftLevels = NaN;
        end
        
        % 3D array for storing the sounds for each trial.
        sounds = NaN(nTrialsPerSession,nAudioChannels,audioSamplingRate*t_stimDuration);
        
    %                           External speakers
    %----------------------------------------------------------------------
    elseif strcmp(audioMode,'speaker')
        % Length of the ramping up and down in the beginning and the end of
        % the sounds in seconds.
        SPEAKER.ramp = 0.005;
        
        % Mapping of the channels to physical locations (in degrees), i.e.
        % physical location of channel i is SPEAKER.channelMap(i)
        % Speakers are mapped physically as follows: 
        % Channel 1: close right
        % Channel 2: close left
        % Channel 3: far right
        % Channel 4: far left
        SPEAKER.channelMap = auditoryLocationLevels([3,2,4,1]);
        
        % Volume compensation for each channel to guarantee equal output
        % volume (must be measured externally). 
        SPEAKER.channelComp = [1.15,1,1.07,1];
        
        % The locations of the auditory cue
        SPEAKER.auditoryCueLocations = auditoryLocationLevels(2:3);
        
        % 3D array for storing the sounds for each trial.
        sounds = NaN(nTrialsPerSession,1,audioSamplingRate*t_stimDuration);
        
    %                           Audio recordings
    %----------------------------------------------------------------------    
    elseif strcmp(audioMode,'rec')
        
        % The number of how many times each set of recording is presented. 
        repNumSet = (nTrialsPerSession/numel(REC.stimLocations))/numel(REC.recIDs);
        if mod(repNumSet,1) ~= 0
            error('Each sound recording should be presented equal times! ');
        end
        
        % Generating a matrix 
        recIDsForSession = repmat(1:numel(REC.recIDs),numel(REC.stimLocations),repNumSet);
        for i = 1:size(recIDsForSession,1)
            recIDsForSession(i,:) = recIDsForSession(i,randperm(size(recIDsForSession,2)));
        end
        
        % This vector keeps track of how many times a given stimulus
        % location was used for generating stimuli. 
        numOfusedRecs = zeros(size(REC.stimLocations));
        
        % 3D array for storing the sounds for each trial. 
        sounds = NaN(nTrialsPerSession,nAudioChannels,size(REC.recordedSounds{1,1},1));
        
    end
    
    % 3D array for storing the spheres coordinates of the gaussian clouds 
    % for each trial.
    spheres = NaN(nTrialsPerSession,2,numOfSpheres);
    
    
    %       Saving the visual and auditory stimuli for each trial. 
    %----------------------------------------------------------------------
    for iTrial = 1:nTrialsPerSession
        
        if ~isnan(locationVisual(iTrial))
            % Generating spheres with desired xstd and ystd, centered 
            % around x = 0, y = 0. For each run a new set of gaussian
            % clouds are generated. 
            tempSpheres = generatespheres(reliabilityVisual(iTrial),...
                          vertSTDVisual,numOfSpheres,diameterSphere,0);
            % Shifting the generated shperes along the x axis to the
            % visual stimulus location. 
            tempSpheres(:,1) = tempSpheres(:,1) + locationVisual(iTrial);
            % Reshaping is required to fit the spheres in the destination
            % 3D matrix (dimensions: nTrialsPerRun x 2 x numOfSpheres).  
            spheres(iTrial,:,:) = shiftdim(tempSpheres',-1);
        end
        
        if ~isnan(locationAuditory(iTrial))
            
            if strcmp(audioMode,'itd')
                
                % Generate Gaussian white noise with mean = 0 and SD = 1/3 (so the
                % 3*SD < 1.
                noise = randn(round(audioSamplingRate * t_stimDuration),1)/3;
                % Cutoff the values which are more extreme than 3*SD.
                noise(noise > 1) = 1;
                noise(noise < -1) = -1;
                % Setting the ramping on's and off's at the beginning and the end
                % of the sounds.
                nSampleRamp = round(audioSamplingRate * ITD.ramp);
                % ramp-on
                noise(1:nSampleRamp) = noise(1:nSampleRamp).*(1:nSampleRamp)'/nSampleRamp;
                % ramp-off
                noise(end-nSampleRamp+1:end) = noise(end-nSampleRamp+1:end).*(nSampleRamp:-1:1)'/nSampleRamp;
                % Getting the appropriate amount of samples for time shift for the
                % particular stimulus location. Notice, that the number of samples
                % for the time shift is not signed anymore, as the the absolute
                % value of the
                nSampleTimeShift = round(audioSamplingRate * abs(ITD.timeShiftLevels(ITD.azimuthLevels == locationAuditory(iTrial))));
                % Shifting the appropriate channel with the given amount to
                % generate the ITD.
                if locationAuditory(iTrial) > 0
                    % In this case the sound should come from the right, hence it
                    % reaches the left ear later, so we shift the left channel with
                    % the needed amount.
                    tempSound(1,:) = [zeros(nSampleTimeShift,1); noise(1:end-nSampleTimeShift)];
                    tempSound(2,:) = noise;
                elseif locationAuditory(iTrial) < 0
                    tempSound(1,:) = noise;
                    tempSound(2,:) = [zeros(nSampleTimeShift,1); noise(1:end-nSampleTimeShift)];
                else
                    tempSound(1,:) = noise;
                    tempSound(2,:) = noise;
                end
                % Saving the generated sound.
                sounds(iTrial,:,:) = shiftdim(tempSound,-1);
                
            elseif strcmp(audioMode,'speaker')
                
                % Generate Gaussian white noise with mean = 0 and SD = 1/3 (so the
                % 3*SD < 1.
                noise = randn(round(audioSamplingRate * t_stimDuration),1)/3;
                % Cutoff the values which are more extreme than 3*SD.
                noise(noise > 1) = 1;
                noise(noise < -1) = -1;
                % Setting the ramping on's and off's at the beginning and the end
                % of the sounds.
                nSampleRamp = round(audioSamplingRate * SPEAKER.ramp);
                % ramp-on
                noise(1:nSampleRamp) = noise(1:nSampleRamp).*(1:nSampleRamp)'/nSampleRamp;
                % ramp-off
                noise(end-nSampleRamp+1:end) = noise(end-nSampleRamp+1:end).*(nSampleRamp:-1:1)'/nSampleRamp;
                                
                % Saving the generated sound.
                sounds(iTrial,:,:) = shiftdim(noise,-1);
                
            elseif strcmp(audioMode,'rec')
                
                % Find the index of the actual location.
                actLocIdx = find(REC.stimLocations == locationAuditory(iTrial));
                % Increment the corresponding numerator of the used locations.
                numOfusedRecs(actLocIdx) = numOfusedRecs(actLocIdx) + 1;
                % Find the number of recroding needed from the randomized pool.
                actRecIdx = recIDsForSession(actLocIdx,numOfusedRecs(actLocIdx));
                % Write the choosen sound to the sound array.
                sounds(iTrial,:,:) = shiftdim(REC.recordedSounds{actLocIdx,actRecIdx}',-1);
                
            end
            
        end
        
    end
    
    %              Saving information into the input struct. 
    %    This struct contains the necessary input data for one session.
    %----------------------------------------------------------------------
    % Audio parameters
    % Number of audio channels to use
    S.audioParam.nChannels = nAudioChannels;
    % Sampling rate
    S.audioParam.sr = audioSamplingRate;
    % Audio Mode
    S.audioParam.mode = audioMode;
    % Details
    if strcmp(audioMode,'itd')
        S.audioParam.details = ITD;
    elseif strcmp(audioMode,'speaker')
        S.audioParam.details = SPEAKER;
    elseif strcmp(audioMode,'rec')
        S.audioParam.details = REC;
    end
    % Auditory location levels
    S.auditoryLocationLevels = auditoryLocationLevels;
    % The auditory stimuli for the trials in the actual run. 
    S.auditoryStimuli = sounds;
    % Stimulus conditions for each trial in the actural run. 
    S.condition = condition;
    % Diameter of the spheres making up the gaussian clouds/visual stimuli. 
    S.diameterSphere = diameterSphere;
    % Experiment mode.
    S.expMode = expMode;
    % Codename of the experiment:
    S.expName = 'DEC_2';
    % Experiment stage
    S.expStage = expStage;
    % Eyetracker settings
    if testMode || syncMode
        S.eyelink.enable = false;
    else
        S.eyelink.enable = true;
    end
    S.eyelink.dummy = false;
    % Flag whether to give audio feedback
    if strcmp(expMode,'train') && iSessionInDay == 1
        % The first training session is auditory with feedback. 
        S.feedback = true;
    else
        S.feedback = false;
    end
    % Hand to use to answer
    S.hand = hand;
    % The number of the session.
    S.iDay = iDay;
    % The number of the actual run. 
    S.iSessionInDay = iSessionInDay;
    % The serial number of the session overall
    S.iSessionOverall = currSession;
    % Locations of auditory stimuli in the actual run.
    S.locationAuditory = locationAuditory;
    % Locations of visual stiumli in the actual run. 
    S.locationVisual = locationVisual;
    % Maximum visual angle
    S.maxVisAngle = maxVisAngle;
    % The number of trials in the actual run.
    S.nTrials = nTrialsPerSession;
    % The reliabilities of the stimuli in the actual run. 
    S.reliabilityVisual = reliabilityVisual;
    % The name of the RunScript
    S.runScriptName = runScriptName;
    % The resolution of the screen. 
    S.screenRes = reqScreenRes;
    % The order of sessions in the whole experiment. 
    if syncMode
        S.sessionOrder = NaN;
    else
        S.sessionOrder = sessions;
    end
    % The name of the computer on which the experiment was run. 
    S.setupSpec = setupSpec;
    % The actual subject's ID number.
    S.subID = subID;
    % Flag for synchronization measurement mode
    S.syncMode = syncMode;
    % Task
    S.task = task;
    % Flag for testing mode
    S.testMode = testMode;
    % Timing parameters.
    S.timing.t_endResponse = t_endResponse;
    S.timing.t_feedbackDuration = t_feedbackDuration;
    S.timing.t_respCueDuration = t_respCueDuration;
    S.timing.t_stimDuration = t_stimDuration;
    S.timing.t_toResponse = t_toResponse;
    S.timing.t_toStimOnset = t_toStimOnset;
    % The experimenter's username. 
    S.user = getenv('username');
    % Viewing distance in mm. 
    S.viewDist = viewDist;
    % Visual location levels
    S.visualLocationLevels = visualLocationLevels;
    % The visual stimuli in the actual run. These are the positions of the 
    % spheres of the gaussian clouds.  
    S.visualStimuli = spheres;
    
    fprintf('\n \n Day: %i \n Session: %i \n Press any key to go on \n Press esc  to quit \n \n',iDay,iSessionInDay)
    KbWait;
    
    [~,~,keyCode] = KbCheck;
    if keyCode(quit)
        break;
    end
    
    info = DEC_2_Present(S);
    
    % Saving experiment data
    savedfname = fullfile(DEC_2_setupdir(expStage,'data_behav_sub',subID),...
        [S.subID,'_',S.expMode,'_',num2str(S.iDay),'_',num2str(S.iSessionInDay),...
        '_','all_runscript','_',datestr(now,'ddmmyyyy_HHMM'),'.mat']);
    save(savedfname);
    
    quickanalyse(info);
    
    % Clearing the java heap memory to avoid java OutOfMemory exception
    jheapcl;
    
end

% Removing MatlabGarbageCollector.jar from the dynamic java path
javaaddpath(fullfile(DEC_2_setupdir(expStage,'pres'),'MatlabGarbageCollector.jar'));
