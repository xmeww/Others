function varargout = runeyelink(mode,varargin)
% Performs various operations with the eyelink system. 
% 
% USAGE: 
%   S = runeyelink('setup',S);
%       runeyelink('record',S);
%       runeyelink('close',S);
% 
% INPUT: 
%   mode: argument to set the mode in which to operate the function.
%       Possible values are: 'setup','record','close'. 
%   S: structure with all necessary parameters from the experiment   
%      contains a logical, whether or not to use dummy mode for 
%      initialization
%
% OUTPUT: 
%   S: with additional values in it. 

%% Parsing input
p = inputParser;

validModes = {'close','record','setup'};
checkMode = @(x) any(validatestring(x,validModes));

addRequired(p,'mode',checkMode);
addRequired(p,'S');

parse(p,mode,varargin{:});

mode = p.Results.mode;
S = p.Results.S;

if any(strcmp(mode,{'close','setup'})) && isempty(S)
    error('''S'' must be specified for modes ''close'' and ''setup''! ');
end

    
%% 
switch mode
    
    case 'setup'
        
        varargout{1} = setupeyelink(S);
        
    case 'record'
        
        recordeyelink(S);
        
    case 'close'
        
        closeeyelink(S);
        
end


end

function S = setupeyelink(S)
% Sets up the Eyelink eyetracker. 
% 
% USAGE: 
%   S = setupeyelink(S)
% 
% INPUT: 
%   S: structure with all necessary parameters from the experiment   
%      contains a logical, whether or not to use dummy mode for 
%      initialization
%
% OUTPUT: 
%   S.eyelink.el:
%       structure containing specifications, and used by eyelink
%       functions
% 
% DETAILS: 
%   Tries to connect to an EyeLink eyetracker. Performs the initialization 
%   of eyelink, connection to the eyelink data file, as well as calibration
%   and drift correction. If this works out the eye tracker can be used 
%   during stimulus presentation. 
% 
%   Sample data:
%   Keyword     Data Type
%   LEFT,       Sets the intended tracking eye (usually include both LEFT and...
%   RIGHT       RIGHT)
%   GAZE        includes screen gaze position data
%   GAZERES     includes units-per-degree screen resolution at point of gaze
%   HREF        head-referenced eye position data 
%   HTARGET     target distance and X/Y position (EyeLink Remote only)
%   PUPIL       raw pupil coordinates 
%   AREA        pupil size data (diameter or area)
%   BUTTON      buttons 1-8
%   STATUS      state and change flags warning and error flags (not yet supported)
%   INPUT       input port data lines (not yet supported)
% 
%   Event data:
%   Keyword     Effect
%   GAZE        includes display (gaze) position data.
%   GAZERES     includes units-per-degree screen resolution (for start, end of event)
%   HREF        includes head-referenced eye position
%   AREA        includes pupil area or diameter
%   VELOCITY    includes velocity of parsed position-type (average, peak, start and end)
%   STATUS      includes warning and error flags, aggregated across event (not yet supported)
%   FIXAVG      include ONLY averages in fixation end events, to reduce file size
%   NOSTART     start events have no data other than timestamp
% 
%   GAZE,GAZERES,AREA,HREF,VELOCITY  - default: all useful data
%   GAZE,GAZERES,AREA,FIXAVG,NOSTART - reduced data for fixations
%   GAZE,AREA,FIXAVG,NOSTART         - minimal data
%   

%% Parsing input parameters
p = inputParser;

fieldsReqS = {'disp','expMode','expName','expStage','eyelink','iDay',...
    'iSessionInDay','iSessionOverall','subID'}; 
fieldsReqDisp = {'grey','white','win'};
fieldsReqEyelink = {'dispRect','dummy','enable'};

checkS = @(x) all(isfield(x,fieldsReqS)) && ...
              all(isfield(x.disp,fieldsReqDisp)) && ...
              all(isfield(x.eyelink,fieldsReqEyelink));

addRequired(p,'S',checkS);

parse(p,S);

S = p.Results.S;

%%

% Do we even want to record eyetracking data?
if ~S.eyelink.enable
    warning('Eyelink disabled! ');
    S.eyelink.online = 0;
    return;
end

% Initialization of the connection with the eyetracker.
S.eyelink.online = EyelinkInit_Mate(S.eyelink.dummy, 1);
if ~S.eyelink.online
    warning('Eyelink initialization aborted! ');
    cleanup;
    return;
end

% Setting up the default parameters for the eyelink
S.eyelink.el = EyelinkInitDefaults(S.disp.win);

% Customizing parameters
S.eyelink.el.backgroundcolour        = S.disp.grey;
S.eyelink.el.msgfontcolour           = S.disp.white;
S.eyelink.el.imgtitlecolour          = S.disp.white;
S.eyelink.el.calibrationtargetcolour = S.disp.white;
S.eyelink.el.calibrationtargetsize   = 1;
S.eyelink.el.calibrationtargetwidth  = 0.5;
S.eyelink.el.targetbeep              = 0;
S.eyelink.el.feedbackbeep            = 0;
S.eyelink.el.displayCalResults       = 1;
S.eyelink.el.eye_used                = S.eyelink.el.RIGHT_EYE;
% s.el.allowlocaltrigger       = 0; % allow user to trigger him or herself
% s.el.allowlocalcontrol       = 0; % allow control from subject-computer

EyelinkUpdateDefaults(S.eyelink.el);

% Open file on host computer
if Eyelink('Openfile','data.edf')
    warning('Cannot create EDF file ''%s'' ','data.edf');
    cleanup;
    return;
end
Eyelink('Command','add_file_preamble_text = "Experiment %s, %s, subject %s, day %d, ''%s'' session %d',...
    S.expName,S.expStage,S.subID,S.iDay,S.expMode,S.iSessionInDay);


% make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~S.eyelink.dummy
    warning('Eyelink is not available anymore! ');
    cleanup;
    return;
end

% This Command is crucial to map the gaze positions from the tracker to
% screen pixel positions to determine fixation
Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld',...
    S.eyelink.dispRect(1),S.eyelink.dispRect(2),S.eyelink.dispRect(3)-1,S.eyelink.dispRect(4)-1);
Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld',...
    S.eyelink.dispRect(1),S.eyelink.dispRect(2),S.eyelink.dispRect(3)-1,S.eyelink.dispRect(4)-1);

% Use conservative online saccade detection (cognitive setting)
% Equivalent with these settings:
% Eyelink('Command','recording_parse_type = GAZE');
% Eyelink('Command','saccade_velocity_threshold = 30');
% Eyelink('Command','saccade_acceleration_threshold = 9500');
% Eyelink('Command','saccade_motion_threshold = 0.15');
% Eyelink('Command','saccade_pursuit_fixup = 60');
% Eyelink('Command','fixation_update_interval = 0');
Eyelink('Command','select_parser_configuration = 0');

% Other tracker configurations
Eyelink('Command','calibration_type = HV9');
Eyelink('Command','generate_default_targets = YES');
Eyelink('Command','enable_automatic_calibration = YES');
Eyelink('Command','automatic_calibration_pacing = 1000');
Eyelink('Command','binocular_enabled = NO');
Eyelink('Command','use_ellipse_fitter = NO');
Eyelink('Command','sample_rate = 2000');
% illumination, 1 = 100%, 2 = 75%, 3 = 50%
Eyelink('Command','elcl_tt_power = %d', 3);

% Automatic sequencing during calibration?
[~,reply] = Eyelink('ReadFromTracker','enable_automatic_calibration');
if reply 
    fprintf('Automatic sequencing ON\n');
else
    fprintf('Automatic sequencing OFF\n');
end

% set edf data
Eyelink('Command','file_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,STATUS,INPUT');
Eyelink('Command','file_event_data  = GAZE,GAZERES,AREA,HREF,VELOCITY');
Eyelink('Command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');

% % set link data (can be used to react to events online)
% Eyelink('Command','link_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,FIXUPDATE,INPUT');
% Eyelink('Command','link_sample_data  = LEFT,GAZE,GAZERES,AREA,STATUS,INPUT');

% Calibrate the eye tracker
EyelinkDoTrackerSetup_Mate(S.eyelink.el);

% Cleanup routine:
function cleanup
    % Shutdown Eyelink:
    Eyelink('Shutdown');
    S.eyelink.online = 0;
end

end

function recordeyelink(S)
%% Parsing input parameters
p = inputParser;

checkS = @(x) isfield(x.eyelink,'online');

addRequired(p,'S',checkS);

parse(p,S);

S = p.Results.S;

%%
if S.eyelink.online
    status = Eyelink('StartRecording');
    if status ~= 0
        error('Eyelink StartRecording error');
    end
end

end

function closeeyelink(S)
%% Parsing input parameters
p = inputParser;

checkS = @(x) all(isfield(x.eyelink,{'edfFileName','online'})) && ischar(x.eyelink.edfFileName);

addRequired(p,'S',checkS);

parse(p,S);

S = p.Results.S;

%%

if S.eyelink.online
    
    % Eyelink file from input
    [~,name,ext] = fileparts(S.eyelink.edfFileName);
    
    % Stop recording and close file
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    
    % Receive file from host computer and move to the experimental data folder
    try
        fprintf('Receiving data file ''%s''\n',[name,ext]);
        status = Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n',status);
        end
    catch
        warning('Problem receiving data file ''%s''',[name,ext]);
    end
    if exist('data.edf','file')
        movefile('data.edf',S.eyelink.edfFileName); % rename and move file to destination
        fprintf('Data file moved to subject folder.\n');
    end
    
    % Converting edf file to asc
    if exist(S.eyelink.edfFileName,'file')
        try
            fprintf('Converting file %s\n',[name,ext]);
            [status,~] = system(sprintf('edf2asc %s',S.eyelink.edfFileName));
            
%             if status
%                 warning('Can not convert file %s ',[name,ext]);
%             end
        catch
            warning('Can not convert file %s ',[name,ext]);
        end
    end
    
    % Shut down the eyetracker
    Eyelink('Shutdown');
        
end

end

