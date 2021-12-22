function outPath = DEC_2_setupdir(expStage,dirID,varargin)
% Set up folders for the experiment. 
% 
% USAGE: 
%   out = DEC_2_setupdir(expStage,dirID);
%   out = DEC_2_setupdir(expStage,dirID,'Name',Value);
%
% INPUT:
%   dirID: 
%   Name-Value input arguments
%       setupID:
%       subID:
% OUTPUT: 
%   outPath: path of the required directory
%
% SIDE EFFECTS:
%   Some folders are created if not found. 
%

%% Parsing input. 
p = inputParser;

validExpStages = {'pilot_1','pilot_2','pilot_3','pilot_4','final'};

validDirIDs = {'anal_behav','anal_eeg','anal_behav_scripts','anal_behav_sub',...
    'anal_eeg_scripts','anal_eeg_sub','anal_eeg_sub_mvpa','anal_eeg_sub_mvpa_preproc',...
    'anal_eeg_sub_erp','anal_eeg_sub_tf','data_behav','data_eeg','data_sound',...
    'data_behav_sub','data_eeg_sub','data_sound_sub','pres','study_root',...
    'utils'};

checkExpStage = @(x) any(validatestring(x,validExpStages));
checkDirID = @(x) any(validatestring(x,validDirIDs));

addRequired(p,'expStage',checkExpStage);
addRequired(p,'dirID',checkDirID);
addOptional(p,'subID','',@ischar);

parse(p,expStage,dirID,varargin{:});

expStage = p.Results.expStage;
dirID = p.Results.dirID;
subID = p.Results.subID;

%% Setting up the basic directories if necessary. 
[~,setupID] = system('hostname');
setupID = regexp(setupID,'\w*-\w*','match');
setupID = setupID{1};

if strcmp(computer,'PCWIN64') || strcmp(computer,'PCWIN')
    [~,userID] = system('echo %username%');
    userID = regexp(userID,'\w*','match');
    userID = userID{1};
elseif strcmp(computer,'GLNXA64') || strcmp(computer,'GLNX86')
    [~,userID] = system('id -n -u');
    userID = regexp(userID,'\w*','match');
    userID = userID{1};
else
    error('Can''t find user name!');
end

% base folder depending on the setup
if strcmpi(setupID,'PSYCHL-132433')
    baseDir = 'E:';
    mode = 'home';
elseif strcmpi(setupID,'COLLES-151401')
    baseDir = fullfile('/home',userID);
    mode = 'analysis';
elseif strcmpi(setupID,'COLLES-d152225')
    baseDir = fullfile('/home',userID);
    mode = 'analysis';
elseif strcmpi(setupID,'COLLES-140591')
    baseDir = fullfile('C:','Users',userID);
    mode = 'presentation';
else
    error('Unidentified setup!')
end

mypath.study_root = fullfile(baseDir,'Studies','VE_EEG_Decode_BayesianIntegration');
if ~exist(mypath.study_root,'dir')
    mkdir(mypath.study_root);
end

mypath.anal_behav = fullfile(mypath.study_root,'behavioural_analysis',expStage);
if ~exist(mypath.anal_behav,'dir') && strcmp(mode,'home')
    mkdir(mypath.anal_behav);
end

mypath.anal_eeg = fullfile(mypath.study_root,'EEG_analysis',expStage);
if ~exist(mypath.anal_eeg,'dir') && any(strcmp(mode,{'home','analysis'}))
    mkdir(mypath.anal_eeg);
end

mypath.data_behav = fullfile(mypath.study_root,'behavioural_data',expStage);
if ~exist(mypath.data_behav,'dir') && any(strcmp(mode,{'home','presentation'}))
    mkdir(mypath.data_behav);
end

mypath.data_eeg = fullfile(mypath.study_root,'EEG_data',expStage);
if ~exist(mypath.data_eeg,'dir') && any(strcmp(mode,{'home','presentation'}))
    mkdir(mypath.data_eeg);
end

mypath.data_sound = fullfile(mypath.study_root,'sound_database',expStage);
if ~exist(mypath.data_sound,'dir') && any(strcmp(mode,{'home','presentation'}))
    mkdir(mypath.data_sound);
end

mypath.pres = fullfile(mypath.study_root,'presentation',expStage);
if ~exist(mypath.pres,'dir') && any(strcmp(mode,{'home','presentation'}))
    mkdir(mypath.pres);
end

%% Returning the required path
% If not found either an error is thrown or the folder is created. 

if strcmp(dirID,'anal_behav')
    outPath = mypath.anal_behav;
elseif strcmp(dirID,'anal_eeg')
    outPath = mypath.anal_eeg;
elseif strcmp(dirID,'anal_behav_scripts')
    outPath = fullfile(mypath.anal_behav,'scripts');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_behav_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_behav,subID);
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_scripts')
    outPath = fullfile(mypath.anal_eeg,'scripts');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_eeg,subID);
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_sub_mvpa')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_eeg,subID,'MVPA');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_sub_mvpa_preproc')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_eeg,subID,'MVPA','preproc');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_sub_erp')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_eeg,subID,'ERP');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'anal_eeg_sub_tf')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.anal_eeg,subID,'TF');
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'data_behav')
    outPath = mypath.data_behav;
elseif strcmp(dirID,'data_eeg')
    outPath = mypath.data_eeg;
elseif strcmp(dirID,'data_behav_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.data_behav,subID);
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'data_eeg_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.data_eeg,subID);
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
elseif strcmp(dirID,'data_sound_sub')
    if strcmp(subID,'')
        error('Subject ID must be specified!');
    end
    outPath = fullfile(mypath.data_sound,subID);
    if ~exist(outPath,'dir')
        error('The specified folder does not exist!');
    end
elseif strcmp(dirID,'pres')
    outPath = mypath.pres;
elseif strcmp(dirID,'study_root')
    outPath = mypath.study_root;
elseif strcmp(dirID,'utils')
    outPath = fullfile(baseDir,'MATLAB');
    if ~exist(outPath,'dir')
        error('The specified folder does not exist!');
    end
else
    outPath = '';
end


end