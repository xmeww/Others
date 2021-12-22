clear all;
clc;
% error beeps off
beep off


%%                         TRIGGER
% =========================================================================
% Important! Open IOPort for sending Triggers
USBTTL('Open', 'usbserial-ftEHGLIF');



%%                       Subject Details
% =========================================================================
S.SubID       = input('input the subjects ID: ', 's');
S.SubAge      = input('input the subjects age: ');
S.handedness  = input('input the subjects handedness: ', 's');
S.proffession = input('input the subjects proffession: ', 's');
S.proffession = input('Any hearing disabilities? ', 's');



%%                       Saving
% =========================================================================
folder        = input('folderName: ', 's');
sdir          = fullfile(['/Users/anette/Documents/MATLAB/Experiment2/Subdata/' folder]);
% do a sanity check that folder exists and no data are in there that could
% be overwritten
filesInFolder = dir(sdir);
if isdir(sdir) == 1
    fprintf([sdir '\n  \n'])
else
    error('!!!directory not found!!!')
end

if length(filesInFolder)>2
    warning('MATLAB:rmpath:DirNotEmpty', ' \n \n !!!directory not empty!!! \n ')
    filesInFolder.name
    goOn = input('\n \n directory is not empty - are you sure to continue?: 1 - yes, 2 - no \n');
    if goOn == 2
        error('!!!directory not empty!!!')
    end
end
S.saveDir = sdir;



%%                      start program
% =========================================================================

% -------------------------------------------------------------------------
% GENERAL
% number of sessions
% 16 sessions: 8 high noise vs 8 low noise
nSessions = 5;
% always needs to be put if using KbName!!!!!
KbName('UnifyKeyNames');
% define names of keys
quit = KbName('q');    
S.debug = 0;

% -------------------------------------------------------------------------
% determine loudness per frequency
load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/Audiometry.mat'])
dBdiff    = 2;
% loudness of masking frequencies
audThresM = abs(dB(2,:))/dBdiff;
for iTimes = 1:10
    S.audScaling(iTimes) = (1/(10^(dBdiff/20)))^audThresM(iTimes);
end

% loudness of target Frequencies
audThresT = abs(dBTar(2,:))/dBdiff;
for iTimes = 1:5
    S.audScalingTar(iTimes) = (1/(10^(dBdiff/20)))^audThresT(iTimes);
end

% loudness masking frequencies should be presented at
loudM    = 54; % loudness in dB
audLoudM = loudM/dBdiff;
S.loudDB = (10^(dBdiff/20))^audLoudM;

% loudness target frequencies should be presented at
loudT       = 52;
audLoudT    = loudT/dBdiff;
S.loudDBTar = (10^(dBdiff/20))^audLoudT;


% -------------------------------------------------------------------------
% AUDITORY
% percent aud signal strength
S.SNRaud = 1/30;      

% percent aud noise strength
S.SNRaudNoise = 1/30; 


% -------------------------------------------------------------------------
% VISUAL
% percent vis signal strength
S.SNRvis = 0.05; 

% 1 = sine-wave grating, 2 = gabor grating, 3 = square
S.visStim = 1;  

% defines whether the visual stimulus wil be presented with the first
% (=1), second (=2) auditory target or both (= 3)
S.timeAV = 3;


% -------------------------------------------------------------------------
% TRIAL & TARGETS
% nr of Tragets
S.nrTrials = 240;  

% duration of one stim in seconds
S.durTarAud = 0.3;   

% duration (ms) of the visual target
S.durTarVis = 0.15;
        
% play masking sound or not
S.training = 0;    

% mean SOA between targets
S.SOA = 1.05; % 1050 ms

% Which target frequencies
S.nrFreqs = 1:5;


% -------------------------------------------------------------------------
% MEG
% MEG recording = 1
S.MEG = 1;

% play mask
S.mask = 1;

% -------------------------------------------------------------------------
% start program
for iSession = 1:nSessions
 
    S.Session = iSession;
    
    if iSession == 5
        S.mask = 0;
    end
    
    fprintf('\n press any button to go on to session: %i \n press /q/ to quit \n \n', iSession)
    KbWait;
    % stopping the loop if user presses 'q'
    [~,~,keyCode] = KbCheck;
    if keyCode(quit)
        break;
    end
    
    % start program
    o = Exp2_IM_10(S);
    savedfname=[S.saveDir, '/STARTSubject_' S.SubID,'_', num2str(iSession), '.mat'];
    save(savedfname)
    
    Start_Exp2_IM_Analyse;
   
end

USBTTL('Close')
fprintf('\n \n Done :) \n \n ')
% Done :)