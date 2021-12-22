function quickanalyse(info)
% Quick analysis for assessing data quality of a given run. 

%% Parsing input
p = inputParser;

addRequired(p,'info',@isstruct);

parse(p,info);

info = p.Results.info;


%% Preprocessing data
% Loading files specifying parameters
trigDef = load('trigDef.mat');
trigDef = trigDef.trigDef;
setupSpec = load('setup_spec.mat');
setupSpec = setupSpec.setup_spec;
trialRejCrit = load('trialRejCrit.mat');
trialRejCrit = trialRejCrit.trialRejCrit;

% Preprocessing behavioural data
behavData = preprocbehav(info);
if isempty(behavData)
    warning('Too few responses! ');
    return;
end
% Preprocessing eyetracking data
% [~,fileName] = fileparts(info.inputParameter.eyelink.edfFileName);
% ascFilePath = fullfile(DEC_2_setupdir(info.inputParameter.expStage,'data_behav_sub',info.subID),[fileName,'.asc']);
ascFilePath = strrep(info.inputParameter.eyelink.edfFileName,'.edf','.asc');
eyeData = preproceye(ascFilePath,trigDef,setupSpec);

%% Checking trials

saccadeThresholds = [trialRejCrit.saccade.ampl,trialRejCrit.saccade.pv,trialRejCrit.saccade.dur];

% Defining fixation criteria
fixation = struct();
% x coordinate of fixation cross in pixels
fixation.x = (setupSpec.screen_coords_pix(3)-setupSpec.screen_coords_pix(1))/2;
% y coordinate of fixation cross in pixels
fixation.y = (setupSpec.screen_coords_pix(4)-setupSpec.screen_coords_pix(2))/2;
% radius of the allowed area around the fixation cross in pixels.
% (converted from 2 degrees)
fixation.r = round(posdeg2pix(trialRejCrit.fixationRadiusDeg,...
    setupSpec.view_dist,setupSpec.screen_coords_pix(3)/setupSpec.screen_coords_mm(3)));
fixationCriteria = [fixation.x,fixation.y,fixation.r];

% Marking bad trials
% Last column of trialinfo gives the following info:
% 0 - good
% 1 - no response
% 2 - early response
% 3 - wrong hand
% 5 - missing eyetracker data
% 6 - eyeblink
% 7 - saccade
% 8 - wrong fixation location

badTrials = checktrials(behavData.data,eyeData,trigDef,trialRejCrit.winOfInt,trialRejCrit.minRespTime,saccadeThresholds,fixationCriteria);

%% Plotting results
plotbehav(behavData,info);

plotbadtrials(badTrials,info);

close all;

end


function dataBehav = preprocbehav(info)

% Collecting relevant data from the run. 
session = repmat(info.iSessionInDay,info.inputParameter.nTrials,1);
iTrialInSession = (1:info.inputParameter.nTrials)';
condition = info.inputParameter.condition;
task = info.inputParameter.task;
locVis = info.inputParameter.locationVisual;
locAud = info.inputParameter.locationAuditory;
relVis = info.inputParameter.reliabilityVisual;
resp = info.resp.locResp;
respTime = info.resp.locRespTimes - info.timing.actAuditoryStartTime;
wrongHand = info.resp.wrongHand;

if sum(~isnan(resp)) < size(resp,1)/2
    dataBehav = [];
    return;
end

dataBehav.data = table(session,iTrialInSession,condition,task,locAud,locVis,relVis,resp,respTime,wrongHand);

dataBehav.subID = info.subID;

end


function dataEye = preproceye(ascFilePath,trigDef,setupSpec)

stimTriggers = trigDef.trig_eye(trigDef.type == 'stim');
trig_visonset_corr_eyelink = setupSpec.trig_visonset_corr_eyelink;

dataEye = struct('event',[],'Fs',[]);

[event,hdr] = read_eyelink_event(ascFilePath);

evValues = {event.value}';
evTypes = {event.type}';
evStartSamples = [event.sample]';
% Finding the onset samples of stimulus triggers
targEvStartSamples = evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,stimTriggers));
% Correcting stimulus triggers for trigger-visual onset asynchrony
targEvStartSamples = targEvStartSamples+(trig_visonset_corr_eyelink*hdr.Fs);
% Saving corrected values into the original structure
evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,stimTriggers)) = targEvStartSamples;
evStartSamples = num2cell(evStartSamples);
[event.sample] = evStartSamples{:};

dataEye.event = event;
dataEye.Fs = hdr.Fs;

end


function badTrials = checktrials(behavData,eyeData,trigDef,winOfInt,minRespTime,saccadeThresholds,fixationCriteria)

nTrials = size(behavData,1);
% Array for collecting rejection info
badTrials = zeros(nTrials,1);

for actTrialInSession = 1:nTrials
    %% Checking response
    % Is there a response? 
    if isnan(behavData.resp(behavData.iTrialInSession == actTrialInSession))
        badTrials(actTrialInSession) = 1;
        continue;
    end
    
    % Is the response too early?
    if behavData.respTime(behavData.iTrialInSession == actTrialInSession) < minRespTime
        badTrials(actTrialInSession) = 2;
        continue;
    end
    
    % Was the response made with the wrong hand?
    if behavData.wrongHand(behavData.iTrialInSession == actTrialInSession)
        badTrials(actTrialInSession) = 3;
        continue;
    end
    
    %% Checking eye events
    % Finding eye events corresponding the current trial
    eyeEventsActTrial = selecteyeevents(eyeData,actTrialInSession,winOfInt,trigDef);
    if isempty(eyeEventsActTrial)
%         warning('No eyetracking data for trial %d',actTrialInSession);
        badTrials(actTrialInSession) = 5;
        continue;
    end
    
    % Is there a blink? 
    if any(ismember({eyeEventsActTrial.type},'Blink'))
        badTrials(actTrialInSession) = 6;
        continue;
    end
    
    % Is there a saccade? 
    if any(ismember({eyeEventsActTrial.type},'Saccade'))
        
%         badTrials(actTrialInSession) = 7;
%         continue;
        
        saccades = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Saccade'));
        % Are any of the saccades above threshold?
        if checksaccades(saccades,saccadeThresholds,eyeData.Fs)
            badTrials(actTrialInSession) = 7;
            continue;
        end
        
    end
    
    % Is fixation average location correct?
    if all(~ismember({eyeEventsActTrial.type},'Fixation'))
        badTrials(actTrialInSession) = 8;
        continue;
    else
        fixations = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Fixation'));
        % Are any of the saccades above threshold?
        if checkfixations(fixations,fixationCriteria)
            badTrials(actTrialInSession) = 8;
        end
        
    end
    
end


end


function selectedEvents = selecteyeevents(eyeData,actTrial,winOfInt,trigDef)
% Selects eyetracker events corresponding to a particular trial
% 
% INPUT: 
%   eyeData: structure array of eyetracker data of the session
%       corresponding to the trial of interest
%   actTrial: serial number of the trial of interest within its session
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   trigDef: trigger definition table
% 
% OUTPUT:
%   selectedEvents: structure array of events corresponding to the trial of
%       interest. Empty structure if no events were found
% 

evValues = {eyeData.event.value}';
evTypes = {eyeData.event.type}';
evStartSamples = [eyeData.event.sample]';
evEndSamples = evStartSamples+[eyeData.event.duration]';
targEvStartSamples = evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,trigDef.trig_eye(trigDef.type == 'stim')));

% If the trial is not found return an empty structure
if ~ismember(1:size(targEvStartSamples,1),actTrial)
    selectedEvents = struct([]);
    return;
end

% Find window of interest around target event
startSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(1)*eyeData.Fs);
endSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(2)*eyeData.Fs);

% Find events which:
% begin within the window of interest
crit1 = evStartSamples >= startSampleWinOfInt & evStartSamples <= endSampleWinOfInt;
% end within the winow of interest
crit2 = evEndSamples >= startSampleWinOfInt & evEndSamples <= endSampleWinOfInt;
% cover the whole window of interest
crit3 = evStartSamples <= startSampleWinOfInt & evEndSamples >= endSampleWinOfInt;
    
% Choose those events which fall under any of the above criteria
selectedEvents = eyeData.event(crit1 | crit2 | crit3);

end


function aboveThreshold = checksaccades(saccades,thresholds,Fs)
% Checks whether saccade parameters are above the given threshold
% 
% INPUT: 
%   saccades: structure array containing saccade events
%   thresholds: 1 x 3 vector of the threshold values for defining a saccade
%               thresholds(1) = minimum amplitude in degrees
%               thresholds(2) = minimum peak velocity in degrees/seconds
%               thresholds(3) = minimum duration in seconds
%   Fs: sampling frequency
% 
% OUTPUT:
%   aboveThreshold: true if any of the saccade events are above all three 
%       thresholds. 
% 

if any(~ismember({saccades.type},'Saccade'))
    error('Event of wrong type passed to checksaccades!');
end

temp = regexp({saccades.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
ampl = temp(:,1);
pv = temp(:,2);
% Converting duration from samples to seconds
dur = [saccades.duration]'./Fs;

if ampl >= thresholds(1) % & pv >= thresholds(2) & dur >= thresholds(3))
    aboveThreshold = true;
else
    aboveThreshold = false;
end

end


function failedCriteria = checkfixations(fixations,criteria)
% Checks whether fixation parameters meet the given criteria
% 
% INPUT: 
%   fixations: structure array containing fixation events
%   criteria: vector of the criteria for accepting a fixation event
%               criteria(1) = fixation cross x coordinate (in pixels)
%               criteria(2) = fixation cross y coordinate (in pixels)
%               criteria(3) = radius of the accepted area around the
%                             fixation cross (in pixels)
% 
% OUTPUT:
%   failedCriteria: true if any of the fixation events fail to meet the
%       criteria
% 

if any(~ismember({fixations.type},'Fixation'))
    error('Event of wrong type passed to checkfixations!');
end

temp = regexp({fixations.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
% This should produce an N x 2 matrix, where N is the number of fixation
% events and the two columns are the average x and y positions 
% respectively. 

if any(~isPointInCircle(temp,repmat(criteria,size(temp,1),1)))
    failedCriteria = true;
else
    failedCriteria = false;
end

end


function plotbehav(behavData,info)

if unique(behavData.data.task) == 1 % Auditory localization
    
    if all(isnan(behavData.data.locVis))
        behavData.accuracies = varfun(@calcaccuracy,behavData.data(~isnan(behavData.data.resp),:),...
            'InputVariables','resp','GroupingVariables',{'locAud'});
        behavData.unisens = true;
        behavData.figTitle = 'Auditory localization, Unisensory stimuli';
    else
        behavData.accuracies = varfun(@calcaccuracy,behavData.data(~isnan(behavData.data.resp),:),...
            'InputVariables','resp','GroupingVariables',{'locVis','locAud','relVis'});
        behavData.unisens = false;
        behavData.figTitle = 'Auditory localization, Bisensory stimuli';
    end
    
else % Visual localization
    
    if all(isnan(behavData.data.locAud))
        behavData.accuracies = varfun(@calcaccuracy,behavData.data(~isnan(behavData.data.resp),:),...
            'InputVariables','resp','GroupingVariables',{'locVis','relVis'});
        behavData.unisens = true;
        behavData.figTitle = 'Visual localization, Unisensory stimuli';
    else
        behavData.accuracies = varfun(@calcaccuracy,behavData.data(~isnan(behavData.data.resp),:),...
            'InputVariables','resp','GroupingVariables',{'locAud','locVis','relVis'});
        behavData.unisens = false;
        behavData.figTitle = 'Visual localization, Bisensory stimuli';
    end
    
end


% Drawing a figure with the relevant info. 
if behavData.unisens
    
    % Unisensory visual
    if unique(behavData.data.task) == 2 
        
        figure();
        set(gcf,'Units','normalized','Position',[0.05 0.25  0.5 0.25],...
            'PaperPositionMode','auto');
        
        taskLevels = unique(behavData.accuracies.locVis);
        taskVarName = 'locVis';
                
        iFig = 1;
        
        for i = 1:length(taskLevels)
            
            acc = behavData.accuracies.calcaccuracy_resp(behavData.accuracies.(taskVarName) == taskLevels(i),:)';
            
            subplot(1,numel(taskLevels),iFig);
            plot(repmat(taskLevels,1,size(acc,2)),acc,'LineWidth',1.5,...
                'Marker','o');
            ylim([0 1]);
            if i == 1
                title(sprintf('\n'));
                legend({'R+','R-'},'Location','NorthEast');
            else
                set(gca,'YTickLabel',[]);
            end
            set(gca,'XTick',taskLevels);
            
            iFig = iFig + 1;
            
        end
        
        suplabel(behavData.figTitle,'t',[.08 .08 .84 .84]);
        
        savedfname = fullfile(DEC_2_setupdir(info.inputParameter.expStage,'data_behav_sub',info.subID),...
            [info.subID,'_',info.inputParameter.expMode,'_',num2str(info.iDay),'_',...
            num2str(info.iSessionInDay),'_','qanalyse_loc','_',datestr(now,'ddmmyyyy_HHMM')]);
        
        saveas(gcf,savedfname,'png');
        
    % Unisensory auditory
    else
        
        figure();
        set(gcf,'Units','normalized','Position',[0.05 0.25  0.5 0.25],...
            'PaperPositionMode','auto');
        
        taskLevels = unique(behavData.accuracies.locAud);
        taskVarName = 'locAud';
        
        iFig = 1;
        
        for i = 1:length(taskLevels)
            acc = behavData.accuracies.calcaccuracy_resp(behavData.accuracies.(taskVarName) == taskLevels(i),:)';
            
            subplot(1,numel(taskLevels),iFig);
            plot(taskLevels,acc,'LineWidth',1.5,'Marker','o');
            ylim([0 1]);
            if i == 1
                title(sprintf('\n'));
            else
                set(gca,'YTickLabel',[]);
            end
            set(gca,'XTick',taskLevels);
            
            iFig = iFig + 1;
        end
        
        suplabel(behavData.figTitle,'t',[.08 .08 .84 .84]);
        
        savedfname = fullfile(DEC_2_setupdir(info.inputParameter.expStage,'data_behav_sub',info.subID),...
            [info.subID,'_',info.inputParameter.expMode,'_',num2str(info.iDay),'_',...
            num2str(info.iSessionInDay),'_','qanalyse_loc','_',datestr(now,'ddmmyyyy_HHMM')]);
        
        saveas(gcf,savedfname,'png');
    end
    
% Bisensory
else
    
    figure();
    set(gcf,'Units','normalized','Position',[0.05,0.2,0.6,0.6],...
        'PaperPositionMode','auto');
    
    if unique(behavData.data.task) == 1 % auditory localization
        taskLevels = unique(behavData.accuracies.locAud);
        nonTaskLevels = unique(behavData.accuracies.locVis);
        taskVarName = 'locAud';
        nonTaskVarName = 'locVis';
        taskLocStr = 'Auditory location';
        nonTaskLocStr = 'Visual location';
    else
        taskLevels = unique(behavData.accuracies.locVis);
        nonTaskLevels = unique(behavData.accuracies.locAud);
        taskVarName = 'locVis';
        nonTaskVarName = 'locAud';
        taskLocStr = 'Visual location';
        nonTaskLocStr = 'Auditory location';
    end
    
    visLocLevels = unique(behavData.accuracies.locVis);
    audLocLevels = unique(behavData.accuracies.locAud);
    
    iFig = 1;
    for i = 1:numel(nonTaskLevels)
        for j = 1:numel(taskLevels)
            
            acc = behavData.accuracies.calcaccuracy_resp(...
                behavData.accuracies.(nonTaskVarName) == nonTaskLevels(i) & ...
                behavData.accuracies.(taskVarName) == taskLevels(j),:)';
            
            subplot(numel(visLocLevels),numel(audLocLevels),iFig);
            plot(repmat(audLocLevels,1,size(acc,2)),acc,'LineWidth',1.5,...
                'Marker','o');
            
            ylim([0 1]);
            if i == 1
                set(gca,'XTick',[]);
            elseif i == numel(nonTaskLevels)
                xlabel(sprintf('\n%.2f',taskLevels(j)));
                set(gca,'XTick',taskLevels);
            else
                set(gca,'XTick',[]);
            end
            
            if j == 1
                ylabel(sprintf('%.2f\n',nonTaskLevels(i)));
            else
                set(gca,'YTickLabel',[]);
            end
            
            if iFig == 1
                legend({'R+','R-'},'Location','NorthEast');
            end
            
            iFig = iFig + 1;
            
        end
    end
    
    suplabel(behavData.figTitle,'t');
    suplabel(nonTaskLocStr,'y');
    suplabel(taskLocStr,'x');
    
    savedfname = fullfile(DEC_2_setupdir(info.inputParameter.expStage,'data_behav_sub',info.subID),...
        [info.subID,'_',info.inputParameter.expMode,'_',num2str(info.iDay),'_',...
        num2str(info.iSessionInDay),'_','qanalyse_loc','_',datestr(now,'ddmmyyyy_HHMM')]);
    
    saveas(gcf,savedfname,'png');
    
end

end

function plotbadtrials(badTrials,info)

exclLevelsAll = 0:8;
exclLevelsBad = [NaN,1:8];

freqOfLevelsAll = calcfreqoflevels(badTrials,exclLevelsAll);
freqOfLevelsBad = calcfreqoflevels(badTrials(badTrials ~= 0),exclLevelsBad);

figure();
set(gcf,'Units','normalized','Position',[0.7,0.3,0.225,0.3],...
    'PaperPositionMode','auto');
    
h = bar(1:2,[freqOfLevelsAll;freqOfLevelsBad],'stacked');
ylabel('Proportion of trials');
ylim([0,1]);
set(gca,'XTick',[],'Box','off');
axesPosition = get(gca,'Position');
axesPosition(3) = 0.5*axesPosition(3);
set(gca,'Position',axesPosition);
if freqOfLevelsAll(1) < 1
    yMaxRightYaxis = 1-freqOfLevelsAll(1);
else
    yMaxRightYaxis = 1;
end
hNewAxes = axes('Position',axesPosition,'Color','none','YLim',[0 yMaxRightYaxis],...
                'YAxisLocation','right','XTick',[],'Box','off');
hLeg = legend(h,{'good','no resp','early resp','wrong hand','EEG artf','no EYE data','blink','saccade','wrong fixation'},...
    'Location','SouthEastOutside');
% Offsetting legend position
legendPosition = get(hLeg,'Position');
legendPosition(1) = axesPosition(1)+(axesPosition(3)*1.2);
set(hLeg,'Position',legendPosition);

savedfname = fullfile(DEC_2_setupdir(info.inputParameter.expStage,'data_behav_sub',info.subID),...
    [info.subID,'_',info.inputParameter.expMode,'_',num2str(info.iDay),'_',...
    num2str(info.iSessionInDay),'_','qanalyse_qual','_',datestr(now,'ddmmyyyy_HHMM')]);

saveas(gcf,savedfname,'png');

end

function out = calcfreqoflevels(X,levels)

resp = X(:,1);
nLevels = numel(levels);
out = zeros(1,nLevels);

for i = 1:nLevels
    out(i) = sum(resp == levels(i))/size(resp,1);
end

end