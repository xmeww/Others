function rec = getaudio(D)
%  [subj, audio] = getaudio(subj, mypath, audio, audiotype)
%  getaudio retrieves audio and subject data either from the CIPIC database or recordings
%  Input:
%       D: struct with following fields: 
%           subjID: subject ID. 
%           stimLications: vector of stimulus locations. 
%           recIDs: vector of recording serial numbers to be included. 
%           cutoff: the cutoff frequency of the low-pass filtering. 
%           folderPath: full path of the folder which contains the recorded
%               sounds. 
%           stimDuration
%  Output:
%       rec: cell array (nummber of stim locations x number of recordings)
%           containing the recorded sounds. 

rec = cell(length(D.stimLocations),length(D.recIDs));

% Read audio recordings
for i = 1:length(D.stimLocations)
    for k = 1:length(D.recIDs)
        fname = fullfile(D.folderPath,sprintf('%s_%.1f_%d_low_%s.wav',D.subID,D.stimLocations(i),D.recIDs(k),D.cutoff));

        [y,Fs] = audioread(fname);

        if Fs ~= D.audioFreq
            error('Audio sampling frequency is used inconsistently.')
        end
        
        rec{i,k} = y;
        
        if abs(length(y)-round(D.stimDuration * D.audioFreq)) > 10
            error('Audio duration is used inconsistently.')
        end
    end
end
