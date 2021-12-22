clear all;
clc;

% subject ID
SubID = input('input the subjects ID: ', 's');

louder = 50/2;
threshold = 54/2;
% chose the SAME!!! attenuation as for the final experiment
A.attenuation = 1/30 * (1/(10^(2/20)))^threshold * (10^(2/20))^louder;
% duration of each stimulus
A.durStim     = 0.3;
% select the difference in dB
A.dBdiff      = 2;
% attenuation to reduce noise
A.offset      = -threshold*2;

% load relevant information
load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' SubID '/SNDinfo.mat']);

% frequencies to be tested
freqs = sort([250 500 1000 2000 4000 Sd.tarFreqs]);
nFreqs = length(freqs);

% start the audiometry
dB = zeros(2,nFreqs);

for iFreq=1:nFreqs
    
    fprintf('\n \n Next test frequency: %i \n \n', freqs(iFreq))
    % chose frequency
    A.fc = freqs(iFreq);
    dB(1,iFreq) = freqs(iFreq);
    % do audiometry
    dB(2,iFreq) = AudiometryComparison(A);
    
    % wait until button was released
    KbWait([], 1);
    
end

% select target frequencies only (for later use)
dBTar = zeros(2,5);
for iTarFreq = 1:5
    dBTar(:,iTarFreq) = dB(:,(dB(1,:)==Sd.tarFreqs(1,iTarFreq)));
end

% calculate the number of attenuation (for later use)
dB2    = abs(dB(2,:)/A.dBdiff);
dBTar2 = abs(dBTar(2,:)/A.dBdiff);

eval(['save /Users/anette/Documents/MATLAB/Experiment2/Tones/' SubID '/AudiometryComparison' ])