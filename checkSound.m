clear all;
clc;

S.SubID = 'Sub13';


% -------------------------------------------------------------------------
% determine loudness per frequency
load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/Audiometry.mat'])
dBdiff    = 2;
% loudness of masking frequencies
audThresM  = abs(dB(2,:))/dBdiff;
for iTimes = 1:10
    S.audScaling(iTimes) = (1/(10^(dBdiff/20)))^audThresM(iTimes);
end

% loudness of target Frequencies
audThresT = abs(dBTar(2,:))/dBdiff;
for iTimes = 1:5
    S.audScalingTar(iTimes) = (1/(10^(dBdiff/20)))^audThresT(iTimes);
end

% loudness masking frequencies should be presented at
loudM    = 56; % loudness in dB
audLoudM = loudM/dBdiff;
S.loudDB = (10^(dBdiff/20))^audLoudM;

% loudness target frequencies should be presented at
loudT       = 50;
audLoudT    = loudT/dBdiff;
S.loudDBTar = (10^(dBdiff/20))^audLoudT;


% -------------------------------------------------------------------------
% AUDITORY
% percent aud signal strength
S.SNRaud = 1/30;      

% percent aud noise strength
S.SNRaudNoise = 1/30; 


% -------------------------------------------------------------------------
% load mask sound if necessary
load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/SNDinfo.mat']);

% SOUND 
%----------------------------------------------------------------------
nMaskFreq = length(Sd.PossFreq);
for iT = 1:5
    fr(iT) = find(Sd.PossFreq == Sd.tarFreqs(iT));
    Sd.beep(iT,:) = Sd.beep(iT,:) * S.audScalingTar(iT) * S.loudDBTar * S.SNRaud;
end

load(['/Users/anette/Documents/MATLAB/Experiment2/Tones/' S.SubID '/S1/maskSND_T1.mat']);


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


for iFreq = 1:5
    tarFreqIdx = iFreq;
    tarFreq    = fr(tarFreqIdx);
    maskOnNext = ones(1,nMaskFreq);
    maskOnNext(tarFreq-3:tarFreq+3) = 0;
    
    mxSnd(iFreq) = sum(max(maskSND(find(maskOnNext),:),[],2)) + max(Sd.beep(iFreq,:),[],2);
    
end

