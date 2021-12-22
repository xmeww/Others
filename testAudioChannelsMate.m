function testAudioChannelsMate(deviceindex,nChannels,channelComp)
%testAudioChannels(nChannels [, deviceID])
%Plays white noise over the channel specified by presssing a number on the
%numeric pad. 
%Specify total number of audio channels. Can optionally also
%assign a deviceindex - use PsychPortAudio('GetDevices') to list audio
%devices. Uses default if none is assigned.


% Get screen id
screen = max(Screen('Screens')); %Output to external monitor by default
black = BlackIndex(screen); %Get black index for current screen. Probs ~0

[win, screenrect] = Screen('OpenWindow', screen, black); %Open window
Screen('Preference', 'SkipSyncTests', 0)

sr = 44100;

% Initialize PsychPortAudio sound driver
InitializePsychSound;

% Open audio device
paMaster = PsychPortAudio('Open',deviceindex,1+8,[],sr,nChannels);
% Start master immediately:
PsychPortAudio('Start',paMaster,0,0,1);

% Output volume
volume = 0.025;

paSlave = NaN(1,nChannels);
for i = 1:nChannels
    paSlave(i) = PsychPortAudio('OpenSlave',paMaster,1,1,i);
    PsychPortAudio('Volume',paSlave(i),1.0,volume*channelComp(i));
end

nSamples = 1*sr;

audio.loc = zeros(nSamples,nChannels,nChannels);

whitenoise = randn(1,nSamples); %Gaussian distribution white noise
whitenoise(whitenoise>1) = 1; %Cuts off anything outside 1/-1
whitenoise(whitenoise<-1) = -1;



for i=1:nChannels %For each channel:
    audio.loc(:,i,i) = whitenoise; %Fill the corresponding column of that channels matrix slice with white noise
%     audio.loc(:,:,i) = audio.loc(:,:,i) / max(max(abs(audio.loc(:,:,i)))); % normalize to be within [-1 1]
end

keyNames = cell(1,nChannels);
listKeyCodes = zeros(1,nChannels);

if nChannels < 10
    for i = 1:nChannels
        keyNames{i} = num2str(i);
        listKeyCodes(i) = KbName(keyNames{i});
    end
else
    for i = 1:9
        keyNames{i} = num2str(i);
        listKeyCodes(i) = KbName(keyNames{i});
    end
end
while 1
    [pressed, secs, keycode, deltasecs] = KbCheck; %Check for keypresses and collect data about them
    keyid = find(keycode); %Get a code for the key pressed (from the array of 1s and zeroes)
    if pressed && length(keyid) == 1 %If only one key was pressed then:
        if keycode(KbName('q'))
            break
        elseif ismember(keyid, listKeyCodes) %Else if it was one of our defined response keys:
            % Buffer the appropriate sound
            PsychPortAudio('FillBuffer', paSlave(listKeyCodes == keyid), whitenoise); %Buffer the sound in the audio device buffer
            % Present audio data
            PsychPortAudio('Start',paSlave(listKeyCodes == keyid)); %Start the audio device (equivalent of screen 'flip'), for one repetition, at the predicted onset time
        elseif nChannels > 9 && keycode(KbName('p'))
            % Buffer the appropriate sound
            PsychPortAudio('FillBuffer', paSlave(listKeyCodes == keyid), audio.loc(:,:,11)'); %Buffer the sound in the audio device buffer
            % Present audio data
            PsychPortAudio('Start', paSlave(listKeyCodes == keyid)); %Start the audio device (equivalent of screen 'flip'), for one repetition, at the predicted onset time
        end
    end
end
sca;

%Shut down the audio device
PsychPortAudio('Close'); 