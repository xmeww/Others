% generates all frequencies that are equally spaced on a log scale

% spacing
logSpace = 1.13795;

iFreq = 1;
f(iFreq) = 200; % minimum frequency
mxF = max(f);
while mxF < 5000 % maximal frequency
    f(iFreq+1)=round(f(iFreq)*logSpace); % generates frequencies
    iFreq = iFreq + 1;
    mxF = max(f);
end

% target frequencies
ftar = f(15:2:23);