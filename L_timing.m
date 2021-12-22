

fixation white
Screen(w,'Flip');
WaitSecs(2);

fixation 
% Perform some initial Flip to get us in sync with retrace:
% tvbl is the timestamp (system time in seconds) when the retrace
% started. We need it as a reference value for our display timing:
tvbl = Screen('Flip',w);

t_startTrial(iTrial) = (tvbl+ifi);
numifis =1; % 0 start immediately (next possible retrace )
t_deadline = tvbl + numifis * ifi - 0.5 * ifi;
while 1 
Screen('Flip',w,t_deadline)
end 