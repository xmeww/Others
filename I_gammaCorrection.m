function gammaCorrection(S)



% -------------------------------------------------------------------------
% INPUT STRUCT
%
% if no input is given to the system, S will be created and all default
% settings will be used
if ~exist('S', 'var'),           S = [];               end              

if ~isfield(S, 'debug'),         S.debug  = 0;          end

if ~isfield(S, 'nSteps'),        S.nSteps = 18;          end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              START                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    
    %%                 START PASYCHTOOLBOX
    % =====================================================================
    % Make sure we're running on PTB-3
    AssertOpenGL;
    
    % always needs to be put if using KbName!!!!!
    KbName('UnifyKeyNames');   
    % define names of keys
    quit = KbName('q');
    
    
    % ---------------------------------------------------------------------
    % DEBUG
    % check for breakpoints
    brkp = dbstatus;
    if (S.debug == 1) || ~isempty(brkp)
        % simplyfies debugging by making windows transparent
        PsychDebugWindowConfiguration;
    end
 
    
    

    
    %%                 OPEN & PREPARE SCREENS
    % =====================================================================
    
    % Get the list of Screens and choose the one with the highest screen
    % number. Screen 0 is, by definition, the display with the menu bar.
    % Often when two monitors are connected the one without the menu bar is
    % used as the stimulus display.  Chosing the display with the highest
    % dislay number is a best guess about where you want the stimulus
    % displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    
    % Open double-buffered window: Optionally enable stereo output if
    % stereo > 0. 1 = Frame-sequential, 4 = Stereoscope (dual-display), 8 =
    % Anaglyph Red-Blue.
    PsychImaging('PrepareConfiguration');
    
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'CheckOnly');
    
    
    backgroundvalue = 0;
    % Prepare setup of imaging pipeline for onscreen window. This is the
    % first step in the sequence of configuration steps
    [w, winRect]=PsychImaging('OpenWindow', screenNumber, 255*backgroundvalue);
    
    % Using this blendfunction the information which is presented to the
    % screen is non-transperent, i.e. info coming later will overwrite info
    % from before. If you want an additive screen function you can do the
    % following:
    % Screen(w,'BlendFunction',GL_ONE, GL_ONE);
    % BlendFunctions can rapidly be changed all over the script.
    Screen(w,'BlendFunction',GL_ONE, GL_ZERO);
    % PsychColorCorrection('SetColorClampingRange', w, 0.01, 0.99)
    
    % define window height
    win_w = (winRect(3)-winRect(1));
    win_h = (winRect(4)-winRect(2));
    
    
    

    
 %%                 TEXTURE
    % =====================================================================
    
    grayImg = ones(win_h, win_w);
    grayTex = Screen('MakeTexture', w, grayImg, [], [], 2);

    
    
    
    
    %%                MEASURING THE REFRESH RATE
    % =====================================================================
    
    % Finding out the RefreshRate of the monitor Query nominal framerate as
    % returned by Operating system: If OS returns 0, then we assume that we
    % run on a flat-panel with fixed 60 Hz refresh interval.
    framerate=Screen('NominalFramerate', w);
    if (framerate==0)
        framerate=60;
    end;
    
    % gives ms of one time frame eg. 1/60 = 16.6667 ms
    ifinominal=1 / framerate;
    fprintf('The refresh interval reported by the operating system is %2.5f ms.\n', ifinominal*1000);
    fprintf('The measured framerate is %2.5f ms.\n', framerate);
    
    % Perform a calibration loop to determine the "real" interframe
    % interval for the given gfx-card + monitor combination:
    %
    %
    % Measure monitor refresh interval again, just for fun... This will
    % trigger a calibration loop of minimum 100 valid samples and return
    % the estimated ifi in 'ifi': We require an accuracy of 0.05 ms ==
    % 0.00005 secs. If this level of accuracy can't be reached, we time out
    % after 20 seconds...
    [ifi nvalid stddev ]= Screen('GetFlipInterval', w, 100, 0.00005, 5);
    fprintf('Measured refresh interval, as reported by "GetFlipInterval" is %2.5f ms. (nsamples = %i, stddev = %2.5f ms)\n', ifi*1000, nvalid, stddev*1000);
    
    
    
    
    
    
    
    %%                INITIALIZE VARIABLES
    % =====================================================================    
    %----------------------------------------------------------------------
    % INITIALIZING
    stepValue  = linspace(0, 255, S.nSteps);
    finalSteps = (length(stepValue));
    


    
    
    %%              START PRESENTATION
    for iFlip = 1:finalSteps
        
               
        % -----------------------------------------------------------------
        % VISUAL STIMULUS        
        % Draw noise
        text = ['actualValue presented is: ' num2str(stepValue(iFlip))];        
        Screen('DrawTexture', w, grayTex, [], [], [], [], [], stepValue(iFlip));
        Screen('TextFont', w, 'Arial');
        Screen('TextSize', w, 18);
        Screen('DrawText', w, text, win_w/2,  win_h/2, [255 0 0]);
        Screen('DrawingFinished', w, 0);
        
        
        % -----------------------------------------------------------------
        % FLIP
        %
        % Flip: The clearmode argument specifies if flip should clear the
        % drawing buffer after flip (=0 - default), keep it "as is" for
        % incremental drawing/updating of stims (=1) or don't do anything
        % to the framebuffer at all (=2). We return the timestamp, when VBL
        % starts in tvbl: This is when the front- and back drawing surfaces
        % get exchanged and it is the crucial reference value for computing
        % the 'tDeadline' presentation deadline for the next 'Flip'
        % command. The rasterbeam-position (scanline) when the measurement
        % was taken is returned in beampos(iFlip), the time when flip
        % returned to Matlab is returned in flipfin(iFlip), estimated
        % stimulus onset time aka end of VBL is returned in so(iFlip).
        %
        % The first value "tvbl" is needed for tDeadline calculation if one
        % wants to emulate WaitBlanking of old PTB - see formula above.
        % beampos > screen height means that flip returned during the VBL
        % interval. Small values << screen height are also ok, they just
        % indicate either a slower machine or some types of flat-panels...
        Screen('Flip', w);
        WaitSecs(1);
        kbWait;
        
        % -----------------------------------------------------------------
        % QUIT
        [~,~,keyCode] = KbCheck;
        % stopping the loop if user presses 'q'
        if keyCode(quit)
            break;
        end
        
        
        
    end % Draw next frame
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          END LOOP                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%                         CLOSE
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % SCREENS
    %
    % Close display: If we skipped/missed any presentation deadline during
    % Flip, Psychtoolbox will automatically display some warning message on
    % the Matlab console:
    Screen('CloseAll');
    
    % ---------------------------------------------------------------------
    % Shutdown realtime scheduling:
    Priority(0);
    
    

    
    
catch  %#ok<*CTCH>
    % This "catch" section executes in case of an error in the "try"
    % section above. Importantly, it closes the onscreen window if its open
    % and shuts down realtime-scheduling of Matlab
    % ---------------------------------------------------------------------
    % SCREENS
    %
    % Close display: If we skipped/missed any presentation deadline during
    % Flip, Psychtoolbox will automatically display some warning message on
    % the Matlab console:
    Screen('CloseAll');    
    
    % ---------------------------------------------------------------------
    % Shutdown realtime scheduling:
    Priority(0);  
    
    
    % ---------------------------------------------------------------------
    % output error that occured
    psychrethrow(psychlasterror); 
end

% Done!!!! :)