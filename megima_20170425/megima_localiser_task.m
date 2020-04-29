function Localizer(Parameters, emulate)

%   Localizer(Parameters, emulate)
% 
%   Letter task with low-contrast gratings for classifier training.
%
%   Arguments in:
%
%       Parameters  =   Struct containing various parameters
%       emulate 	=   0 = Triggered by scanner
%                       1 = Trigger by keypress
%
%   Written by SEB in June 2012

%% Behavioural data
Behaviour                   = struct;
Behaviour.eventTime         = [];
Behaviour.response          = [];
Behaviour.responseTime      = [];

%% Initialize randomness & keycodes
setupRand;
setupKeyCodes;

%% Stimulus conditions 
if emulate % In manual start there are no dummies
    Parameters.dummies = 0;
    Parameters.overrun = 0;
end

%% Initialise experiment
try 
HideCursor; 

if ~emulate
    b1 = Bitsi_Scanner('/dev/ttyS1');
    b2 = Bitsi_Scanner('/dev/ttyS2');
end

[wPtr, rect] = Screen('OpenWindow', Parameters.screen, Parameters.background, Parameters.resolution); % Opens fullsize window with window pointer 'wPtr' and background colour

Screen('TextStyle', wPtr, 1)
Screen('Preference', 'DefaultFontSize', Parameters.fontSize);
Screen('Preference', 'DefaultFontName', Parameters.fontName);
Screen('TextSize', wPtr, Parameters.fontSize);
Screen('TextFont', wPtr, Parameters.fontName);

hardwareCLUT = Screen('LoadCLUT', Parameters.screen);
Screen('BlendFunction', wPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% if ~emulate
    Screen('LoadCLUT', wPtr, Parameters.CLUT); % Loads gamma-corrected CLUT
% end

%% Loading screen
txt = 'Creating stimuli';			
txtLoc = [50 50];
Screen('FillRect', wPtr, Parameters.background, Parameters.resolution);
Screen('DrawText', wPtr, txt, txtLoc(1), txtLoc(2), 255); % Adds text to fixation point

startCreatingStimuli = Screen('Flip', wPtr); % Draws fixation point

%% Create fixation bullseye
fixTexture = Screen('MakeTexture', wPtr, Parameters.fixBackground);
Screen('FillArc',fixTexture,1,CenterRect([0 0 Parameters.fixCenterDiameter Parameters.fixCenterDiameter],Parameters.fixRect),0,360);
Screen('FrameArc',fixTexture,1,CenterRect([0 0 Parameters.fixCenterDiameter*3 Parameters.fixCenterDiameter*3],Parameters.fixRect),0,360,Parameters.fixCenterDiameter/2,Parameters.fixCenterDiameter/2);

%% Create event cue
cueTexture = Screen('MakeTexture', wPtr, Parameters.fixBackground);
Screen('FillArc',cueTexture,1,CenterRect([0 0 Parameters.fixCenterDiameter Parameters.fixCenterDiameter],Parameters.fixRect),0,360);

%% Create stimuli
Imglist = zeros(2);

img1 = Parameters.stim1.*Parameters.mask;
img1 = round(img1*Parameters.ampColourIdx/2+Parameters.meanColourIdx);
Imglist(1) = Screen('MakeTexture', wPtr, img1);

img2 = Parameters.stim2.*Parameters.mask;
img2 = round(img2*Parameters.ampColourIdx/2+Parameters.meanColourIdx);
Imglist(2) = Screen('MakeTexture', wPtr, img2);

%% Standby screen
txt = 'Waiting for trigger to begin';
txtLoc = [50 50];
Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
Screen('DrawText', wPtr, txt, txtLoc(1), txtLoc(2), Parameters.foreground2); % Adds text to fixation point

creatingStimuliDone = Screen('Flip', wPtr); % Draws text message, reader for experiment

%% Wait for start of experiment
if emulate
    KbWait([],2);
    
    Screen('FillRect', wPtr, Parameters.background, rect);
    Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect));
    Screen('Flip', wPtr);

    WaitSecs(Parameters.TR*Parameters.dummies);    

else
    b2.clearResponses();
    firstScan = 0;
    while firstScan == 0
        while b2.numberOfResponses() == 0
            WaitSecs(0.001);
        end
        [resp] = b2.getResponse(0.001, true);
        if resp == 97
            firstScan = 1;
        end
    end
    
    if Parameters.dummies > 0
        Screen('FillRect', wPtr, Parameters.background, rect);
        Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect));
        Screen('Flip', wPtr);
        WaitSecs(Parameters.TR*Parameters.dummies);
    end
end

%% Begin main experiment 
refreshDur = Screen('GetFlipInterval',wPtr);
slack = refreshDur / 2;

%% Fixation block before start of run
if Parameters.fix > 0;
    Screen('FillRect', wPtr, Parameters.background, rect);
    Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect));
    Screen('Flip', wPtr);
    WaitSecs(Parameters.TR*Parameters.fix);
end

Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
currTime = Screen('Flip', wPtr); % Draws fixation point
T.startTime = currTime;

for iBlock = 1:Parameters.nPerRun 
    
    Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
    currTime = Screen('Flip', wPtr); % Fixation presentation
    T.blockStart(iBlock) = currTime - T.startTime;
 
    if mod(iBlock,2)==1 % If iBlock is odd
        
        while (GetSecs - T.startTime < Parameters.blockSwitch(iBlock) - slack) % Subtracting half flip interval gives plenty of slack to catch the desired frame for change
            WaitSecs(0.001);
        end
    
        Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
        currTime = Screen('Flip', wPtr); % Fixation presentation
        T.blockEnd(iBlock) = currTime - T.startTime;
        
    elseif mod(iBlock,2)==0 % If iBlock is even
        
        stimulus = Randi(2); % First stimulus is randomly chosen
        prevKeypress = 0;  % If previously key was pressed
        stimNum = 0;
        isEvent = 0;
        wasEvent = 0;
        if ~emulate
            b1.clearResponses();
        end
        
        while (GetSecs - T.startTime < Parameters.blockSwitch(iBlock) - slack) % Subtracting half flip interval gives plenty of slack to catch the desired frame for change            

            if rand < Parameters.probOfEvent && wasEvent~=1;
                isEvent = 1;
            end
            
            if isEvent==1 && wasEvent==0
                
                stimNum = stimNum + 1;
                
                % Draw event
                Screen('DrawTexture', wPtr, Imglist(stimulus), [], [Parameters.stimRect]');
                Screen('DrawTexture', wPtr, cueTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
                currTime = Screen('Flip', wPtr);
                T.check(iBlock,stimNum) = currTime - T.startTime;
                T.event(iBlock,stimNum) = 1;
                Behaviour.eventTime(iBlock,stimNum) = currTime - T.startTime;
                
                stimulus = 3-stimulus;
                isEvent = 0;
                wasEvent = 1;                
            
            elseif isEvent==0 && wasEvent==1
                
                stimNum = stimNum + 1;
                
                % Draw event
                Screen('DrawTexture', wPtr, Imglist(stimulus), [], [Parameters.stimRect]');
                Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
                currTime = Screen('Flip', wPtr);
                T.check(iBlock,stimNum) = currTime - T.startTime;
                T.event(iBlock,stimNum) = 1;
                
                stimulus = 3-stimulus;
                isEvent = 0;
                wasEvent = 0;                

            elseif isEvent==0 && wasEvent==0
                
                stimNum = stimNum + 1;
                
                % Prepare next video page...
                Screen('DrawTexture', wPtr, Imglist(stimulus), [], [Parameters.stimRect]');
                Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
                currTime = Screen('Flip', wPtr);
                T.check(iBlock,stimNum) = currTime - T.startTime;
                T.event(iBlock,stimNum) = 0;
                
                stimulus = 3-stimulus;
                
            end
            
            while (GetSecs - currTime < Parameters.stimDur - slack) % Subtracting half flip interval...
                %% Behavioural response
                if ~emulate
                    
                    timeout = 0.001;
                    [resp time] = b1.getResponse(timeout, 'true');
                    if resp == 97
                        Behaviour.response = [Behaviour.response; resp];    
                        Behaviour.responseTime = [Behaviour.responseTime; time - T.startTime];
                    end
                    WaitSecs(0.001);
                    
                elseif emulate
                    
                    [keyPress keyTime key] = KbCheck;
                    if key(KeyCodes.Escape)
                        % Abort screen
                        Screen('FillRect', wPtr, Parameters.background, rect);
                        DrawFormattedText(wPtr, 'Experiment was aborted!', 'center', 'center', Parameters.foreground2);
                        Screen('Flip',wPtr);
                        WaitSecs(0.5);
                        ShowCursor;
                        Screen('LoadCLUT', wPtr, hardwareCLUT);
                        Screen('CloseAll');
                        disp(' ');
                        disp('Experiment aborted by user!');
                        disp(' ');
                        return
                    end
                    if keyPress
                        if ~prevKeypress
                            prevKeypress = 1;
                            Behaviour.response = [Behaviour.response; find(Key)];    % True & false will be converted to 1 & 0, empty will mean nothing is added
                            Behaviour.responseTime = [Behaviour.responseTime; keyTime - T.startTime];
                        end
                    else
                        if prevKeypress
                            prevKeypress = 0;
                        end
                    end
                end

            end
          
            Screen('DrawTexture', wPtr, fixTexture, Parameters.fixRect, CenterRect(Parameters.fixRect, rect)); 
            currTime = Screen('Flip', wPtr); % Fixation presentation
            T.blockEnd(iBlock) = currTime - T.startTime;
            
        end % end while even block
    end % end if odd or even block
end % end iBlock: block per run
    
%% Last screen
currTime = Screen('Flip', wPtr); % One last flip to properly record time of end of last replay_trial, end of last runnum, and end of run.
T.endTime = currTime;

%% Fixation block after of run
if Parameters.fix > 0
    Screen('DrawTexture', wPtr, fixTexture);
    Screen('Flip', wPtr);
    WaitSecs(Parameters.TR*Parameters.fix);
end

WaitSecs(Parameters.TR * Parameters.overrun);

%% Clean up
if ~emulate
    close(b1);
    close(b2);
    delete(instrfind);
end

ShowCursor;
Screen('LoadCLUT', wPtr, hardwareCLUT); 
Screen('CloseAll');

%% Experiment duration
newLine;
disp('Experiment done');
T.experimentDuration = T.endTime - T.startTime;
durationMinutes = floor(T.experimentDuration/60);
durationSeconds = mod(T.experimentDuration, 60);
disp(['Cycling lasted ' num2str(durationMinutes) ' minutes, ' num2str(durationSeconds) ' seconds']);
newLine;

catch 

%     Screen('LoadCLUT', 0, hardwareCLUT);
    if ~emulate
        close(b1);
        close(b2);
        delete(instrfind);
    end

    Screen('CloseAll');
    ShowCursor;
    psychrethrow(psychlasterror);
end % try - catch

%% Save workspace
save(['data' filesep Parameters.sessionName]);

fprintf('It took %0.6f secs to load and create the images',creatingStimuliDone-startCreatingStimuli);
fprintf('\n');
