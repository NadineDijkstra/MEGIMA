function megima_init(subject, session, location)

%   megima_init(subject, session, location)
%   
%   meg imagery task initialization code. Calls the function megima_task
%   with the desired set of P.
%
%   subject:       String with subject ID  
%   session:       Session integer
%   location:      'debug'        = 0 = Laptop
%                  'behavioural'  = 1 = Behavioural lab
%                  'meg'          = 2 = MEG scanner
%                  'mri'          = 3 = MRI scanner
%
%   If no other input arguments are given, the code will assume it needs 
%   to run in demonstration mode. If only situation isn't entered, the default
%   assumption is that you're running the code for the scanner.
%
%   Written by SEB in FEB 2017, adapted by Nadine Dijkstra APR 2017

addpath('functions');
P = struct;

%% Initialize randomness & keycodes
P.rngSeed = rng('shuffle', 'twister');
KbName('UnifyKeyNames');
RestrictKeysForKbCheck([]); % reenable all keys for KbCheck
Screen('Preference', 'SkipSyncTests', 1); 

%% User-defined input
if nargin == 0
    subject = 'S00';
    session = 1;
    location = 'debug';
end

switch location
    case 'debug'
        P.situation = 0; % 0 = Desktop, 1 = Behavioural, 2 = Trio, 3 = Random computer
        P.windows   = true;
        P.dataPath  = 'debugData/';
        P.leftKey   = 'LeftArrow';
        P.rightKey  = 'RightArrow';
        P.keys      = {'a','s','d','f'};
    case 'behavioural'
        P.situation = 1; % 0 = Desktop, 1 = Behavioural, 2 = Trio, 3 = Random computer
        P.windows   = false;
        P.dataPath  = 'behaviouralData/';
        P.leftKey   = 97;
        P.rightKey  = 98;
        P.keys      = [97,98,99,100]; % CHECK THIS!!! <--------
        P.leftKeyOff = 65;
        P.rightKeyOff = 66;
    case 'meg'
        P.situation = 2;
        P.windows   = false;
        P.dataPath  = 'megData/';
        P.leftKey   = 97;
        P.rightKey  = 98;
        P.keys      = [97,98,99,100];
        P.leftKeyOff = 65;
        P.rightKeyOff = 66;
    case 'mri'
        P.situation = 3;
        P.windows   = false;
        P.dataPath  = 'mriData/';
        P.leftKey   = 97;
        P.rightKey  = 98;
        P.leftKeyOff = 65;
        P.rightKeyOff = 66;
end

if ~exist(P.dataPath, 'dir'), mkdir(P.dataPath); end

P.subject = subject;
P.theDate = datestr(now, 'yyyymmdd');

% Determine current session
[P.session, P.sessionName] = current_session([P.theDate '_megima_' subject], session, P.dataPath); 
if strcmp(P.session, 'abort'), error('No response given, exiting.'); end

%% Screen parameters
[frameHz,pixPerDeg,calibrationFile] = get_monitor_info(location);

P.frameHz          = frameHz;
P.pixPerDeg        = pixPerDeg;
P.calibrationFile  = calibrationFile;
P.screen           = 0; % main screen
P.resolution       = Screen('Rect', P.screen); %[750 50 1250 550]; %Screen('Rect', P.screen); %
P.backgroundColour = 127;
P.fontName         = 'Arial';
P.fontSize         = 20;

%% Gamma-corrected CLUT
P.meanLum = 0.5;
P.contrast = 1;
P.stimContrast = 0.2;
P.amp = P.meanLum * P.contrast;

nColours = 255;	% number of gray levels to use in mpcmaplist, should be uneven
mpcMapList = zeros(256,3);	% color look-up table of 256 RGB values, RANGE 0-1
tempTrial = linspace(P.meanLum-P.amp, P.meanLum+P.amp, nColours)';	% make grayscale gradient

P.meanColourIdx = ceil((nColours)/2) -1; % mean colour number (0-255) of stimulus (used to index rows 1:256 in mpcmaplist)
P.ampColourIdx  = floor(P.stimContrast*(nColours-1)); % amplitude of colour variation for stimulus

if P.situation ~= 0 && ~isempty(P.calibrationFile)
    load(P.calibrationFile);	% function loads inverse gamma table and screen dacsize from most recent calibration file
    
    mpcMapList(1:nColours,:) = repmat(tempTrial, [1 3]);
    mpcMapList(256,1:3) = 1;
    mpcMapList = round(map2map(mpcMapList,gamInverse));
    
    P.CLUT = mpcMapList;
end

%% Experiment parameters
P.nStim = 2; % per trial
P.nStimuli = 16; % actual number of stimuli
if strcmp(location,'debug')
    P.nTrialsPerRun = 5;
else
    P.nTrialsPerRun = 24;
end

P.runIdx = ((P.session-1)*P.nTrialsPerRun+1):(P.session*P.nTrialsPerRun);

load(sprintf('trialMatrix_%s.mat',subject))
P.trialMatrix = trialMatrix(P.runIdx,:);

%% Stimulus parameters
P.stimSizeInDegree = 5;
P.distFromFixInDeg = 1;
P.distFromFixInPix = round(P.pixPerDeg * P.distFromFixInDeg);
P.stimSizeInPix = round([P.stimSizeInDegree*P.pixPerDeg P.stimSizeInDegree*P.pixPerDeg]);	% Width, height of stimulus
P.baseRect = [0 0 160 160]; % stimulus frame

%% Time Parameters                                                        
P.fix               = 2;
P.preCueDur         = 0.3;
P.cueDur            = 0.4;
P.postCueDur        = 0.3;
P.stimDur           = 0.8;
P.ISIrange          = [0.4, 0.6];
P.maskDur           = 0.5;
P.postStimRange     = [0.4, 0.6];
P.retroCueDur       = 0.5;
P.imageryDur        = 3.5; 
P.responseDur       = 3; % so slow, then continue
P.ITIrange          = [1.8, 2.2]; 
P.confDur           = 3; % time the vividness rating will be on screen for (=response epoch)
P.confDim           = P.confDur - 1; % how far into the conf task should the grating start to dim?
P.postConfDur       = 0.2;

%% Fixation parameters
P.fixRadiusInDegree = .6;

%% Retrocue parameters
P.retroCueID = {'1', '2'};
P.retroCueColour = [255 255 255];

%% Vividness rating parameters
P.lineLengthPix            = 300;
P.sliderLengthPix          = 40;
P.matchMean                = P.backgroundColour;
P.matchContrast            = 0.4;
P.matchAmp                 = P.matchMean * P.matchContrast;
P.lineDiamPix              = 5;
P.midX                     = (P.resolution(3)-P.resolution(1))/2+P.resolution(1);
P.midY                     = (P.resolution(4)-P.resolution(2))/2+P.resolution(2);
P.lineRect                 = [P.midX-P.lineLengthPix/2 P.midY-P.lineLengthPix/2 P.midX+P.lineLengthPix/2 P.midY+P.lineLengthPix/2];
P.sliderRect               = [P.midX-P.sliderLengthPix/2 P.midY-P.sliderLengthPix/2 P.midX+P.sliderLengthPix/2 P.midY+P.sliderLengthPix/2];
P.yOffset                  = 120;

%% Trigger parameters
P.triggerStart = 66; % start trial
P.triggerRetrocue = 70; % + identity of the cue 
P.triggerProbe = 80; 
P.triggerResponse = 90; % + identity of the response
% perception = first-or-second stim   | stim identity e.g. 110 for 1st and stim 10
% imagery =    first-or-second stim+2 | stim identity e.g. 301 or 402

%% Run sthe expefriment
megima_task(P);

end 