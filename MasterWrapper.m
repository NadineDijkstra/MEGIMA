% Walkthrough of the steps that led to the results in the paper
% 'Differential temporal dynamics of visual imagery and perception' by N.
% Dijkstra, P. Mostert, F. de Lange, S.E. Bosch & M.A.J. van Gerven
%
% For clarity, only from after preprocessing and preparation of source
% info, for full details see paper

%% Some parameters and settings

addpath('/vol/ccnlab-scratch1/naddij/shareIMAMEG/Analyses')
root     = '/vol/ccnlab-scratch1/naddij/shareIMAMEG/Data';
analyses = '/vol/ccnlab-scratch1/naddij/shareIMAMEG/Analyses';
subjects = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'};
nSubjects = length(subjects);

for sub = 1:nSubjects
    
    % get the subject info
    subjectID              = subjects{sub}';  
    subjectDir             = fullfile(root,subjectID);
    dataDir                = fullfile(subjectDir,'CleanData');
    
    %% Multivariate analyses   
    addpath(fullfile(analyses,'MultivariateAnalyses'))
    
    cfg                    = [];
    cfg.root               = root;
    cfg.nFold              = 5; % number of folds for cross-validation
    cfg.channels           = 'MEG';
    
    cfg.decfg.gamma        = 0.05;
    cfg.decfg.samplemethod = 'downsample'; % to balance trial numbers
    cfg.decfg.nMeanS       = 9; % average over x samples, for smoothing
    
    % calculate the temporal generalization matrices
    cfg.outputDir          = fullfile(subjectDir,'DecodingTG');
    
    % for perception
    cfg.dataFiles{1}       = fullfile(dataDir,'dataFP'); % training data
    cfg.dataFiles{2}       = fullfile(dataDir,'dataFP'); % test data
    
    cfg.conIdx{1,1}        = [101:108,201:208]; % the markers for each class
    cfg.conIdx{1,2}        = [109:116,209:216];
    cfg.conIdx{2,1}        = cfg.conIdx{1,1};
    cfg.conIdx{2,2}        = cfg.conIdx{1,2};
    
    DecodingCross(cfg)    
    
    % for imagery
    cfg.dataFiles{1}       = fullfile(dataDir,'dataIM');
    cfg.dataFiles{2}       = fullfile(dataDir,'dataIM'); % test data
    
    cfg.conIdx{1,1}        = [121:128,221:228];
    cfg.conIdx{1,2}        = [129:136,229:236];
    cfg.conIdx{2,1}        = cfg.conIdx{1,1};
    cfg.conIdx{2,2}        = cfg.conIdx{1,2};
    
    DecodingCross(cfg)
    
    
    % calculate the cross-decoding matrices
    cfg.outputDir          = fullfile(subjectDir,'DecodingCROSS');
    
    % form perception to imagery    
    cfg.dataFiles{1}       = fullfile(dataDir,'dataFP'); % training data
    cfg.dataFiles{2}       = fullfile(dataDir,'dataIM'); % test data
    
    cfg.conIdx{1,1}        = [101:108,201:208];
    cfg.conIdx{1,2}        = [109:116,209:216];
    cfg.conIdx{2,1}        = [121:128,221:228];
    cfg.conIdx{2,2}        = [129:136,229:236];
    
    DecodingCross(cfg)
    
    % from imagery to perception 
    cfg.dataFiles{1}       = fullfile(dataDir,'dataIM'); % training data
    cfg.dataFiles{2}       = fullfile(dataDir,'dataFP'); % test data
    
    cfg.conIdx{1,1}        = [121:128,221:228];
    cfg.conIdx{1,2}        = [129:136,229:236];
    cfg.conIdx{2,1}        = [101:108,201:208];
    cfg.conIdx{2,2}        = [109:116,209:216];
    
    DecodingCross(cfg)    
        
    %% Source reconstruction
    addpath(fullfile(analyses,'SourceReconstruction'))    
    
    % calculate source activity difference ERF
    cfg             = [];
    cfg.movie       = 'no';
    cfg.dataDir     = dataDir;
    cfg.output      = fullfile(subjectDir,'SourceReconstruction');
    cfg.subjectID   = subjectID;
    cfg.time        = [];%[0 1];
    cfg.lpfilter    = 'yes';
    cfg.dataName    = 'dataIM';
    cfg.conIdx{1}   = [101:108,121:128,221:228]; cfg.conIdx{2} = [109:116,129:136,229:236];
    SourcerecConDiff(cfg)
    
    cfg.dataName    = 'dataFP';
    cfg.time        = [];
    SourcerecConDiff(cfg)
    
    % seperate grid points into atlas regions for later averaging
    cfg             = [];
    cfg.plot        = 'no';
    cfg.diff        = 'yes'; % look at the difference or average over conditions
    cfg.dataSet     = 'dataIM';
    cfg.subjectID   = subjectID;
    cfg.sourceRecDir = fullfile(subjectDir,'SourceReconstruction');
    cfg.root        = root;
    SourcerecAtlas(cfg)
    
    cfg.dataSet     = 'dataFP';    
    SourcerecAtlas(cfg)       
    
end

%% Group analyses decoding results
% perception temporal generalization
cfg.dataName = 'dataFP';
cfg.dirName  = 'DecodingTG';
cfg.root     = root;
cfg.nPerm    = 1000; % number of permutations
cfg.alpha    = 0.05;
cfg.subjects = subjects;
load(fullfile(root,subjects{1},cfg.dirName,cfg.dataName),'time') % steal a time variabel
cfg.time     = time;

DecodingGA(cfg); 

% imagery temporal generalization
cfg.dataName = 'dataIM';
load(fullfile(root,subjects{1},cfg.dirName,cfg.dataName),'time') % steal a time variabel
cfg.time     = time;

DecodingGA(cfg);

% cross-decoding from perception to imagery
cfg.dataName = 'dataFP to dataIM';
cfg.dirName  = 'DecodingCROSS';
load(fullfile(root,subjects{1},cfg.dirName,cfg.dataName),'time') % steal a time variabel
cfg.time     = time;

DecodingGA(cfg);

% cross-decoding from imagery to perception
cfg.dataName = 'dataIM to dataFP';
load(fullfile(root,subjects{1},cfg.dirName,cfg.dataName),'time') % steal a time variabel
cfg.time     = time;

DecodingGA(cfg);

%% Group averaging source reconstructions
cfg = [];
cfg.dataName = 'dataFP_atlasRegions.mat';
cfg.movie = 'yes';
cfg.root = root;

SourcerecNorm(cfg,subjects)


cfg.dataName = 'dataIM_atlasRegions.mat';
SourcerecNorm(cfg,subjects)

%% Make figures
addpath(fullfile(analyses,'MakeFigures'))    

% Figure 'Decoding performance of perception and imagery over time'

% diagonal decoding
cfg            = [];
cfg.dataName   = 'dataFP';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingTG');
cfg.time       = linspace(-0.2,1,360);
cfg.ylim       = [0.4 1];
cfg.fillcolor  = 'red'; % color of the shaded area
MPdecodingDiag(cfg)

cfg.dataName   = 'dataIM';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingTG');
cfg.time       = linspace(-0.2,4.5,1410);
cfg.ylim       = [0.46 0.65];
cfg.fillcolor  = 'red'; % color of the shaded area
MPdecodingDiag(cfg); xlim([-0.2 1]); 


% temporal generalization matrices
cfg            = [];
cfg.dataName   = 'dataFP to dataFP';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingTG');
cfg.time{1}    = linspace(-0.2,1,360);
cfg.clim       = [0.175 0.825];
cfg.steps      = 0.05;
MPdecodingMatrix(cfg)

cfg.dataName   = 'dataIM to dataIM';
cfg.time{1}    = linspace(-0.2,4.5,1410);
cfg.clim       = [0.37 0.63];
cfg.steps      = 0.02;
MPdecodingMatrix(cfg)
xlim([-0.2 4]); ylim([-0.2 1]) 


% Temporal specificity
cfg            = [];
cfg.dataName   = 'dataFP to dataFP';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingTG');
cfg.time{1}    = linspace(-0.2,1,360);
GeneralisationMeasure(cfg); 


cfg.dataName   = 'dataIM to dataIM';
cfg.time{1}    = linspace(-0.2,4.5,1410);
GeneralisationMeasure(cfg); xlim([0.5 1])


% Source activations
% 90 ms sourceplot perception
cfg            = [];
cfg.dataName   = 'dataFP_atlasRegions';
cfg.time       = 0.090;
cfg.dataDir    = fullfile(root,'GroupResults','SourceReconstruction');
cfg.clim       = [-7 7];
MPsourcePlot(cfg)

% 160 ms sourceplot
cfg.time       = 0.160;
MPsourcePlot(cfg)

% 700 ms sourceplot imagery
cfg.dataName   = 'dataIM_atlasRegions';
cfg.time       = 0.2; % here 0 is still offset instead of onset retro-cue
cfg.clim       = [-2 2];
MPsourcePlot(cfg)

% 1073 ms sourceplot
cfg.time       = 0.573;
MPsourcePlot(cfg)


% Figure 'Generalization between perception and imagery over time'

% from perception to imagery
cfg            = [];
cfg.dataName   = 'dataFP to dataIM';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingCROSS');
cfg.time{1}    = linspace(-0.2,1,360);
cfg.time{2}    = linspace(-0.2,4.5,1410);
cfg.clim       = [0.37 0.63];
cfg.steps      = 0.02;
MPdecodingMatrix(cfg); xlim([-0.2 4])

% from imagery to perception
cfg.dataName   = 'dataIM to dataFP';
cfg.time{2}    = linspace(-0.2,1,360);
cfg.time{1}    = linspace(-0.2,4.5,1410);
cfg.clim       = [0.37 0.63];
cfg.steps      = 0.02;
MPdecodingMatrix(cfg); ylim([-0.2 4])

% PI generalization at 90 ms
cfg            = [];
cfg.dataName   = 'dataFP to dataIM';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingCROSS');
cfg.timePoint  = 0.090; % time of training point to plot
cfg.time{1}    = linspace(-0.2,1,360);
cfg.time{2}    = linspace(-0.2,4.5,1410);
cfg.ylim       = [0.45 0.60];
cfg.fillcolor  = 'red'; % color of the shaded area
MPdecodingGeneralize(cfg); xlim([-0.2 4]);

% PI generalization at 160 ms
cfg.timePoint  = 0.160;
MPdecodingGeneralize(cfg); xlim([-0.2 4]);

% PI generalization at 210 ms
cfg.timePoint  = 0.210;
MPdecodingGeneralize(cfg); xlim([-0.2 4]);

% PI generalization at 500 ms
cfg.timePoint  = 0.500;
MPdecodingGeneralize(cfg); xlim([-0.2 4]);

% IP generalization at 400 ms
cfg            = [];
cfg.dataName   = 'dataIM to dataFP';
cfg.dataDir    = fullfile(root,'GroupResults','DecodingCROSSsmooth');
cfg.timePoint  = 0.400; % time of training point to plot
cfg.time{1}    = linspace(-0.2,4.5,1410);
cfg.time{2}    = linspace(-0.2,1,360);
cfg.ylim       = [0.45 0.65];
cfg.fillcolor  = 'red'; % color of the shaded area
MPdecodingGeneralize(cfg); 

% IP generalization at 700 ms
cfg.timePoint  = 0.700; % time of training point to plot
MPdecodingGeneralize(cfg); 

% IP generalization at 1073 ms 
cfg.timePoint  = 1.073; % time of training point to plot
MPdecodingGeneralize(cfg); 


