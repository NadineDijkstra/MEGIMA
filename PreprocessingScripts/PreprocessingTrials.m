function PreprocessingTrials(cfg0,subjectdata)
% function PreprocessingTrials(subjectID)
% subjectdata is a structure containing all the relevant info about the
% subject such as outputDir, datadir etc. 

% set defaults
ft_defaults

% cd to root
root = '/vol/ccnlab-scratch1/naddij/ImaMEG';
cd(root)

%% Some settings
saveDir                     = fullfile(root,subjectdata.outputDir,'PreprocData');
if ~exist(saveDir,'dir'); mkdir(saveDir); end
cfgS                        = [];
cfgS.dataset                = fullfile(subjectdata.subjectdir,subjectdata.datadir);
cfgS.continuous             = 'yes';
cfgS.dftfilter              = 'yes';
cfgS.demean                 = 'no'; % baseline correction comes later
cfgS.padding                = 10;

% for the different segments, define information
eventvalue                  = cfg0.eventvalue;
prestim                     = cfg0.prestim;
poststim                    = cfg0.poststim;
saveName                    = cfg0.saveName;

%% Loop over segments
for s = 1:cfg0.nSeg
    
    fprintf('\t Getting the data for segment %d \n',s)
    % define it
    cfg                         = cfgS;
    cfg.trialfun                = 'ft_trialfun_general';
    cfg.trialdef.eventtype      = 'UPPT001';
    cfg.trialdef.eventvalue     = eventvalue{s}; % stimulus 1
    cfg.trialdef.prestim        = prestim{s};
    cfg.trialdef.poststim       = poststim{s};
    cfg                         = ft_definetrial(cfg);
    
    % get it
    data                        = ft_preprocessing(cfg);
    
    % downsample
    cfg.resamplefs              = 300;
    data                        = ft_resampledata(cfg, data); % resample the data
    
    % fix sample info and add trialnumbers for later 
    data = fixsampleinfo(data);    
    data = rmfield(data, 'cfg');   
    data.trialnumbers = (1:240)';
        
    % get baseline if baseline
    if strcmp(saveName{s},'dataBase')
        nTrials = length(data.trial);
        nChannels = length(data.label);
        meanBase = zeros(nTrials,nChannels);
        for t = 1:nTrials % get the mean over time per trial
           meanBase(t,:) = mean(data.trial{t},2);            
        end       
    end
    
    % remove mean baseline from data
    fprintf('\t Removing the baseline from %s \n',saveName{s})
    nTrials = length(data.trial);
    nSamples = size(data.trial{1},2);
    for t = 1:nTrials        
        for sn = 1:nSamples
        data.trial{t}(:,sn) = data.trial{t}(:,sn)-meanBase(t,:)';
        end
    end

    % save and clean up
    save(fullfile(saveDir,saveName{s}),'data','-v7.3')
    clear cfg data
    
end
