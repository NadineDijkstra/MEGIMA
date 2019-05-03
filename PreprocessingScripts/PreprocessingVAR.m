function PreprocessingVAR(subjectdata)
% function PreprocessingVAR(subjectdata)
% Visual artefact rejection as part of the preprocessing progress. 

% set defaults
ft_defaults

% cd to root
root = '/vol/ccnlab-scratch1/naddij/ImaMEG';
cd(root)

% some settings
outputDir               = fullfile(root,subjectdata.outputDir,'VARData');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
PreprocData             = fullfile(root,subjectdata.outputDir,'PreprocData');

%% Loop over data segments
dataCell                = str2fullfile(PreprocData,'data*.mat');
nSeg                    = numel(dataCell);
saveName                = cell(nSeg,2);
saveName{1,1} = 'artBase'; saveName{2,1} = 'artFP';
saveName{3,1} = 'artIM'; saveName{4,1} = 'artSP';
saveName{1,2} = 'dataBase'; saveName{2,2} = 'dataFP';
saveName{3,2} = 'dataIM'; saveName{4,2} = 'dataSP';
stimOn                  = [0 0;0 0.8; 0 4;0 0.8];

for s = 1:nSeg 

% load the data
load(dataCell{s})

% get some basic values
time                    = data.time{1};
nTrials                 = length(data.trial);

% select only the MEG channels
cfg                     = [];
cfg.channel             = 'MEG';
tmp_data                = ft_selectdata(cfg,data);
tmp_data.trialinfo      = [tmp_data.trialinfo, (1:nTrials)'];  

% Overall artifacts
cfg                     = [];
cfg.method              = 'summary';
tmp_data_overall        = ft_rejectvisual(cfg, tmp_data); % save them
removed_n_overall       = setdiff(1:nTrials, tmp_data_overall.trialinfo(:, end));
clear tmp_data_overall;

% Muscle artifacts
iNan = cell(nTrials);
for in = 1:length(tmp_data.trial)
    iNan{in} = isnan(tmp_data.trial{in});
    tmp_data.trial{in}(iNan{in}) = 0;
end

cfg                     = [];
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 100;
tmp_data_filtered       = ft_preprocessing(cfg, tmp_data);

for in = 1:length(tmp_data.trial)
    tmp_data_filtered.trial{in}(iNan{in}) = nan;
end

tmp_data_muscle         = ft_rejectvisual([], tmp_data_filtered);
removed_n_muscle        = setdiff(1:nTrials, tmp_data_muscle.trialinfo(:, end));
clear tmp_data_muscle tmp_data_filtered;

% Blinks during stimulus
cfg                     = [];
cfg.channel             = {'EEG057', 'UADC007'};
tmp_data                = ft_selectdata(cfg, data);
tmp_data.trialinfo      = [tmp_data.trialinfo, (1:nTrials)'];

X = reshape(cell2mat(tmp_data.trial), [2, length(tmp_data.time{1}), length(tmp_data.trial)]);
X_EOGv = squeeze(X(1, :, :))';
X_EOGv = (X_EOGv - mean(X_EOGv(:))) ./ std(X_EOGv(:));
X_pupDil = squeeze(X(2, :, :))';
X_pupDil = (X_pupDil - mean(X_pupDil(:))) ./ std(X_pupDil(:));

cutoff_z_EOGv = 5;
cutoff_z_pupDil = 3.5;

blinkMask               = (abs(X_EOGv) > cutoff_z_EOGv) | (abs(X_pupDil) > cutoff_z_pupDil);
tmp                     = zeros(1,length(time)); tmp(1,find(time==stimOn(s,1)):find(time==stimOn(s,2)))=1;
stimOnMask              = repmat(tmp,[nTrials,1]);
removed_n_blinks        = tmp_data.trialinfo(any(stimOnMask & blinkMask, 2), end);

clear tmp_data blinkMask stimOnMask

% Inspect potentially contaminated trials
cfg                     = [];
cfg.channel             = {'MEG', 'EEG057', 'EEG058'};
cfg.megscale = 1e8;
cfg.artfctdef.overall.artifact = data.sampleinfo(removed_n_overall, :);
cfg.artfctdef.muscle.artifact = data.sampleinfo(removed_n_muscle, :);
cfg.artfctdef.blinks.artifact = data.sampleinfo(removed_n_blinks, :);
clear removed_n_overall removed_n_muscle removed_n_blinks
cfg.preproc.hpfilter    = 'no';
cfg.preproc.hpfreq      = 100;
cfg.renderer            = 'painters';
cfgart                  = ft_databrowser(cfg, data);

% Get baseline artifacts
if s == 1
    artifactsBase = cat(1,cfgart.artfctdef.visual.artifact,cfgart.artfctdef.overall.artifact,cfgart.artfctdef.muscle.artifact,cfgart.artfctdef.blinks.artifact);
    nArts = length(artifactsBase);
    baseArts = zeros(nArts,1);
    for a = 1:nArts % begin artifact later than begin-sample and earlier than end-sample
       baseArts(a,1) = find(artifactsBase(a,1) >= data.sampleinfo(:,1) & artifactsBase(a,1) <= data.sampleinfo(:,2));
    end
    baseArts = unique(baseArts);
elseif s > 1
    artifactsBase = data.sampleinfo(baseArts,:); % get the sample data
    cfgart.artfctdef.baseline.artifact = artifactsBase;
end
clear nArts

% Save the artifacts
save(fullfile(outputDir,saveName{s,1}),'cfgart')

% Reject the artifacts
cfgart.artfctdef.reject = 'complete';
data                    = ft_rejectartifact(cfgart, data);

% Save the data and clean up
save(fullfile(outputDir,saveName{s,2}),'data','-v7.3')
clear data cfgart  

end

