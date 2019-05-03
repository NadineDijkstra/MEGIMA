function SourcerecConDiff(cfg0)
% function SourcerecConDiff(cfg)

root = '/vol/ccnlab-scratch1/naddij/ImaMEG';
sourceRecDir = fullfile(root,cfg0.subjectID,'SourceReconstruction');

% fill out the defaults
if ~isfield(cfg0,'time'); cfg0.time = []; end
if ~isfield(cfg0,'lpfilter'); cfg0.lpfilter = 'no'; end

%% Prepare the MEG data

% get the data
load(fullfile(root,cfg0.subjectID,'CleanData',cfg0.dataName));

% average over time for imagery
if ~isempty(cfg0.time)
    fprintf('Averaging over time window %.2f to %.2f \n',cfg0.time(1),cfg0.time(2))
    cfg         = [];
    cfg.latency = cfg0.time;
    cfg.avgovertime = 'yes';
    data = ft_selectdata(cfg,data);
elseif strcmp(cfg0.lpfilter,'yes')
    cfg          = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 30;
    data = ft_preprocessing(cfg,data);
end

% timelock analysis
cfg  = [];
cfg.channel    = {'MEG'};
cfg.keeptrials = 'yes';
tlckAll        = ft_timelockanalysis(cfg,data);

%% Forward solution
load(fullfile(sourceRecDir,'headmodel'),'vol')
load(fullfile(sourceRecDir,'sourcespace'))
vol    = ft_convert_units(vol,'mm');

cfg           = [];
cfg.grad      = tlckAll.grad;
cfg.grid.pos  = sourcespace.pos;
cfg.grid.inside = 1:size(sourcespace.pos,1);
if strcmp(cfg0.subjectID,'S15'); cfg.channel = {'MEG','-MLT31'};
else; cfg.channel   = {'MEG'}; end
cfg.headmodel = vol;
cfg.backproject = 'no';
leadfield     = ft_prepare_leadfield(cfg);

%% Inverse solution
cfg = [];
cfg.lambda = 0.01;                 % Regularization; 0.01 corresponds to '1%' in FieldTrip
cfg.numPerm = 1000;                 % 100 is fine for testing, but you want more for the real analysis
cfg.feedback = 2;                  % provide feedback (soonest) every X seconds
cfg.allowNan = 'no';              % Use nanmean, nancov, etc.?
cfg.lf = leadfield;
cfg.cov.method = 'fixed';
if strcmp(cfg0.dataName,'dataIM'); cfg.cov.window = 'all';
elseif strcmp(cfg0.dataName,'dataFP'); cfg.cov.window = [0.05 1]; end
cfg.returnSqrt = 'no';
cfg.locationRescale = 'no';

design = ismember(data.trialinfo,cfg0.conIdx{1});
sourceDiff = sourceAnalysisTwoConditions(cfg, tlckAll, design);
sourceDiff.tri = sourcespace.tri;

% correct for positivity and depth bias
sourceDiff.avg.pow2 = (sourceDiff.avg.pow - sourceDiff.mPerm)./sourceDiff.mPerm;
sourceDiff.avg.pow2(sourceDiff.avg.pow2 < 0) = 0;
sourceDiff.avg.pow2 = sqrt(sourceDiff.avg.pow2);
if ~isempty(cfg0.time)
    save(fullfile(sourceRecDir,sprintf('%s_avg_%.2fto%.2f.mat',cfg0.dataName,cfg0.time(1),cfg0.time(2))),'sourceDiff','cfg0')
else
    save(fullfile(sourceRecDir,cfg0.dataName),'sourceDiff','cfg0')
end


%% Visualization
if isfield(cfg0,'plot')
    bnd.pnt = sourceDiff.pos;
    bnd.tri = sourceDiff.tri;
    sample = find(sourceDiff.time > 0.14);
    m = sourceDiff.avg.pow2(:,sample(1));
    
    figure; ft_plot_mesh(bnd,'vertexcolor',m); colorbar;
    title(cfg0.dataName); view([-90 0])
    pause(3); close
end

if strcmp(cfg0.movie,'yes')
    cfg                = [];
    cfg.funparameter   = 'pow2';
    ft_sourcemovie(cfg,sourceDiff);
end

