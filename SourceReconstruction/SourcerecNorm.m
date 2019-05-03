function SourcerecNorm(cfg0,subjects)
% function SourcerecNorm(cfg,subjects)
%
%  subjects = cell array of subject names to average over
%  cfg.dataName = which data to average
%  cfg.movie = 'yes' or 'no' (default = 'no')
%  cfg.ROIs  = areas to plot the time course of in numbers

% add your own freesurfer path here!
addpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));

% fill out default
if ~isfield(cfg0,'movie'); cfg0.movie = 'no'; end

root =  cfg0.root;
outputDir = fullfile(root,'GroupResults','SourceReconstruction');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% use S01's brain for plotting
load(fullfile(root,'S01','SourceReconstruction',cfg0.dataName),'labIdxL','labIdxR','triL','triR')
load(fullfile(root,'S01','SourceReconstruction',cfg0.dataName(1:6)),'sourceDiff')
% get area labels
[~,~,colortable] = read_annotation(fullfile(root,'FreesurferOutput', 'S01', '/label','lh.aparc.a2009s.annot'));

% load data for all subjects
nSubjects = length(subjects);
if ~exist(fullfile(outputDir,cfg0.dataName),'file')
for sub = 1:nSubjects    
    load(fullfile(root,subjects{sub},'SourceReconstruction',cfg0.dataName),'actPerAreaL','actPerAreaR');   
    if sub == 1
        actLeft = actPerAreaL;
        actRight = actPerAreaR;
    else 
        actLeft = cat(3,actLeft,actPerAreaL);
        actRight = cat(3,actRight,actPerAreaR);
    end  
    clear actPerAreaL actPerAreaR
end
save(fullfile(outputDir,cfg0.dataName), 'actLeft', 'actRight')
else
    load(fullfile(outputDir,cfg0.dataName))    
end

% average over participants
actLeft = squeeze(nanmean(actLeft,3));
actRight = squeeze(nanmean(actRight,3));

% plot one time window
if isfield(cfg0,'time') % if this is the case we first need to average over the appropriate time window before zero clipping and square root
    % get samples
    samples = sourceDiff.time >= cfg0.time(1) & sourceDiff.time <= cfg0.time(2);
    
    % put it into a plotable structure
    mAct = zeros(length(sourceDiff.inside),1);
    for a = 1:length(sourceDiff.inside)
        if a <= length(labIdxL) % left
            mAct(a,:) = mean(actLeft(labIdxL(a), samples),2);
        else % right
            mAct(a,:) = mean(actRight(labIdxR(a-length(labIdxL)), samples),2);
        end
    end
    
    % zero clipping and square root
    mAct(mAct < 0) = 0; mAct = sqrt(mAct);
    
    % plot
    figure(1);   
    % left side
    bnd.pnt = sourceDiff.pos(1:length(labIdxL),:); bnd.tri = triL; subplot(2,2,1);
    ft_plot_mesh(bnd,'vertexcolor',mAct(1:length(labIdxL)),'edgecolor','black');
    view([180 0]); subplot(2,2,3); ft_plot_mesh(bnd,'vertexcolor',mAct(1:length(labIdxL)),'edgecolor','black');
    view([0 0]); subplot(2,2,2); % right side
    bnd.pnt = sourceDiff.pos(length(labIdxL)+1:end,:); bnd.tri = triR;
    ft_plot_mesh(bnd,'vertexcolor',mAct(length(labIdxL)+1:end),'edgecolor','black');
    view([0 0]); subplot(2,2,4); ft_plot_mesh(bnd,'vertexcolor',mAct(length(labIdxL)+1:end),'edgecolor','black');
    view([180 0]); h = colorbar('south'); set(h,'Position',[0.25 0.07 0.5 0.035]);
    subtitle(sprintf('%s Time %.2f to %.2f',cfg0.dataName(1:6),cfg0.time(1),cfg0.time(2)));

    figure(2)
    bnd.pnt = sourceDiff.pos; bnd.tri = sourceDiff.tri;
    ft_plot_mesh(bnd,'vertexcolor',mAct,'edgecolor','black')
    
end

% do zero clip and take the square root
actLeft(actLeft < 0) = 0; actLeft = sqrt(actLeft);
actRight(actRight < 0) = 0; actRight = sqrt(actRight);

% put it into a plotable structure
if ~exist('sourceActAtlas','var')
    sourceActAtlas = sourceDiff;
    sourceActAtlas.avg.pow = [];
    for a = 1:length(sourceDiff.inside)
        if a <= length(labIdxL) % left
            sourceActAtlas.avg.pow(a,:) = actLeft(labIdxL(a), :);
        else % right
            sourceActAtlas.avg.pow(a,:) = actRight(labIdxR(a-length(labIdxL)), :);
        end
    end
    % save
    save(fullfile(outputDir,cfg0.dataName), 'sourceActAtlas', '-append')
end

% if no time dimension
if size(actLeft,2) == 1    
    % plot the data
    bnd.pnt = sourceActAtlas.pos;
    bnd.tri = sourceActAtlas.tri;
    m = sourceActAtlas.avg.pow;
    
    figure; ft_plot_mesh(bnd,'vertexcolor',m,'edgecolor','black'); 
    colorbar;
    title(cfg0.dataName);
else


% plot movie of activation
if strcmp(cfg0.movie,'yes')
cfg                = [];
cfg.maskparameter  = [];
cfg.funparameter   = 'pow';
ft_sourcemovie(cfg,sourceActAtlas);
end

% plot time course of activation
if isfield(cfg0,'ROIs')
figure; nROIs = numel(cfg0.ROIs);
for r = 1:nROIs
    subplot(2,1,1); plot(sourceActAtlas.time,squeeze(mean(actLeft(cfg0.ROIs(r),:,:),3)));
    hold on; title('Activation left'); xlim([sourceActAtlas.time(1) sourceActAtlas.time(end)])

    subplot(2,1,2); plot(sourceActAtlas.time,squeeze(mean(actRight(cfg0.ROIs(r),:,:),3)));
    hold on; title('Activation right')
end
legend(colortable.struct_names{cfg0.ROIs})
xlim([sourceActAtlas.time(1) sourceActAtlas.time(end)])
end
end