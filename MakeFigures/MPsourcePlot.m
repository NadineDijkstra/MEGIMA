function MPsourcePlot(cfg)
% function MPsourcePlot(cfg)
%
% cfg.dataName = name of the group effect
% cfg.dataDir  = dir of the the data file
% cfg.time     = time point to plot activity of 
% cfg.clim     = limits of the color axis

addpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));

% load the data 
load(fullfile(cfg.dataDir,cfg.dataName))

% average over participants
actLeft = squeeze(nanmean(actLeft,3));
actRight = squeeze(nanmean(actRight,3));
actLeft(actLeft < 0) = 0; actLeft = sqrt(actLeft);
actRight(actRight < 0) = 0; actRight = sqrt(actRight);

% to ease processing
time = sourceActAtlas.time;

% get area labels
[~,~,colortable] = read_annotation(fullfile(fileparts(fileparts(cfg.dataDir)),'FreesurferOutput', 'S01', '/label','lh.aparc.a2009s.annot'));

% find samples
samples = find(round(time,3) == cfg.time);
if isempty(samples) % if there is no exact time point, select the ones around it
tmp = find(time < cfg.time); samples(1) = tmp(end);
tmp = find(time > cfg.time); samples(2) = tmp(1);
end

% select data from that time point
source = sourceActAtlas;
source.avg.pow = mean(source.avg.pow(:,samples),2);
source = rmfield(source,{'time','mPerm'});
source.avg = rmfield(source.avg,'pow2');
actLeft = mean(actLeft(:,samples),2); actRight = mean(actRight(:,samples),2);

% plot  it
figure;
bnd.pnt = source.pos;
bnd.tri = source.tri; 
ft_plot_mesh(bnd,'vertexcolor',source.avg.pow,'edgecolor','none');%,'facecolor',[150,150,150]/255);
view([-25 -10])
camlight;
lighting gouraud;
material dull;
% change the colormap 
map = makeColorMaps('redblue');
colormap(flipud(map)); colorbar
if ~isempty(cfg.clim); caxis(cfg.clim); end


% show the names of active areas
leftSig = find(actLeft > 0);
nAreas = numel(leftSig);
fprintf('Active areas left are: \n')
for a = 1:nAreas    
   fprintf('\t %s \n',colortable.struct_names{leftSig(a)})    
end

rightSig = find(actRight > 0);
nAreas = numel(rightSig);
fprintf('Active areas right are: \n')
for a = 1:nAreas    
   fprintf('\t %s \n',colortable.struct_names{rightSig(a)})    
end