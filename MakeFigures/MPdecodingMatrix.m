function MPdecodingMatrix(cfg)
% function MPdecodingMatrix(cfg)
%
% cfg.dataName = name of the group effect
% cfg.dataDir  = dir of the the data file
% cfg.time = 2x1 cell with time vectors for x and y axes
% cfg.clim = limits clim


% load the data
load(fullfile(cfg.dataDir,cfg.dataName))

% to ease processing
time = cfg.time;
if numel(time) == 1; time{2} = time{1}; end
sigA = clusterPvalsA < 0.05;
mAcc = squeeze(mean(acc,3));

% discretize signal for nicer plots
bins = [cfg.clim(1):cfg.steps:cfg.clim(2)];
nBins = length(bins);
mAccBin = mAcc;
for b = 1:nBins-1
   mAccBin(mAcc >= bins(b) & mAcc < bins(b+1)) = (bins(b)+bins(b+1))/2;
    
end

% do the main plotting
figure;
imagesc(time{2},time{1},mAccBin); hold on; 
axis xy; axis image; xlabel('Time (s)'); ylabel('Time (s)'); 

% add contour for significance 
contour(time{2},time{1},sigA','color','k','LineWidth',1)

% change colorbar
h = colorbar; ylabel(h, 'Accuracy'); 
caxis(cfg.clim)

% fix axes
xlim([time{2}(1) time{2}(end-5)])
ylim([time{1}(1) time{1}(end-5)])

% change the colormap 
map = makeColorMaps('redblue');
colormap(flipud(map));