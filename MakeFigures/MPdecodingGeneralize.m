function MPdecodingGeneralize(cfg)
% function MPdecodingGeneralize(cfg)
%
% cfg.dataName = name of the group effect
% cfg.dataDir  = dir of the the data file
% cfg.time     = 2x1 cell with time vectors for x and y axes
% cfg.ylim     = ylim
% cfg.fillcolor= 'color of the shaded area
% cfg.timePoint= time of train point to plot accuracy of

% load the data
load(fullfile(cfg.dataDir,cfg.dataName))

% to ease processing
time = cfg.time;
if numel(time) == 1; time{2} = time{1}; end

% find samples
samples = find(round(time{1},3) == cfg.timePoint);
if isempty(samples) % if there is no exact time point, select the ones around it
tmp = find(time{1} < cfg.timePoint); samples(1) = tmp(end);
tmp = find(time{1} > cfg.timePoint); samples(2) = tmp(1);
end

% select accuracy for that sample
tmp = acc;
acc = squeeze(mean(acc(samples,:,:),1));
mAcc = mean(acc,2);

% select significance of that samples
clusterPvalsA = clusterPvalsA';
if numel(samples) > 1 % both time points around need to be significant
    sig = clusterPvalsA(samples(1),:) < 0.05 & clusterPvalsA(samples(1),:) < 0.05;
else
sig = clusterPvalsA(samples,:) < 0.05;
end

% plot everything
% get the SEM
n = size(acc,2);
SEM = std(acc,[],2)/sqrt(n); 
z1 = mAcc + SEM;
z2 = mAcc - SEM;

% plotting main structure
figure; 
plot(time{2},mAcc,time{2},z1,time{2},z2);    
ylim(cfg.ylim); xlim([time{2}(1) time{2}(end-5)]); 

% add SEM shaded area  
a1 = area(time{2},z1,cfg.ylim(1));
hold on;
set(a1,'LineStyle','none');     set(a1,'FaceColor',cfg.fillcolor);
alpha 0.5
a2 = area(time{2},z2,cfg.ylim(1)); 
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
hold on; plot(time{2},mAcc,'black');
    
if ~isempty(find(sig,1))
    % get significant clusters
    % determine the clusters for the accuracy plot
    s = 1; c= 0;
    while s < length(sig)
        while sig(s) == 0 && s < length(sig)% look for a 1
            s = s + 1;
        end
        c = c + 1;
        clusterA(c,1) = s;
        while sig(s) == 1 && s < length(sig)% keep running through the cluster until end
            s = s + 1;
            if s > length(sig); s = s -1; end
        end
        clusterA(c,2) = s;
        
        % remove clusters of only 1 datapoint
        if clusterA(c,2) == clusterA(c,1)+1; clusterA(c,:) = []; end
    end
    x = clusterA(:,1)==length(sig); % delete extra row
    clusterA(x,:) = []; clear x;
    clusterA(clusterA(:,1)==0,:) = []; nClustersA = size(clusterA,1);
    
    % run over significant clusters to make thicker line and color in
    for c = 1:nClustersA
        % make line thicker where significant
        hold on; plot(time{2}(clusterA(c,1):clusterA(c,2)),mAcc(clusterA(c,1):clusterA(c,2)),'LineWidth',4,'Color','black')
        
        % for positive clusters
        if mAcc(clusterA(c,1):clusterA(c,2)) > cfg.chance
            % color in significant area
            hold on; a3 = area(time{2}(clusterA(c,1):clusterA(c,2)),mAcc(clusterA(c,1):clusterA(c,2)),cfg.ylim(1));
            set(a3,'LineStyle','none');     set(a3,'FaceColor',cfg.fillcolor);
            
            % whiten below chance
            border = zeros(length(time{2}),1)+cfg.chance; border(z2<cfg.chance) = z2(z2<cfg.chance);
            hold on; a4 = area(time{2},border,cfg.ylim(1));
            set(a4,'LineStyle','none');     set(a4,'FaceColor',[1 1 1]);
            clear a3 a4
        else % for negative clusters
            % color in significant area
            border = zeros(length(time{2}))+cfg.chance;
            lower  = mAcc(clusterA(c,1):clusterA(c,2));
            upper  = border(clusterA(c,1):clusterA(c,2));
            
            hold on; ciplot(lower,upper,time{2}(clusterA(c,1):clusterA(c,2)),'red')
            
        end
    end
end

% add chance level line
hold on; plot([time{2}(1) time{2}(end)],[0.5 0.5],'black:')


% fix axes
xlim([time{2}(1) time{2}(end-5)]); ylim(cfg.ylim)
