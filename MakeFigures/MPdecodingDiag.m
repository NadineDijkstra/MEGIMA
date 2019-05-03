function MPdecodingDiag(cfg)
% function MPdecodingWith(cfg)
%
% cfg.dataName = name of the group effect
% cfg.dataDir  = dir of the the data file
% cfg.time     = time vector
% cfg.ylim     = ylim
% cfg.fillcolor= 'color of the shaded area

% load the data
load(fullfile(cfg.dataDir,cfg.dataName))

% to ease processing
acctmp = zeros(size(acc,3),size(acc,2));
for sub = 1:size(acc,3)
    acctmp(sub,:) = diag(squeeze(acc(:,:,sub)));
end
acc = squeeze(acctmp); clear acctmp
sig = diag(clusterPvalsA) < 0.05;
time = cfg.time{1};
mAcc = mean(acc,3);

% get the SEM
n = size(acc,2);
SEM = std(acc,[],2)/sqrt(n); 
z1 = mAcc + SEM;
z2 = mAcc - SEM;

% plotting main structure
figure; 
plot(time,mAcc,time,z1,time,z2);    
ylim(cfg.ylim); xlim([time(1) time(end-5)]); 

% add SEM shaded area  
a1 = area(time,z1,cfg.ylim(1));
hold on;
set(a1,'LineStyle','none');     set(a1,'FaceColor',cfg.fillcolor);
alpha 0.5
a2 = area(time,z2,cfg.ylim(1)); 
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
hold on; plot(time,mAcc,'black');

% add chance level line
hold on; plot([time(1) time(end)],[0.5 0.5],'black:')

if ~isempty(find(sig,1))
% make line thicker where significant
hold on; plot(time(sig),mAcc(sig),'LineWidth',4,'Color','black')

% color in significant area
hold on; a3 = area(time(sig),mAcc(sig),cfg.ylim(1));
set(a3,'LineStyle','none');     set(a3,'FaceColor',cfg.fillcolor);
border = zeros(length(time),1)+0.5; border(z2<0.5) = z2(z2<0.5);
hold on; a4 = area(time,border,cfg.ylim(1)); 
set(a4,'LineStyle','none');     set(a4,'FaceColor',[1 1 1]);
end

% fix axes
xlim([time(1) time(end-5)]); ylim(cfg.ylim)