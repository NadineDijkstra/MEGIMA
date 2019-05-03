function DecodingGA(cfg0)
% averages over participants and performs a cluster based permutation test
% at the group level using nPerm pemrutations
%
% cfg.dataName = which data to use
% cfg.dirName  = directory of the data per subject
% cfg.root     = root directory
% cfg.nPerm    = how many permutations
% cfg.alpha    = alpha cut-off (default = 0.05)
%

if ~isfield(cfg0,'alpha'); cfg0.alpha = 0.05; end % defaults

outputDir = fullfile(cfg0.root,'GroupResults',cfg0.dirName);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

subjects = cfg0.subjects;
time = cfg0.time;

% determine brain data or eye data
if contains(cfg0.dirName,'eye','IgnoreCase',1)
    eyeData = true;
else 
    eyeData = false;
end
%% Load the data
nsubjects = length(subjects);
for sub = 1:nsubjects
    if eyeData
        load(fullfile(cfg0.root,subjects{sub},cfg0.dirName,cfg0.dataName),'Accuracy')
        if sub > 1
            acc = cat(3,acc,Accuracy);
        else
            acc = Accuracy;
        end        
    else
        load(fullfile(cfg0.root,subjects{sub},cfg0.dirName,cfg0.dataName),'Accuracy','m0','m1')
        if sub > 1
            acc = cat(3,acc,Accuracy);
            M0  = cat(3,M0,m0);
            M1  = cat(3,M1,m1);
        else
            acc = Accuracy;
            M0 = m0;
            M1 = m1;
        end
    end
end

%% Run the permutation tests
% first check whether they have already been done
if ~exist(fullfile(outputDir,[cfg0.dataName '.mat']),'file')
    cfg = [];
    cfg.iterations = cfg0.nPerm;
    cfg.tail       = 'two';
    clusterPvalsA = cluster_based_permutationND(permute(acc,[3,2,1]),0.5,cfg);
    save(fullfile(outputDir,cfg0.dataName),'clusterPvalsA','acc','cfg0')
    if ~eyeData
        clusterPvalsM = cluster_based_permutationND(permute(M0,[3,2,1]),permute(M1,[3,2,1]),cfg);
        save(fullfile(outputDir,cfg0.dataName),'M0','M1','clusterPvalsM','-append')
    end
else
    load(fullfile(outputDir,cfg0.dataName),'clusterPvalsA','clusterPvalsM')
end
sigA = clusterPvalsA < cfg0.alpha;
if ~eyeData
sigM = clusterPvalsM < cfg0.alpha;
end


%% Do the plotting
acc = squeeze(acc); 
if ~eyeData; M0 = squeeze(M0); M1 = squeeze(M1); end
if numel(time) < 2; time{2} = time{1}; end % determine time vectors

if numel(size(acc)) == 2 % timedecoding
       
    if iscell(time); if cfg0.dataName(5) == 'I'; time = time{2}; else time = time{1}; end; end
    
    % determine the clusters for the accuracy plot
    s = 1; c= 0;
    while s < length(sigA)
        while sigA(s) == 0 && s < length(sigA)% look for a 1
            s = s + 1;
        end
        c = c + 1;
        clusterA(c,1) = s;
        while sigA(s) == 1 && s < length(sigA)% keep running through the cluster until end
            s = s + 1;
            if s > length(sigA); s = s -1; end
        end
        clusterA(c,2) = s;
    end
    x = find(clusterA(:,1)==length(sigA)); % delete extra row
    clusterA(x,:) = []; clear x; nClustersA = size(clusterA,1);
    
    save(fullfile(outputDir,cfg0.dataName),'clusterA','-append')
    
    % determine the clusters for the discriminant channels plots
    if ~eyeData
    s = 1; c= 0;
    while s < length(sigM)
        while sigM(s) == 0 && s < length(sigM)% look for a 1
            s = s + 1;
        end
        c = c + 1;
        clusterM(c,1) = s;
        while sigM(s) == 1 && s < length(sigM)% keep running through the cluster until end
            s = s + 1;
            if s > length(sigM); s = s -1; end
        end
        clusterM(c,2) = s;
    end
    x = find(clusterM(:,1)==length(sigM)); % delete extra row
    clusterM(x,:) = []; clear x; nClustersM = size(clusterM,1);
    end
    
    % do the actual plotting
    figure;
    subplot(2,1,1); mAcc = mean(acc,2); % mean accuracy
    plot(time,mAcc,'b'); hold on;
    for c = 1:nClustersA % plot the clusters seperately
        plot(time(clusterA(c,1):clusterA(c,2)),mAcc(clusterA(c,1):clusterA(c,2)),'b','LineWidth',2)
        hold on;
    end
    hold on; plot([time(1) time(end)],[0.5 0.5],'r');
    xlim([time(1) time(end)]); title(cfg0.dataName); xlabel('Time (s)'); ylabel('Accuracy')
    
    if ~eyeData
    subplot(2,1,2); m0 = mean(M0,2); m1 = mean(M1,2); % m1 and m0
    plot(time,m0,'black'); hold on; plot(time,m1,'red')
    hold on; plot(time,m1,'red'); hold on;
    for c = 1:nClustersM % plot the clusters seperately
        plot(time(clusterM(c,1):clusterM(c,2)),m0(clusterM(c,1):clusterM(c,2)),'black','LineWidth',2)
        hold on;
        plot(time(clusterM(c,1):clusterM(c,2)),m1(clusterM(c,1):clusterM(c,2)),'red','LineWidth',2)
        hold on;
    end
    xlim([time(1) time(end)]); xlabel('Time (s)'); ylabel('Discriminant channel')
    end
    
    sigTime = time(sigA); [~,b] = max(mAcc);
    if ~isempty(sigTime); fprintf('\t Start significace at %.3f s after stim onset \n',sigTime(1)); end
    fprintf('\t Peak accuracy at %.3f s after stim onset \n',time(b))

else % cross decoding
    figure;
    if size(acc,2) > size(acc,1); subplot(2,1,1); else; subplot(1,2,1); end
    mAcc = squeeze(mean(acc,3));   
    imagesc(time{2},time{1},mAcc); hold on; contour(time{2},time{1},sigA','color','k','LineWidth',1)
    axis xy; axis image; xlabel('Time (s)'); ylabel('Time (s)'); title(cfg0.dataName)
    h = colorbar; ylabel(h, 'Accuracy'); colormap(jet(256)); caxis([min(min(mAcc)) max(max(mAcc))])
    if size(acc,2) > size(acc,1); subplot(2,1,2); else; subplot(1,2,2); end
    if ~eyeData
    m0 = squeeze(mean(M0,3)); m1 = squeeze(mean(M1,3));
    imagesc(time{2},time{1},m1-m0); hold on; contour(time{2},time{1},sigM','color','k','LineWidth',2)
    axis xy; axis image; xlabel('Time (s)'); ylabel('Time (s)');
    h = colorbar; ylabel(h, 'm1 - m0'); caxis([min(min(m1-m0)) max(max(m1-m0))])
    end
end


