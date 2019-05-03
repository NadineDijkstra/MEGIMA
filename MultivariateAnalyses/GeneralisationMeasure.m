function GeneralisationMeasure(cfg)
% function GeneralisationMeasure(cfg)

% load the data
load(fullfile(cfg.dataDir,cfg.dataName))

% to ease processing
time = cfg.time;
if numel(time) == 1; time{2} = time{1}; end
nTime = size(acc,2);
if size(acc,1) == size(acc,2) % train and test same data
sigA  = diag(clusterPvalsA) < 0.05;
else
    sigA = clusterPvalsA < 0.05; sigA = any(sigA);
end

%% King et al., 2014 Plos One

% calculate generalisation per time point as number of significant
% predictions (FDR corrected) during significant time window

% select significant acc.
sigAcc  = acc(sigA,sigA,:);
nTime   = size(sigAcc,1);
    
% test acc off-diagonal to acc diagonal to see whether there is actually
% more info at this time point 
pVals = zeros(nTime,nTime);
for t = 1:nTime
    for T = 1:nTime
        [~,pVals(t,T)] = ttest(squeeze(sigAcc(t,t,:)),squeeze(sigAcc(t,T,:)),'Tail','right');
    end
end

% FDR correct
pVals = reshape(pVals,[1,nTime*nTime]);
h = fdr_bh(pVals,.05,'dep');

% count number of significant points
diagDiff = reshape(h,[nTime,nTime]);


% plot the results
figure;
plot(time{1}(sigA),sum(diagDiff)./nTime,'LineWidth',2,'Color','black'); ylim([0 1])

