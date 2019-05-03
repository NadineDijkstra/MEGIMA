function DecoderVividnessOutput(cfg0)
% function DecoderVividnessOutput(cfg)
%
% Transforms Xhat to evidence for correct class and correlates this with
% vividness over trials
%
% cfg.dataDir = path to subject specific dir
% cfg.dataName = which decoding to use
% cfg.split    = whether to split for the two cues (because of vividness
% differences and overlap differences)

[~,name] = fileparts(cfg0.dataName);

outputDir = fullfile(cfg0.dataDir,'Vividness');

% get the decoding output
load(fullfile(cfg0.dataDir,cfg0.dataName),'Xhat','trueClass','trialnumbers','var','idx')

% get the vividness ratings
load(fullfile(cfg0.dataDir,'Behaviour','vividness'))

% get the trialMatrix info
load(fullfile(cfg0.dataDir,'Behaviour','trialMatrix'))

% change to correct-class activation
Evidence  = Xhat;
Evidence(:, :, trueClass==0) = Evidence(:, :, trueClass==0)*-1;
Evidence = squeeze(Evidence);

% correlate to vividness rating
if ~exist('idx','var'); if numel(trialnumbers) < 2; idx = trialnumbers{1};
    else; idx = trialnumbers{2}; end; end
viv  = vividness(idx);

r    = zeros(size(Xhat,1),size(Xhat,2));
pval = zeros(size(Xhat,1),size(Xhat,2));

for s = 1:size(Xhat,1)
    [r(s,:),pval(s,:)] = corr(squeeze(Evidence(s,:,:))',viv);
end

% save the things
save(fullfile(outputDir,[name '_corr.mat']),'r')

