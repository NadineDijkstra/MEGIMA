function BEHanalysis(subjectdata)
% function BEHanalysis(subjectdata)

% cd to root
root = '/vol/ccnlab-scratch1/naddij/ImaMEG';
cd(root)

% outputDir
outputDir = fullfile(root,subjectdata.outputDir,'Behaviour');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% get the behavioural data
[E,B] = ProcessBehaviour(subjectdata);

vividness = B(:,1);

save(fullfile(outputDir,'vividness'),'vividness')
save(fullfile(outputDir,'responses'),'E','B')