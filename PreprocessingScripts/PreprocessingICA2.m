function PreprocessingICA2(subjectdata)
% function PreprocessingICA(subjectdata)

% instead of running ICA for the three segments seperately, this version
% appends the data prior to the ICA in order to improve the estimation of
% the components. After the components are removed, the data is seperated
% again

% set defaults
ft_defaults

% cd to root
root = '/vol/ccnlab-scratch1/naddij/ImaMEG';
cd(root)

% some settings
outputComp          = fullfile(root,subjectdata.outputDir,'ICAData');
outputData          = fullfile(root,subjectdata.outputDir,'CleanData');
if ~exist(outputComp,'dir'); mkdir(outputComp); end
if ~exist(outputData,'dir'); mkdir(outputData); end
VARData             = fullfile(root,subjectdata.outputDir,'VARData');

%% Get the three data-segments
saveName                = cell(3,1);
saveName{1} = 'dataFP'; saveName{2} = 'dataSP';
saveName{3} = 'dataIM';


% load the VA removed data
load(fullfile(VARData,saveName{1,1})); dataFP = data; clear data
load(fullfile(VARData,saveName{2,1})); dataSP = data; clear data
load(fullfile(VARData,saveName{3,1})); dataIM = data; clear data

% append the data from the three segments
cfg                 = [];
appData             = ft_appenddata(cfg,dataFP,dataSP,dataIM);


% check if ICA already done
if ~exist(fullfile(outputComp,'comp.mat'),'file')
    
    % perform the independent component analysis (i.e., decompose the data)
    cfg                 = [];
    cfg.channel         = 'MEG';
    cfg.method          = 'runica';
    cfg.demean          = 'no';
    comp                = ft_componentanalysis(cfg,appData);
    
    % save the components
    save(fullfile(outputComp,'comp'),'comp','-v7.3')
    
elseif ~exist(fullfile(outputData,[saveName{1} '.mat']),'file') % identify EOG and ECG components
    load(fullfile(outputComp,'comp.mat'))
    
    % correlate to EEG
    EEG = cell2mat(appData.trial);
    EEG = EEG(ismember(appData.label, {'EEG057', 'EEG058', 'EEG059'}), :);
    
    r = corr(cell2mat(comp.trial)', EEG');
    [ro, i] = sort(abs(r),'descend');
    
    fprintf('Highest correlations: \n \t EEG057: comp %d [%.4f] \n \t EEG058: comp %d [%.4f] \n \t EEG059: comp %d [%.4f] \n',i(1,1),r(i(1,1),1),i(1,2),r(i(1,2),2),i(1,3),r(i(1,3),3))
    
    if ~isempty(find(ro(1,:)<0.3,1)) % if some are below 0.3
        % manually check the components
        warning('low correlations, manually check components!')
        return
    end
    
    % inspect these components
    figure;
    cfg                = [];
    cfg.component      = i(1,:);
    cfg.layout         = 'CTF275.lay';
    cfg.commment       = 'no';
    ft_topoplotIC(cfg,comp)
    
    % plot the time course
    tmp_comp = cell2mat(comp.trial);
    figure;
    nComps = length(cfg.component);
    for c = 1:nComps
        subplot(nComps,1,c);
        plot(tmp_comp(cfg.component(1,c),1:2000))
        title(sprintf('Component %d \n',cfg.component(1,c)))
    end
    
    % decide which components to remove and save decision
    comp_removed        = cfg.component;
    save(fullfile(outputComp,'comp'),'comp_removed','-append')
    
    % remove them from the data
    cfg                 = [];
    cfg.component       = comp_removed;
    cfg.demean          = 'no';
    data                = ft_rejectcomponent(cfg, comp, appData);
    
    % seperate the segments again
    cfg                 = [];
    cfg.channels        = 'MEG';
    
    cfg.trials          = ismember(data.trialinfo,100:120);
    dataFP              = ft_selectdata(cfg,data);
    
    cfg.trials          = ismember(data.trialinfo,200:220);
    dataSP              = ft_selectdata(cfg,data);
    
    cfg.trials          = ismember(data.trialinfo,[120:140,220:240]);
    dataIM              = ft_selectdata(cfg,data);
    
    % save the CLEAN data
    data = dataFP;
    save(fullfile(outputData,saveName{1,1}),'data','-v7.3')
    data = dataSP;
    save(fullfile(outputData,saveName{2,1}),'data','-v7.3')
    data = dataIM;
    save(fullfile(outputData,saveName{3,1}),'data','-v7.3')
    
    clear comp comp_removed
else
    fprintf('\n Components already removed for segment %d, clean data saved \n',s)
end
clear data


