function SourcerecAtlas(cfg0)
% function SourcerecAtlas(subjectdata)
%
% applies an atlas to the reconstructed source activity for later averaging
% per atlas region

root = cfg0.root;
subjectID = cfg0.subjectID;
sourcerecDir = cfg0.sourceRecDir;

% change this to own freesurfer folder
addpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));


%% Get the atlas info

[verticesL,labelL,colortable] = read_annotation(fullfile(root,'FreesurferOutput', subjectID, '/label','lh.aparc.a2009s.annot'));
[verticesR,labelR,~]          = read_annotation(fullfile(root,'FreesurferOutput', subjectID, '/label','rh.aparc.a2009s.annot')); %Color tables do not differ with hemisphere

% incase lateralization is not an issue
vertices   = [verticesL; verticesR];
label      = [labelL; labelR];

% read gridpoints
load(fullfile(sourcerecDir,'sourcespace'));

% select the grid points and corresponding label values that are used in
% the MNE grid
vertInUse  = vertices(logical(sourcespace.orig.inuse), :);
labInUse   = label(logical(sourcespace.orig.inuse), :);

inUseL     = sourcespace.orig.inuse(1:length(verticesL), :);
vertInUseL = verticesL(logical(inUseL), :);
labInUseL  = labelL(logical(inUseL), :);

inUseR     = sourcespace.orig.inuse(length(verticesL) +1:end, :);
vertInUseR = verticesR(logical(inUseR), :);
labInUseR  = labelR(logical(inUseR), :);

labs       = colortable.table(:,5);

% label names per grid
labelList = cell(length(vertInUse),1);
labIdx = zeros(length(vertInUse),1);
for i = 1:length(labelList)
    temp = find(labInUse(i) == labs);
    if isempty(temp)
        labIdx(i) = 0;
        labelList{i} = 'NoLabelAssigned';
    else
        labIdx(i) = temp;
        labelList{i} = colortable.struct_names{temp};
    end
end

% label names per grid left
labelListL = cell(length(vertInUseL),1);
labIdxL = zeros(length(vertInUseL),1);
for i = 1:length(labelListL)
    temp = find(labInUseL(i) == labs);
    if isempty(temp)
        labIdxL(i) = 0;
        labelListL{i} = 'NoLabelAssigned';
    else
        labIdxL(i) = temp;
        labelListL{i} = colortable.struct_names{temp};
    end
end

% label names per grid right
labelListR = cell(length(vertInUseR),1);
labIdxR = zeros(length(vertInUseR),1);
for i = 1:length(labelListR)
    temp = find(labInUseR(i) == labs);
    if isempty(temp)
        labIdxR(i) = 0;
        labelListR{i} = 'NoLabelAssigned';
    else
        labIdxR(i) = temp;
        labelListR{i} = colortable.struct_names{temp};
    end
end

cntL = 0;
cntR = 0;
cntCrossHemi = 0;

for i = 1:length(sourcespace.tri)
    % select only those trials that have no interhemispherical
    % connections. This means only throwing out CC (which is
    % labeled as unknown anyway), so no info is really lost
    if sourcespace.tri(i,1) <= length(vertInUseL) && sourcespace.tri(i,2) <= length(vertInUseL) && sourcespace.tri(i,3) <= length(vertInUseL)
        cntL = cntL+1;
        triL(cntL,:) = sourcespace.tri(i,:);
    elseif sourcespace.tri(i,1) > length(vertInUseL) && sourcespace.tri(i,2) > length(vertInUseL) && sourcespace.tri(i,3) > length(vertInUseL)
        cntR = cntR + 1;
        triR(cntR,:) = sourcespace.tri(i,:) - length(vertInUseL); % because the point count starts again at 1 for right! So there is no point 4099 (as noted in tri), that is in the pnt data represented as pnt 1!
    else % possible cross-hemispheric triangles
        cntCrossHemi = cntCrossHemi + 1;
        triCrossHemi(cntCrossHemi, :) = sourcespace.tri(i,:);
    end
end

% Check whether the labels seem to make sense
% bnd.pnt = sourcespace.pos(1:length(vertInUseL), :);
% bnd.tri = triL;
% m = labIdxL;
% ft_plot_mesh(bnd, 'vertexcolor', m);
% 
% hold on
% 
% bnd.pnt = sourcespace.pos(length(vertInUseL)+1:end, :);
% bnd.tri = triR;
% m = labIdxR + max(labIdxL); %to give it another color than left
% ft_plot_mesh(bnd, 'vertexcolor', m);

% saveas(figure(1),fullfile(sourcerecDir,'atlasRegions.png'))
%
% close all

% % plot one region
% region = 1;
% bnd.pnt = sourcespace.pos(1:length(vertInUseL), :);
% bnd.tri = triL;
% m = double(labIdxL==region);
% ft_plot_mesh(bnd, 'vertexcolor', m); view([180 0]);
% title(colortable.struct_names{region})

%% Extract activation per atlas region

% get the source reconstructed activity
if strcmp(cfg0.diff,'yes')
    load(fullfile(sourcerecDir,cfg0.dataSet),'sourceDiff'); source = sourceDiff;
else; load(fullfile(sourcerecDir,cfg0.dataSet),'sourceFace','sourceHouse');
    source = sourceFace;
end

% check whether correct source reconstruction
if ~isfield(source.avg,'pow2'); error('\t Wrong source reconstruction for subject %s dataset %s! \n',...
        subjectID, cfg0.dataSet); end

% first apply depth and noise correction
if strcmp(cfg0.diff,'yes')
    source.avg.pow = (source.avg.pow - source.mPerm)./source.mPerm;
else; sourceFace.avg.pow = (sourceFace.avg.pow - sourceFace.mPerm)./sourceFace.mPerm;
    sourceHouse.avg.pow = (sourceHouse.avg.pow - sourceHouse.mPerm)./sourceHouse.mPerm;
    % average over conditions
    source.avg.pow = (sourceFace.avg.pow+sourceHouse.avg.pow)./2;
end

% left-sided areas
actPerAreaL = zeros(max(labIdxL),size(source.avg.pow,2));
for i = 1:max(labIdxL) % for every label, get the averaged activity. Exclude 0. This is something I gave to no label assigned, so this is not informative anyway
    % get grid points with a given label
    selGrids = (labIdxL==i);
    selGrids = logical([selGrids; zeros(length(labIdxR), 1)]); % pad, because the first half is the left side.
    actPerAreaL(i,:) = nanmean(source.avg.pow(selGrids,:),1);
end
clear selGrids

% right-sided areas
actPerAreaR = zeros(max(labIdxL),size(source.avg.pow,2));
for i = 1:max(labIdxR) % for every label, get the averaged activity.  Exclude 0. This is something I gave to no label assigned, so this is not informative anyway
    % get grid points with a given label
    selGrids = (labIdxR==i);
    selGrids = logical([zeros(length(labIdxL), 1); selGrids]); % pad, because the first half is the left side.
    actPerAreaR(i,:) =  nanmean(source.avg.pow(selGrids,:),1);
end
clear selGrids

% save the activations per area
if strcmp(cfg0.diff,'yes')
    save(fullfile(sourcerecDir,[cfg0.dataSet '_atlasRegions.mat']),'actPerAreaL','actPerAreaR','labIdxL','labIdxR','triL','triR');
else; save(fullfile(sourcerecDir,[cfg0.dataSet,'_avg_atlasRegions.mat']),'actPerAreaL','actPerAreaR','labIdxL','labIdxR','triL','triR');
end

% visualise activity over time
if strcmp(cfg0.plot,'yes')
    sourceActAtlas = source;
    
    % for the left side
    for a = 1:length(labIdxL)
        if labIdxL(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
            sourceActAtlas.avg.pow(a,:) = NaN;
        else
            sourceActAtlas.avg.pow(a,:) = actPerAreaL(labIdxL(a), :);
        end
    end
        
    % for the right side
    for a = 1:length(labIdxR)
        if labIdxR(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
            sourceActAtlas.avg.pow(a+length(vertInUseL),:) = NaN;
        else
            sourceActAtlas.avg.pow(a+length(vertInUseL),:) = actPerAreaR(labIdxR(a), :);
        end
    end    
    
    cfg                = [];
    cfg.maskparameter  = [];
    cfg.funparameter   = 'pow';
    ft_sourcemovie(cfg,sourceActAtlas);
end

