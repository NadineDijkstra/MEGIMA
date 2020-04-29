function makeTrialmatrix(subjectID)
% function makeTrialmatrix(subjectID,nTrials)
% subjectID:    string under which the file will 
%               be saved 'trialMatrix_subjectID'            

rng('shuffle', 'twister');

% create trial matrix
nStim = 16;
nCues = 2;
nProbes = 4;
nBlocks = 10;
class1 = 1:8;
class2 = 9:16;

% get all possible permutations
M = permn(1:nStim,2);

% remove repetitions
i = 1;
while i <= size(M,1);
    if M(i,1) == M(i,2)
        M(i,:) = [];
    end
    i = i+1;   
end
nOptions = size(M,1);

% add cues 
cCue = size(M,2)+1; % cue column
cueVec = 1:nCues;
M(:,cCue) = repmat(cueVec',nOptions/nCues,1);

% add reverse rating column
cReverse = size(M,2)+1;
reverseVec = [0;1];
M(:,cReverse) = repmat(reverseVec,nOptions/2,1);

% add random catch trials
cCatch = size(M,2)+1;
ima = zeros(nOptions,1);
for i = 1:nOptions
    ima(i,1) = M(i,M(i,3)); % identity cued stimulus
end
for i = 1:nStim
    tmp = find(ima == i);
    idx = randi(length(tmp));
    M(tmp(idx),cCatch) = 1;    
end

% create probes for catch trials
cProbes = size(M,2)+1;
M(:,cProbes:cProbes+nProbes-1) = 0;
for i = 1:nOptions
    if M(i,cCatch) == 1
        
       imagined = M(i,M(i,3)); 
       % class1 or class2?
       if ismember(imagined,class1) 
       probes = randperm(numel(class1));
       probes(probes == imagined) = [];
       probes = [imagined probes(1:nProbes-1)];
       probes = Shuffle(probes);
       elseif ismember(imagined,class2) 
       probes = randperm(numel(class1))+numel(class2);
       probes(probes == imagined) = [];
       probes = [imagined probes(1:nProbes-1)];
       probes = Shuffle(probes);       
       end
              
       M(i,cProbes:cProbes+nProbes-1) = probes; 
    end
end

% shuffle 
P = randperm(nOptions);
M = M(P,:);
ima = ima(P,:);

% balance classes per block
idx1 = find(ismember(ima,class1)); idx2 = find(ismember(ima,class2));
nTrialsPerBlock = nOptions/nBlocks;
trialMatrix = zeros(size(M));
for b = 1:nBlocks
    % select half class1, half class 2
    tmp = [M(idx1((b-1)*(nTrialsPerBlock/2)+1:b*(nTrialsPerBlock/2),1),:);M(idx2((b-1)*(nTrialsPerBlock/2)+1:b*(nTrialsPerBlock/2),1),:)];
    
    % permute order again
    P = randperm(size(tmp,1));
    
    % save as trial block
    trialMatrix((b-1)*nTrialsPerBlock+1:b*nTrialsPerBlock,:) = tmp(P,:);
   
end


save(sprintf('trialMatrix_%s',subjectID),'trialMatrix');
