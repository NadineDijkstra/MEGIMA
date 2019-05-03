function source = sourceAnalysisTwoConditions(cfg0, data, D)

%% Define useful variables
iInside = find(cfg0.lf.inside);

numF = size(data.trial, 2);                    % Number of sensors
numT = size(data.trial, 3);                             % Number of time points
numN = size(data.trial, 1);                             % Number of trials
numS = length(cfg0.lf.leadfield);              % Number of sources
numW = length(iInside);                        % Number of filters (i.e. sources *inside* the brain)
numOri = size(cfg0.lf.leadfield{iInside(1)}, 2);      % Number of orientations in leadfield

if rank(cfg0.lf.leadfield{iInside(1)}) ~= numOri
    error(sprintf('Leadfields are rank-deficient.\nDid you specify reducerank = ''yes'', but *not* projectback = ''no'' for ft_prepare_leadfield?'));
end

if strcmp(cfg0.allowNan, 'yes')
    covFun = @nancov;
    meanFun = @nanmean;
else
    covFun = @cov;
    meanFun = @mean;
end    

%% Calculate covariance
if strcmp(cfg0.cov.method, 'fixed')
    cfg = [];
    cfg.latency = cfg0.cov.window;    
    covData = ft_selectdata(cfg, data);

    S = covFun(reshape(permute(covData.trial, [2, 1, 3]), numF, [])');
    
    % Regularize
    S = (1-cfg0.lambda)*S + cfg0.lambda*eye(numF)*trace(S)/numF;     
elseif strcmp(cfg0.cov.method, 'moving')
    % Calculate covariance per time point
    S_raw = nan(numF, numF, numT);

    tStart = tic;
    for it = 1:numT
        S_raw(:, :, it) = covFun(data.trial(:, :, it));
        
        % Feedback
        if (toc(tStart) > cfg0.feedback)
            fprintf('%s - finished calculating covariance for time point %g/%g\n', mfilename, it, numT);
            tStart = tic;
        end
    end

    % Create smoothing kernel
    if strcmp(cfg0.cov.window, 'gaussian')       
        numK = floor(cfg0.cov.windowWidth * (data.fsample/1000) * 3)*2 + 1;
        kernel = exp(-((1:numK)' - 0.5 - numK/2).^2 / (cfg0.cov.windowWidth^2));
        kernel = kernel/sum(kernel);
    else
        numK = floor(cfg0.cov.windowWidth * (data.fsample/1000));
        if rem(numK, 2) == 0            % Make odd number for symmetric centering
            numK = numK + 1;
        end
        
        kernel = ones(numK, 1)/numK;
    end

    % Smooth and regularize covariance
    S_raw = reshape(S_raw, [numF*numF, numT]);
    S_raw = [zeros(numF*numF, (numK-1)/2), S_raw, zeros(numF*numF, (numK-1)/2)];

    S = nan(numF, numF, numT);
    
    tStart = tic;
    for it = 1:numT    
        S(:, :, it) = reshape(S_raw(:, (1:numK) - 1 + it)*kernel, [numF, numF]);

        % Regularize
        S(:, :, it) = ...
            (1-cfg0.lambda)*S(:, :, it) ...
            + cfg0.lambda*eye(numF)*trace(S(:, :, it))/numF; 

        % Feedback
        if (toc(tStart) > cfg0.feedback)
            fprintf('%s - finished smoothing and regularizing covariance for time point %g/%g\n', mfilename, it, numT);
            tStart = tic;
        end
    end

    clear S_raw;
end

%% Do source permutation analysis
if strcmp(cfg0.cov.method, 'fixed')   
    % Calculate spatial filters
    W = nan(numF, numOri, numW);
    Sinv = S\eye(numF);
    
    for iw = 1:numW
        L = cfg0.lf.leadfield{iInside(iw)};
                       
        W(:, :, iw) = Sinv*L;
        W(:, :, iw) = W(:, :, iw)/(L'*W(:, :, iw));
    end
    
    W = reshape(W, [numF, numOri*numW]);
    
    % Iterate over permutations
    mPerm = zeros(numW, numT);
    
    tStart = tic;
    for iPerm = 1:cfg0.numPerm
        % Flip trial labels
        Dperm = D(randperm(numN));
        
        % Average
        m = squeeze(meanFun(data.trial(Dperm==1, :, :), 1)) - squeeze(meanFun(data.trial(Dperm==0, :, :), 1));
        if numT == 1; m = m'; end
        
        % Source project and reshape        
        m = reshape(W'*m, [numOri, numW, numT]);
        
        % Calculate power of vector
        m = squeeze(sum(m.^2, 1));
        if numT == 1; m = m'; end
        
        % Store results from this permutation
        mPerm = mPerm + m/cfg0.numPerm;

        % Feedback
        if (toc(tStart) > cfg0.feedback)
            fprintf('%s - finished doing permutation %g/%g\n', mfilename, iPerm, cfg0.numPerm);
            tStart = tic;
        end        
    end

    % Process the real data
    m = squeeze(meanFun(data.trial(D==1, :, :), 1)) - squeeze(meanFun(data.trial(D==0, :, :), 1));
    if numT == 1; m = m'; end    
    m = reshape(W'*m, [numOri, numW, numT]);   
    m = squeeze(sum(m.^2, 1));
    if numT == 1; m = m'; end   
    
    % Subtract permutation mean
    y = nan(numS, numT);
    y(iInside, :) = m - mPerm;  
    
elseif strcmp(cfg0.cov.method, 'moving')    
    % Do source permutation analysis per time point
    tStart = tic;
    for it = 1:numT
        % Calculate spatial filters
        W = nan(numF, numOri, numW);
        Sinv = S(:, :, it)\eye(numF);

        for iw = 1:numW
            L = cfg0.lf.leadfield{iInside(iw)};

            W(:, :, iw) = Sinv*L;
            W(:, :, iw) = W(:, :, iw)/(L'*W(:, :, iw));
        end

        W = reshape(W, [numF, numOri*numW]);

        % Permutation #1 is true data
        % Calculate permutations at sensor level
        yPerm = nan(numF, cfg0.numPerm+1);

        for iPerm = 1:(cfg0.numPerm+1)        
            if (iPerm == 1)
                yPerm(:, iPerm) = meanFun(data.trial(:, :, it), 1)';
            else
                Dperm = D(randperm(numN));

                yPerm(:, iPerm) = squeeze(meanFun(data.trial(Dperm==1, :, it), 1)) - squeeze(meanFun(data.trial(Dperm==0, :, it), 1));;
            end
        end

        % Source projection
        yPerm = reshape(W'*yPerm, [numOri, numW, cfg0.numPerm+1]);

        % Calculate power of vector
        yPerm = squeeze(sum(yPerm.^2, 1));

        % Calculate permutation mean
        mPerm = mean(yPerm(:, 2:(cfg0.numPerm+1)), 2);

        % Subtract permutation mean
        y(iInside, it) = yPerm(:, 1) - mPerm;

        % Feedback
        if (toc(tStart) > cfg0.feedback)
            fprintf('%s - finished doing source analysis for time point %g/%g\n', mfilename, it, numT);
            tStart = tic;
        end    
    end
end

% Calculate scaling factor
alpha = nan(numS, 1);

for iw = 1:numW
    L = cfg0.lf.leadfield{iInside(iw)};
    alpha(iInside(iw)) = 1/trace(eye(numOri)/(L'*L));
end

if strcmp(cfg0.locationRescale, 'yes')
    y = y .* repmat(alpha, [1, numT]);
end

% Take square root
if strcmp(cfg0.returnSqrt, 'yes')
    y(y < 0) = 0;
    y = sqrt(y);
    
    alpha = sqrt(alpha);
end    

% Return as FieldTrip-struct
inside = zeros(numS, 1);
inside(iInside) = 1;

source = [];
source.time = data.time;
source.inside = inside;
source.pos = cfg0.lf.pos;
source.avg.pow = y;
source.mPerm = nan(numS, numT);
source.mPerm(iInside, :) = mPerm;           % also return the permutation mean
source.alpha = alpha;

end






