function res = bf_features_regcov(BF, S)
% Simple covariance computation with regularization
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0 
    lambda = cfg_entry;
    lambda.tag = 'lambda';
    lambda.name = 'Regularisation';
    lambda.strtype = 'r';
    lambda.num = [1 1];
    lambda.val = {0};
    lambda.help = {'Select the regularisation (in %)'};
 
    regcov      = cfg_branch;
    regcov.tag  = 'regcov';
    regcov.name = 'Regularized covariance';
    regcov.val  = {lambda};
    
    res = regcov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

YY = 0;
ns = 0;

ntrials = length(S.trials);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

for i = 1:ntrials
    for j = 1:numel(S.samples)
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));
        
        % FIXME: I think, We need to remove the temporal mean according to the definition of the covariance matrix and not the spatial mean
        % We can do two things here: loading the whole dataset and remove the mean at once or removing the temporal mean of each trial
        % I would recommend the first one and consider the second on as an option, since the data has already been baseline corrected
        % Note that we can write an online version as well, if you dont one to load wole trials at once.    
        
        Y  = detrend(Y', 'constant');
        YY = YY+(Y'*Y);
        ns = ns + length(S.samples{j})-1;
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

C = YY/ns;

lambda = (S.lambda/100) * trace(C)/size(C,1);

C = C + lambda * eye(size(C));

res = C;