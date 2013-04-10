function res = bf_regularise_manual(BF, S)
% Manual specification of the regularisation parameter
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

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
    
    res      = cfg_branch;
    res.tag  = 'manual';
    res.name = 'User-specified regularisation';
    res.val  = {lambda};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

C = BF.features.(S.modality).C;

lambda = (S.lambda/100) * trace(C)/size(C,1);
C      = C + lambda * eye(size(C));
Cinv   = pinv_plus(C);
U      = eye(size(C));

features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;