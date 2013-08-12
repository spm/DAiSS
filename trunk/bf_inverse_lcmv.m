function res = bf_inverse_lcmv(BF, S)
% Computes LCMV filters
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    keeplf        = cfg_menu;
    keeplf.tag    = 'keeplf';
    keeplf.name   = 'Keep oriented leadfields';
    keeplf.labels = {'yes', 'no'};
    keeplf.values = {true, false};
    keeplf.val    = {false};
    
    lcmv      = cfg_branch;
    lcmv.tag  = 'lcmv';
    lcmv.name = 'LCMV';
    lcmv.val  = {keeplf};
    res = lcmv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


invCy =  BF.features.(S.modality).Cinv;
U     = BF.features.(S.modality).U;

reduce_rank = BF.sources.reduce_rank.(S.modality(1:3));

L = S.L;
W = cell(size(L));

nvert = numel(W);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Computing ' S.modality ' filters']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    if ~isnan(L{i})
        lf    = U'*L{i};
        
        % Robert's code
        [u, dum] = svd(real(pinv_plus(lf' * invCy *lf, reduce_rank, 0)),'econ');
        eta = u(:,1);
        lf  = lf * eta;
        
        % construct the spatial filter
        W{i} = lf'*invCy/(lf' * invCy * lf);
        
        if S.keeplf
            L{i} = lf;
        end
    else
        W{i} = NaN;
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

res.W = W;
res.L = L;