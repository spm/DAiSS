function res = bf_inverse_lcmv(BF, S)
% Computes LCMV filters
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    lcmv      = cfg_const;
    lcmv.tag  = 'lcmv';
    lcmv.name = 'LCMV';
    lcmv.val  = {0};    
    res = lcmv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = {'MEG', 'EEG'};

for m = 1:numel(modalities)        
    if isfield(BF.features, modalities{m})
        invCy =  BF.features.(modalities{m}).Cinv;
        U     = BF.features.(modalities{m}).U;
        
        reduce_rank = BF.sources.reduce_rank.(modalities{m});
        
        L = BF.sources.L.(modalities{m});
        W = cell(size(L));
        
        nvert = numel(W);
        
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modalities{m} ' filters']); drawnow;
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
            else
                W{i} = NaN;
            end
            
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
        
        
        spm_progress_bar('Clear');
        inverse.W.(modalities{m}) = W;
        
    end; %% if modality exists
end; % for modalities



res = inverse;