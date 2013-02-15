function res = bf_inverse_lcmv_use_pca_order(BF, S)
% Computes LCMV filters using spm_pca_order to constrain inverse of data
% cov matrix
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
%
% Based on the paper:
% MEG beamforming using Bayesian PCA for adaptive data covariance matrix regularization.
% Woolrich M, Hunt L, Groves A, Barnes G.
% Neuroimage. 2011 Aug 15;57(4)
%
% Mark Woolrich
% $Id: bf_inverse_lcmv_use_pca_order.m 4847 2012-08-16 17:29:23Z vladimir $

%--------------------------------------------------------------------------
if nargin == 0      
    
    pca_order = cfg_entry;
    pca_order.tag = 'pca_order';
    pca_order.name = 'PCA order';        
    pca_order.val = {};
   
    lcmv      = cfg_branch;
    lcmv.tag  = 'lcmv_use_pca_order';
    lcmv.name = 'LCMV use PCA order';
    lcmv.val  = {pca_order};      
    
    res = lcmv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = {'MEG'};

for m = 1:numel(modalities)
    if isfield(BF.features,modalities{m})
        
        reduce_rank=BF.sources.reduce_rank.(modalities{m});
        
        C = BF.features.(modalities{m}).C;
        
        [invCy, pca_order_used] = pinv_plus(C, S.pca_order); % rank maybe not be detected properly by just using pinv - MWW
        
        L = BF.sources.L.(modalities{m});
        
        W = cell(size(L));
        
        nvert = numel(W);
        
        spm('Pointer', 'Watch');drawnow;
        spm_progress_bar('Init', nvert, ['Computing ' modalities{m} ' filters']); drawnow;
        if nvert > 100, Ibar = floor(linspace(1, nvert,100));
        else Ibar = 1:nvert; end
        
        for i = 1:nvert
            if ~isnan(L{i})
                lf    = L{i};
                
                % Robert's code
                tmp=lf' * invCy *lf;                
                [u, ~] = svd(real(pinv_plus(tmp,reduce_rank,0)),'econ'); % this is faster,  - MWW
                                      
                eta = u(:,1);
                lf  = lf * eta;
                          
                % construct the spatial filter
                %W{i} = pinv(lf' * invCy * lf) * lf' * invCy;                
                W{i} = lf'*invCy/(lf' * invCy * lf); % this is faster - MWW
                
            else
                W{i} = NaN;
            end
            
             if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end
        end
        
        
        spm_progress_bar('Clear');
        inverse.W.(modalities{m}) = W;
    end
end

res = inverse;