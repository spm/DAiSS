function res = bf_inverse_lcmv(BF, S)
% Computes LCMV filters
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0      
    lcmv      = cfg_branch;
    lcmv.tag  = 'lcmv';
    lcmv.name = 'LCMV';
    lcmv.val  = {};
    
    res = lcmv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = {'MEG', 'EEG'};

for m = 1:numel(modalities)
    if isfield(BF.features.C, modalities{m})
        C = BF.features.C.(modalities{m});
        
        invCy =  pinv(C); % assuming already regularized
        
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
                
                % FIXME: The original 'pinv' in matlab does not necessarily work here. Roberts uses a pinv with 10 times bigger tol and Mark has his own pinv. 
                % I would recommend to write a new pinv that can be used throughout this software.
                % It is also clear that if you regularise eoungh the covariance matrix there is no need to a new pinv
                
                % Robert's code
                [u, s, v] = svd(real(pinv(lf' * invCy *lf)));
                eta = u(:,1);
                lf  = lf * eta;
                          
                % construct the spatial filter
                W{i} = pinv(lf' * invCy * lf) * lf' * invCy;
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