function res = bf_inverse_lcmvmontage(BF, S)
% Computes LCMV filters
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    usemontage = cfg_menu;
    usemontage.tag = 'usemontage';
    usemontage.name = 'You can use a reduced channel montage if you like';
    usemontage.help = {'reduce effective number of channels based on previous steps'};
    usemontage.labels = {'Re-montage', 'Use original'};
    usemontage.values = {'Re-montage', 'Use original'};
    usemontage.val = {'Re-montage'};
    
    lcmvmontage      = cfg_branch;
    lcmvmontage.tag  = 'lcmvmontage';
    lcmvmontage.name = 'LCMV with montage option';
    lcmvmontage.val  = {usemontage};
    
    
    res = lcmvmontage;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = {'MEG', 'EEG'};

for m = 1:numel(modalities)
    
    
    if isfield(BF.features, modalities{m})
        if ~isfield(BF.features.(modalities{m}),'montage')
            S.usemontage='Use original';
        end;
        switch S.usemontage,
            case 'Re-montage'
                C = BF.features.(modalities{m}).montage.C;
                
                invCy =  BF.features.(modalities{m}).montage.Cinv;
                U=BF.features.(modalities{m}).montage.U; %% new montage
                
            case 'Use original'
                
                C = BF.features.(modalities{m}).C;
                
                invCy =  BF.features.(modalities{m}).Cinv;;
                
                U=eye(size(C,1));
                disp('Not re-montaging');
        end; % switch
    
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
    
    end; %% if modality exists
end; % for modalities



res = inverse;