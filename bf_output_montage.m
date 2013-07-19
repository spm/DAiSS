function mont = bf_output_montage(BF, S)
% Generates a montage for source extraction
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0    
    woi      = cfg_entry;
    woi.tag  = 'woi';
    woi.name = 'Time window of interest';
    woi.strtype = 'r';
    woi.num = [Inf 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time window to optimise PCA (only used for VOI)'};
    
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd'};
    method.val = {'max'};
    method.values = {'max', 'svd'};
    method.help = {'How to summarise sources in the ROI'};
    
    mont = cfg_branch;
    mont.tag = 'montage';
    mont.name = 'Source montage';
    mont.val  = {woi, method};    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGPLANAR'});

for m  = 1:numel(modalities)    
    U        = BF.features.(modalities{m}).U; 
        
    montage          = [];
    montage.labelorg = BF.inverse.(modalities{m}).channels;
    montage.labelorg = montage.labelorg(:);
    montage.tra      = [];
    if isfield(BF.sources, 'voi')
        montage.labelnew = BF.sources.voi.label;
        for v = 1:numel(montage.labelnew)
            ind = find(BF.sources.voi.pos2voi == v);
            W   = cat(1, BF.inverse.(modalities{m}).W{ind});
            
            switch S.method
                case 'max'
                    
                    Wc          = W* BF.features.(modalities{m}).C*W';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    montage.tra = [montage.tra; W(mi, :)*U'];
                    
                case 'svd'
                    %% just take top pca component for now
                    Wc          = W* BF.features.(modalities{m}).C*W'; % bf estimated source covariance matrix
                    [V,dum,dum]=svd(Wc);
                    montage.tra=[montage.tra;V(:,1)'*W*U'];
                    
            end
        end
    else
        mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI, BF.sources.pos);
        for i = 1:size(mnipos, 1)
            w = BF.inverse.(modalities{m}).W{i};
            if ~isnan(w)
                montage.labelnew{i} = sprintf('%d_%d_%d', round(mnipos(i, :)));
                montage.tra = [montage.tra; w*U'];
            end
        end
    end
    
    mont.(modalities{m}) = montage;
end
