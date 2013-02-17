function res = bf_write_spmeeg(BF, S)
% Writes out beamformer results as M/EEG dataset
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    mode         = cfg_menu;
    mode.tag     = 'mode';
    mode.name    = 'Writing mode';
    mode.help    = {'How to generate the output'};
    mode.labels  = {
        'Write new'
        'Online montage original data'
        'Copy + online montage'
        }';
    mode.values  = {
        'write'
        'online'
        'onlinecopy'
        }';
    mode.val = {'write'};
    
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'What modality to output'};
    modality.labels  = {
        'MEG'
        'EEG'
        }';
    modality.values  = {
        'MEG'
        'EEG'
        }';
    modality.val = {'MEG'};
    
    none = cfg_const;
    none.tag = 'none';
    none.name = 'None';
    none.val  = {0};
    
    addchannels      = cfg_choice;
    addchannels.tag  = 'addchannels';
    addchannels.name = 'Extra channels to add';
    addchannels.values  = {none, spm_cfg_eeg_channel_selector};
    addchannels.val  = {none};
    
    prefix         = cfg_entry;
    prefix.tag     = 'prefix';
    prefix.name    = 'Filename Prefix';
    prefix.help    = {'Specify the string to be prepended to the output (if relevant).'};
    prefix.strtype = 's';
    prefix.num     = [1 Inf];
    prefix.val     = {'B'};
    
    spmeeg      = cfg_branch;
    spmeeg.tag  = 'spmeeg';
    spmeeg.name = 'SPM M/EEG dataset';
    spmeeg.val  = {mode, modality, addchannels, prefix};
    
    res = spmeeg;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

if isfield(S.addchannels, 'channels')
    addchannels = D.chanlabels(D.selectchannels(spm_cfg_eeg_channel_selector(S.addchannels.channels)));
else
    addchannels = {};
end

usemontage = 0; % Might want to support writing out data directly from BF.output in the future

if isfield(BF.output, 'montage')
    
    usemontage = 1;
    
    montage = BF.output.montage.(S.modality);
    
    [sel1, sel2] = spm_match_str(montage.labelorg, addchannels);
    for i = 1:length(sel1)
        montage.labelnew(end+1) = addchannels(sel2(i));
        montage.tra(end+1, sel1(i)) = 1;
    end
    addchannels(sel2) = [];
    
    if ~isempty(addchannels)
        montage.labelorg = [montage.labelorg(:); addchannels(:)];
        montage.labelnew = [montage.labelnew(:); addchannels(:)];
        
        na = numel(addchannels);
        montage.tra((end+1):(end+na), (end+1):(end+na)) = eye(na);
    end
end

if usemontage
    S1 = [];
    S1.montage = montage;
    S1.prefix  = S.prefix; % ignored for online
        
    switch S.mode
        case 'write'                                              
            S1.mode    = 'write';
        case 'online'
            S1.mode = 'switch';
        case 'onlinecopy' 
            S1.mode    = 'switch';            
            D = copy(D, [S.prefix D.fname]);
    end        
    
    S1.D = D;
    D = spm_eeg_montage(S1);
end

res = fullfile(D);
