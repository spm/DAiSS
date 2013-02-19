function res = bf_write_gifti(BF, S)
% Writes out beamformer results as GIfTI meshes
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    normalise         = cfg_menu;
    normalise.tag     = 'normalise';
    normalise.name    = 'Global normalisation';
    normalise.help    = {'Normalise image values by the mean'};
    normalise.labels  = {
        'No'
        'Each image separately'
        'Across images'
        }';
    normalise.values  = {
        'no'
        'separate'
        'all'
        }';
    normalise.val = {'separate'};
    
    space         = cfg_menu;
    space.tag     = 'space';
    space.name    = 'Image space';
    space.help    = {'Specify image space'};
    space.labels  = {
        'MNI'
        'Native'
        'MNI-aligned'
        }';
    space.values  = {
        'mni'
        'native'
        'aligned'
        }';
    space.val = {'mni'};
    
    res      = cfg_branch;
    res.tag  = 'gifti';
    res.name = 'GIfTI';
    res.val  = {normalise, space};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

if ~isfield(BF.sources, 'mesh')
    error('Only mesh source space is supported');
end

scale = ones(1, numel(BF.output.image));
switch S.normalise
    case  'separate'
        for i = 1:numel(BF.output.image)
            val = BF.output.image(i).val;
            scale(i) = 1./mean(val(~isnan(val)));
        end
    case  'all'
        val = spm_vec({BF.output.image(:).val});
        scale = scale./mean(val(~isnan(val)));
end

switch S.space
    case 'mni'
        source = BF.sources.mesh.canonical;
    case 'aligned'
        source = BF.sources.mesh.individual;
        source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toMNI_aligned, source.vert);
    case 'native'
        source = BF.sources.mesh.individual;
        source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toNative, source.vert);
end

source = export(gifti(source), 'patch');

nimages = numel(BF.output.image);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nimages , 'Writing out images'); drawnow;
if nimages  > 100, Ibar = floor(linspace(1, nimages ,100));
else Ibar = 1:nimages; end

for i = 1:nimages
    fname = fullfile(pwd, [BF.output.image(i).label '.gii']);
    
    source.cdata = scale(i)*BF.output.image(i).val;
    
    save(gifti(source), fname);
            
    res.files{i, 1} = fname;
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

