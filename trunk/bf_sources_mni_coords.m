function res = bf_sources_mni_coords(BF, S)
% Generate beamforming grid
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% MWW
% $Id: bf_sources_mni_coords.m 4847 2012-08-16 17:29:23Z vladimir $

%--------------------------------------------------------------------------
if nargin == 0 
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'Pos coords';        
    pos.val = {};
     
    mni_coords = cfg_branch;
    mni_coords.tag = 'mni_coords';
    mni_coords.name = 'Mni Coords';
    mni_coords.val = {pos};    
    
    res=mni_coords;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% transform MNI coords in MNI space into space where we are doing the
% beamforming
M = inv(BF.data.transforms.toMNI);

grid.pos   = S.pos;

res = ft_transform_geometry(M, grid);