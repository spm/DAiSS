function voi = bf_sources_voi(BF, S)
% Generate a set of VOIs specified in MNI coordinates
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% $Id$

%--------------------------------------------------------------------------
if nargin == 0 
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the VOI'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates'; 
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the VOI in MNI coordinates'};
    pos.val = {};  
            
    ori = cfg_entry;
    ori.tag = 'ori';
    ori.name = 'Orientation';
    ori.strtype = 'r';
    ori.num = [1 3];
    ori.help = {'Source orientatons (only for single points, leave zeros for unoriented)'};
    ori.val = {[0 0 0]};
    
    voidef = cfg_branch;
    voidef.tag = 'voidef';
    voidef.name = 'VOI';
    voidef.val = {label, pos, ori};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'VOIs';
    vois.num  = [1 Inf];
    vois.values = {voidef};
    vois.val = {voidef};
    
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the VOIs (leave 0 for single point)'};
    
    resolution = cfg_entry;
    resolution.tag = 'resolution';
    resolution.name = 'Resolution';
    resolution.strtype = 'r';
    resolution.num = [1 1];
    resolution.val = {5};
    resolution.help = {'Resolution for placing grid points in each VOI (in mm)'};
    
    voi = cfg_branch;
    voi.tag = 'voi';
    voi.name = 'VOIs in MNI space';
    voi.val = {vois, radius, resolution};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% transform MNI coords in MNI space into space where we are doing the
% beamforming
M = inv(BF.data.transforms.toMNI);

nvoi = numel(S.voidef);
voi = [];
voi.label = {S.voidef(:).label}';
if S.radius > 0
    vec = -S.radius:S.resolution:S.radius;
    [X, Y, Z]  = ndgrid(vec, vec, vec);
    sphere   = [X(:) Y(:) Z(:)];
    sphere(sqrt(X(:).^2 + Y(:).^2 + Z(:).^2) > S.radius, :) = [];
    npnt = size(sphere, 1);
else
    sphere = 0;
    npnt = 1;
    ori = cat(1, S.voidef(:).ori);
    if any(any(ori))
        voi.ori = ori;
    end
end

voi.pos = [];
voi.pos2voi = [];
for s = 1:nvoi
    voi.pos = [voi.pos; sphere+repmat(S.voidef(s).pos, npnt, 1)];
    voi.pos2voi = [voi.pos2voi s*ones(1, npnt)];
end

voi = ft_transform_geometry(M, voi);