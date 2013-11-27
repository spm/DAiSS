function out = bf_data
% Prepares the data and initialises the beamforming pipeline
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the B.mat file containing the beamforming data will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG dataset';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the forward model is stored.'};
val.val = {1};

gradsource = cfg_menu;
gradsource.tag = 'gradsource';
gradsource.name = 'Where to get MEG sensors';
gradsource.help = {'Taking sensors from D.sensors makes it possible to',...
    'use the same head model and coregistration with multiple datasets.',...
    'This relies on the assumption that the sensors are in head coordinates',...
    'and the fiducals are at the same locations'};
gradsource.labels = {'D.inv', 'D.sensors'};
gradsource.values = {'inv', 'sensors'};
gradsource.val = {'inv'};

space = cfg_menu;
space.tag = 'space';
space.name = 'Coordinate system to work in';
space.help = {'Select the coordinate system for the forward model'};
space.labels = {'MNI-aligned', 'Head', 'Native'};
space.values = {'MNI-aligned', 'Head', 'Native'};
space.val = {'MNI-aligned'};

overwrite = cfg_menu;
overwrite.tag = 'overwrite';
overwrite.name = 'Overwrite BF.mat if exists';
overwrite.help = {'Choose whether to overwrite the existing BF.mat file'};
overwrite.labels = {'Yes', 'No'};
overwrite.values = {1, 0};
overwrite.val = {0};

out = cfg_exbranch;
out.tag = 'data';
out.name = 'Prepare data';
out.val = {dir, D, val, gradsource, space, overwrite};
out.help = {'Prepare the input for beamforming'};
out.prog = @bf_data_run;
out.vout = @bf_data_vout;
out.modality = {'EEG'};
end

function  out = bf_data_run(job)

outdir     = job.dir{1};
val        = job.val;
space      = job.space;
gradsource = job.gradsource;
D      = spm_eeg_load(job.D{1});

if ~isfield(D, 'inv')
    error('Please run head model specification.');
end

if numel(D.inv) < val
    error('Invalid inversion index');
end

cd(outdir);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(pwd,'BF.mat'),'file') && ~job.overwrite
    str = {'Current directory contains existing BF file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing BF file)',spm('time'));
        out = []; return
    end
end

BF        = [];

BF.data.D = D;
BF.data.mesh = D.inv{val}.mesh;

eegind = 0;
megind = 0;
for m = 1:numel(D.inv{val}.forward)
    if strncmp('EEG', D.inv{val}.forward(m).modality, 3)
        eegind = m;
    elseif strncmp('MEG', D.inv{val}.forward(m).modality, 3)
        megind = m;
    end
end

istemplate = D.inv{val}.mesh.template;

if megind > 0
    
    siunits  = isfield(D.inv{val}.forward(megind), 'siunits') &&...
        D.inv{val}.forward(megind).siunits;
    
    datareg  = D.inv{val}.datareg(megind);
    forward  = D.inv{val}.forward(megind);
    
    vol      = forward.vol;
    
    if isequal(gradsource, 'inv')
        if siunits
            sens     =  forward.sensors;
            toMNI    =  forward.toMNI;
            to_mm    = diag([1e3 1e3 1e3 1]);
        else
            sens     = datareg.sensors;
            toMNI    = datareg.toMNI;
            to_mm    = eye(4);
        end
    else
        sens     = D.sensors('MEG');
    end      
    
    if isfield(forward, 'mesh_correction')
        BF.data.MEG.mesh_correction = forward.mesh_correction;
    else
        BF.data.MEG.mesh_correction = [];
    end
      
    if siunits
        sens  = ft_convert_units(sens, 'm');
    end
    
    M          = to_mm\toMNI;
    [U, L, V]  = svd(M(1:3, 1:3));
    M(1:3,1:3) = U*V';    
    
    switch space
        case 'MNI-aligned'            
            BF.data.MEG.vol  = ft_transform_vol(M, vol);
            BF.data.MEG.sens = ft_transform_sens(M, sens);            
            
            BF.data.transforms.toMNI         = toMNI/M;
            BF.data.transforms.toMNI_aligned = to_mm;
            BF.data.transforms.toHead        = inv(M);
            BF.data.transforms.toNative      = D.inv{val}.mesh.Affine\BF.data.transforms.toMNI;
        case 'Head'
            BF.data.MEG.vol  = vol;
            BF.data.MEG.sens = sens;     
            
            BF.data.transforms.toMNI         = toMNI;
            BF.data.transforms.toMNI_aligned = to_mm*M;
            BF.data.transforms.toHead        = eye(4);
            BF.data.transforms.toNative      = D.inv{val}.mesh.Affine\BF.data.transforms.toMNI;
        case 'Native'
           error('Native coordinates option is deprecated for MEG.');
    end
end
            

if eegind > 0
    siunits  = isfield(D.inv{val}.forward(eegind), 'siunits') &&...
        D.inv{val}.forward(eegind).siunits;
        
    forward  = D.inv{val}.forward(eegind);    
    datareg  = D.inv{val}.datareg(eegind);
    
    vol      = forward.vol;
    
    if siunits
        sens     = forward.sensors;
        toMNI    = forward.toMNI;
        to_mm    = diag([1e3 1e3 1e3 1]);
    else
        sens     = datareg.sensors;
        toMNI    = datareg.toMNI;
        to_mm    = eye(4);
    end   
    
    BF.data.EEG.vol  = vol;        
    BF.data.EEG.sens = sens;             
                
    if isfield(forward, 'mesh_correction')
        BF.data.EEG.mesh_correction = forward.mesh_correction;
    else
        BF.data.EEG.mesh_correction = [];
    end
    
    if isfield(BF.data, 'transforms')  % With MEG
        if istemplate
            error('Combining EEG and MEG cannot be done with template head for now.');
        else
            if isa(vol, 'char')
                vol = ft_read_vol(vol);
            end
            
            fromNative = inv(BF.data.transforms.toNative);
            
            if siunits
                fromNative = to_mm\fromNative;
            end
            
            BF.data.EEG.vol  = ft_transform_vol(fromNative, vol);
            BF.data.EEG.sens = ft_transform_sens(fromNative, sens);
        end
    else                             % EEG only        
        M          = to_mm\toMNI;
        [U, L, V]  = svd(M(1:3, 1:3));
        M(1:3,1:3) = U*V';
        
        switch space
            case 'Native'
                BF.data.EEG.vol  = vol;
                BF.data.EEG.sens = sens;
                
                BF.data.transforms.toMNI         = toMNI;
                BF.data.transforms.toMNI_aligned = to_mm*M;
                BF.data.transforms.toHead        = eye(4); 
                BF.data.transforms.toNative      = to_mm; 
            case {'MNI-aligned'}
                if isa(vol, 'char')
                    vol = ft_read_vol(vol);
                end
                
                BF.data.EEG.vol  = ft_transform_vol(M, vol);
                BF.data.EEG.sens = ft_transform_sens(M, sens);
                
                BF.data.transforms.toMNI         = toMNI/M;
                BF.data.transforms.toMNI_aligned = to_mm;
                BF.data.transforms.toHead        = inv(M);
                BF.data.transforms.toNative      = to_mm/M;
            case {'Head'}              
               error('Head space is not defined for EEG data');
        end
    end
end

BF.data.space   = space;
BF.data.siunits = siunits;

bf_save(BF, 'overwrite');

out.BF{1} = fullfile(outdir, 'BF.mat');

end

function dep = bf_data_vout(job)
% Output is always in field "BF", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
