function res = bf_group_batch(BF, S)
% Run a DAiSS batch on a group of subjects
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0 
    batchfile = cfg_files;
    batchfile.tag = 'batchfile';
    batchfile.name = 'Batch .mat or .m file';
    batchfile.filter = '(.*.mat$)|(.*.m$)';
    batchfile.num = [1 Inf];
    batchfile.help = {'Select batch specification file.'};
    
    batch = cfg_branch;
    batch.tag = 'batch';
    batch.name = 'Batch';
    batch.val = {batchfile};
    
    res = batch;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

batchfile = char(S.batchfile);

if isequal(spm_file(batchfile, 'ext'), 'm')
    cdir = pwd;
    if ~isempty(spm_file(batchfile, 'path'))
        cd(spm_file(batchfile, 'path'));
    end
    vars = who;
    eval(spm_file(batchfile, 'basename'));
    cd(cdir);
    name = setdiff(who, [vars; {'vars'}]);
    if numel(name)~= 1
        error('Invalid batch specification');
    end
    matlabbatch = eval(char(name));
    if ~isa(matlabbatch, 'cell')
        error('Invalid batch specification');
    end
elseif isequal(spm_file(batchfile, 'ext'), 'mat')
    tmp  = load(batchfile);
    name = fieldnames(tmp);
    if numel(name)~= 1
        error('Invalid batch specification');
    end
    matlabbatch = tmp.(name);
    if ~isa(matlabbatch, 'cell')
        error('Invalid batch specification');
    end
else
    error('Invalid batch specification');
end

try 
    matlabbatch{1}.spm.tools.beamforming;
    [~, matlabbatch] = spm_jobman('harvest', matlabbatch);
catch
    error('Invalid batch specification. DAiSS batch expected.');
end

res = cell(1, numel(BF));
for i = 1:numel(BF)
    if isfield(matlabbatch{1}.spm.tools.beamforming, 'data');
        D = spm_eeg_load(BF{i});
        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(D)};
        res = mkdir(D.path, [S.prefix 'BF']);
        matlabbatch{1}.spm.tools.beamforming.data.dir = {fullfile(D.path, [S.prefix 'BF'])};
        
    else
        module = fieldnames(matlabbatch{1}.spm.tools.beamforming);
        matlabbatch{1}.spm.tools.beamforming(module).BF = BF(i);        
    end
    out = spm_jobman('run', matlabbatch);
    res(i) = out(end);
end
