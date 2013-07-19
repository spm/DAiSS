function out = bf_inverse
% Computes inverse projectors
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

% dir Directory
% ---------------------------------------------------------------------
BF = cfg_files;
BF.tag = 'BF';
BF.name = 'BF.mat file';
BF.filter = '^BF.mat$';
BF.num = [1 1];
BF.help = {'Select BF.mat file.'};

%--------------------------------------------------------------------------
% method
%--------------------------------------------------------------------------
plugin      = cfg_choice;
plugin.tag  = 'plugin';
plugin.name = 'Inverse method';

inverse_funs = spm_select('List', fileparts(mfilename('fullpath')), '^bf_inverse_.*\.m$');
inverse_funs = cellstr(inverse_funs );
for i = 1:numel(inverse_funs)
    plugin.values{i} = feval(spm_file(inverse_funs{i},'basename'));
end


out = cfg_exbranch;
out.tag = 'inverse';
out.name = 'Inverse solution';
out.val = {BF, plugin};
out.help = {'Compute inverse projectors'};
out.prog = @bf_inverse_run;
out.vout = @bf_inverse_vout;
out.modality = {'EEG'};
end

function  out = bf_inverse_run(job)

outdir = spm_file(job.BF{1}, 'fpath');

cd(outdir);

BF = bf_load('BF.mat', {'data', 'sources', 'features'});
plugin_name = cell2mat(fieldnames(job.plugin));
S = job.plugin.(plugin_name);

if ~isa(S, 'struct')
    S = [];
end

D = BF.data.D;

modalities = intersect(fieldnames(BF.features), {'MEG', 'MEGPLANAR', 'EEG'});

for m = 1:numel(modalities)
    S(1).modality = modalities{m};
    
    S.L = cell(size(BF.sources.L.(modalities{m}(1:3))));
    channels = cell(length(BF.features.(modalities{m}).chanind), 1);
    
    [sel1, sel2] = spm_match_str(D.chanlabels(BF.features.(modalities{m}).chanind),...
        BF.sources.channels.(modalities{m}(1:3)));
    
    nlfcolumns = max(cellfun('size', BF.sources.L.(modalities{m}(1:3)), 2));
    
    channels(sel1) = BF.sources.channels.(modalities{m}(1:3))(sel2);
    
    % This is for MEG-EEG fusion
    if isequal(modalities{m}, 'MEG') && length(sel1)<length(BF.features.('MEG').chanind)
        fuse = 1;
        [sel3, sel4] = spm_match_str(D.chanlabels(BF.features.('MEG').chanind),...
            BF.sources.channels.('EEG'));
        
        channels(sel3) = BF.sources.channels.EEG(sel4);
    else
        fuse = 0;
    end
    
    spm_progress_bar('Init', numel(S.L), ['Preparing ' modalities{m} ' leadfields']); drawnow;
    if numel(S.L) > 100, Ibar = floor(linspace(1, numel(S.L),100));
    else Ibar = 1:numel(S.L); end
    
    for i = 1:numel(S.L)
        if ~isnan(BF.sources.L.(modalities{m}(1:3)){i})
            S.L{i} = nan(length(BF.features.(modalities{m}).chanind), nlfcolumns);
            S.L{i}(sel1, :) = BF.sources.L.(modalities{m}(1:3)){i}(sel2, :);
            
            if fuse
                S.L{i}(sel3, :) = BF.sources.L.('EEG'){i}(sel4, :);
            end
        else
            S.L{i} = nan;
        end
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    spm_progress_bar('Clear');
    
    BF.inverse.(modalities{m}) = feval(['bf_inverse_' plugin_name], BF, S);
    BF.inverse.(modalities{m}).channels = channels;
end

bf_save(BF);

out.BF{1} = fullfile(outdir, 'BF.mat');
end


function dep = bf_inverse_vout(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end
