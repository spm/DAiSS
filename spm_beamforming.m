function spm_beamforming
% GUI gateway to Beamforming toolbox
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

pipelines = spm_select('List', fullfile(spm('dir'), 'toolbox', 'Beamforming'), '^bf_pipeline_.*\.m$');
pipelines = cellstr(pipelines);
pipelines = [cell(length(pipelines), 1), pipelines(:)];
for i = 1:size(pipelines, 1)
    [junk, pipelines{i, 2}] = fileparts(pipelines{i, 2});
    pipelines{i, 1} = strrep(pipelines{i, 2}, 'bf_pipeline_', '');
    pipelines{i, 1} = strrep(pipelines{i, 1}, '_', ' ');
end

str = sprintf('%s|', pipelines{:, 1});
str = str(1:(end-1));

fun = spm_input('Beamforming',1,'m', str, char(pipelines(:, 2)));
  
eval(fun);

spm_jobman('interactive', matlabbatch);