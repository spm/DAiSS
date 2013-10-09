function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

if nargin == 1 && exist(fullfile(pwd, 'BF.mat'), 'file')
    save('BF.mat', '-struct', 'BF', '-append', '-v7.3');
else
    save('BF.mat', '-struct', 'BF', '-v7.3');
end