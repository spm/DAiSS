function res = bf_features_csd(BF, S)
% Compute cross-spectral density matrix for DICS
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0   

    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency band of interest';
    foi.strtype = 'r';
    foi.num = [1 2];
    foi.help = {'Frequency window within which to compute CSD over (Hz)'};
    
    taper = cfg_menu;
    taper.tag = 'taper';
    taper.name = 'Taper';
    taper.help = {'Save taper as well as power'};
    taper.labels = {'Hanning', 'Rectangular', 'DPSS', 'Sine'};
    taper.values = {'hanning', 'rectwin', 'dpss', 'sine'};
    taper.val = {'dpss'};
    
    lambda = cfg_entry;
    lambda.tag = 'lambda';
    lambda.name = 'Regularisation';
    lambda.strtype = 'r';
    lambda.num = [1 1];
    lambda.val = {0};
    lambda.help = {'Select the regularisation (in %)'};
            
    keepreal = cfg_menu;
    keepreal.tag = 'keepreal';
    keepreal.name = 'Keep real';
    keepreal.labels = {'Yes', 'No'};
    keepreal.val = {0};
    keepreal.values = {1, 0};
    keepreal.help = {'Keep only the real part of the CSD'};
    
    csd      = cfg_branch;
    csd.tag  = 'csd';
    csd.name = 'Cross-spectral density';
    csd.val  = {foi, taper, lambda, keepreal};
    
    res = csd;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

ntrials    = length(S.trials);
centerfreq = mean(S.foi);
tapsmofrq  = 0.5*(abs(diff(S.foi)));

Cf = 0;

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing CSD'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

for i = 1:ntrials
    for j = 1:numel(S.samples)
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));   

        
        [fourier, ntap] = ft_specest_mtmfft(Y, D.time(S.samples{j}), 'freqoi', centerfreq, ...
            'tapsmofrq', tapsmofrq, 'taper', S.taper, 'verbose', 0);
        
        
        dat  = transpose(fourier);
        Cf = Cf + (dat * ctranspose(dat)) ./ ntap;
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

Cf = Cf/ntrials;

lambda = (S.lambda/100) * trace(Cf)/size(Cf,1);

if S.keepreal
    % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
    % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
    invCf = pinv(real(Cf) + lambda * eye(size(Cf)));
else
    invCf = pinv(Cf + lambda * eye(size(Cf)));
end

features=[];

features.C    = Cf;
features.Cinv = invCf;

res = features;