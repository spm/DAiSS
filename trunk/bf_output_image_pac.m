function res = bf_output_image_pac(BF, S)
% Computes phase-amplitude coupling
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Bernadette van Wijk, Vladimir Litvak
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    all = cfg_const;
    all.tag = 'all';
    all.name = 'All';
    all.val  = {1};
    
    condlabel = cfg_entry;
    condlabel.tag = 'condlabel';
    condlabel.name = 'Condition label';
    condlabel.strtype = 's';
    condlabel.val = {''};
    
    conditions = cfg_repeat;
    conditions.tag = 'conditions';
    conditions.name = 'Conditions';
    conditions.help = {'Specify the labels of the conditions to be included in the inversion'};
    conditions.num  = [1 Inf];
    conditions.values  = {condlabel};
    conditions.val = {condlabel};
    
    whatconditions = cfg_choice;
    whatconditions.tag = 'whatconditions';
    whatconditions.name = 'What conditions to include?';
    whatconditions.values = {all, conditions};
    whatconditions.val = {all};
    
    sametrials = cfg_menu;
    sametrials.tag = 'sametrials';
    sametrials.name = 'Trials same as for filters';
    sametrials.labels = {'yes', 'no'};
    sametrials.values = {true, false};
    sametrials.val = {false};
    sametrials.help = {'Take the same trials as used for filter computation',...
        'This is useful for bootstrap.'};
    
    woi = cfg_entry;
    woi.tag = 'woi';
    woi.name = 'Time window of interest';
    woi.strtype = 'r';
    woi.num = [1 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time windows (in ms)'};       
        
    phasefreq         = cfg_entry;
    phasefreq.tag     = 'phasefreq';
    phasefreq.name    = 'Phase frequencies';
    phasefreq.strtype = 'r';
    phasefreq.num     = [1 Inf];
    phasefreq.val     = {5:3:30};
    phasefreq.help    = {'Frequencies to compute phase for (as a vector)'};    
    
    
    phaseres         = cfg_entry;
    phaseres.tag     = 'phaseres';
    phaseres.name    = 'Phase resolution';
    phaseres.strtype = 'r';
    phaseres.num     = [1 1];
    phaseres.val     = {2};
    phaseres.help    = {'Frequency resolution for phase computation'};
    
    ampfreq         = cfg_entry;
    ampfreq.tag     = 'ampfreq';
    ampfreq.name    = 'Amplitude frequencies';
    ampfreq.strtype = 'r';
    ampfreq.num     = [1 Inf];
    ampfreq.val     = {30:5:100};
    ampfreq.help    = {'Frequencies to compute amplitude for (as a vector)'};
    
    ampres         = cfg_entry;
    ampres.tag     = 'ampres';
    ampres.name    = 'Amplitude resolution';
    ampres.strtype = 'r';
    ampres.num     = [1 1];
    ampres.val     = {15};
    ampres.help    = {'Frequency resolution for amplitude computation'};
    
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'Specify modality'};
    modality.labels  = {
        'MEG'
        'EEG'
        }';
    
    modality.values  = {
        'MEG'
        'EEG'
        }';
    modality.val = {'MEG'};
    
    image_pac      = cfg_branch;
    image_pac.tag  = 'image_pac';
    image_pac.name = 'PAC image';
    image_pac.val  = {whatconditions, sametrials, woi, phasefreq, phaseres, ampfreq, ampres, modality};
    
    res = image_pac;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

S.woi = 1e-3*S.woi; % ms -> s

samples =  D.indsample(S.woi(1)):D.indsample(S.woi(2));
nsamples = length(samples);

if isfield(S.whatconditions, 'all')
    S.whatconditions.condlabel = D.condlist;
end

for i = 1:numel(S.whatconditions.condlabel)
    if S.sametrials
        trials{i} = BF.features.trials(strmatch(S.whatconditions.condlabel{i},...
            D.conditions(BF.features.trials)));
    else
        trials{i} = D.indtrial(S.whatconditions.condlabel{i}, 'GOOD');
    end
    
    if isempty(trials{i})
        error('No trials matched the selection.');
    end
    
end

if isempty(trials)
    error('No trials matched the selection, check the specified condition labels');
end


channels = BF.features.(S.modality).chanind;
U        = BF.features.(S.modality).U;
nchan    = size(U, 2);

alltrials = spm_vec(trials);
ntrials   = length(alltrials);

nphase  = length(S.phasefreq);
namp    = length(S.ampfreq);

W = BF.inverse.W.(S.modality);
nvert = numel(W);

spm_progress_bar('Init', nvert, ...
    sprintf('Scanning grid points image')); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

pac   = nan(nphase, namp, nvert);
phase = nan(nphase, nsamples*ntrials);
amp   = nan(namp, nsamples*ntrials);

for i = 1:nvert
    if ~isnan(W{i})
        w    = W{i};
        
        Y  = w*U'*reshape(D(channels, samples, alltrials), nchan, []);
        
        Y  = reshape(Y, nsamples, ntrials)';
        
        for f = 1:nphase
             fY = ft_preproc_bandpassfilter(Y, D.fsample, S.phasefreq(f) + S.phaseres*[-1 1], 2);
             phase(f, :) = pi + reshape(angle(hilbert(fY')), 1, []);
        end
        
        for f = 1:namp
            fY = ft_preproc_bandpassfilter(Y, D.fsample, S.ampfreq(f) + S.ampres*[-1 1], 2);
            amp(f, :) =reshape(abs(hilbert(fY')), 1, []);
        end
        
        for f = 1:nphase
            for g = 1:namp
                pac(f,g,i) = abs(sum(amp(g, :).*exp(sqrt(-1)*phase(f, :)))/(nsamples*ntrials));
                pac(f,g,i) = pac(f,g,i)./sqrt(amp(g, :)*amp(g, :)'/(nsamples*ntrials));
            end
        end
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

image(1).val     = squeeze(sum(sum(pac, 2), 1));
image(1).label   = ['pac_total_'  spm_file(D.fname, 'basename')];

c = 2;
for f = 1:nphase
    for g = 1:namp
        image(c).val     = squeeze(pac(f, g, :));
        image(c).label   = ['pac_phase_' num2str(S.phasefreq(f)) 'Hz_amp_'...
            num2str(S.ampfreq(g)) 'Hz_' spm_file(D.fname, 'basename')];
        c = c+1;
    end
end

res = image;