function res = bf_features_freq(BF, S)
% Simple band limited covariance computation with regularization
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    lambda = cfg_entry;
    lambda.tag = 'lambda';
    lambda.name = 'Regularisation';
    lambda.strtype = 'r';
    lambda.num = [1 1];
    lambda.val = {0};
    lambda.help = {'Select the regularisation (in %)'};
    
    bayespca = cfg_menu;
    bayespca.tag = 'bayespca';
    bayespca.name = 'Method of Bayes PCA';
    bayespca.help = {'Select Bayesian method to get regularization'};
    bayespca.labels = {'Minka trunc', 'Minca reg'};
    bayespca.values = {'Minka trunc', 'Minca reg'};
    bayespca.val = {'Minka trunc'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.val = {[-Inf Inf]};
    foi.help = {'Frequency windows within which to compute covariance over (sec)'};
    
    
    
    
    
    regmethod = cfg_choice;
    regmethod.tag = 'regmethod';
    regmethod.name = 'Regularization scheme';
    regmethod.help = {'Choose whether to load custom design'};
    regmethod.values = {lambda, bayespca};
    regmethod.val = {bayespca};
    
    freq      = cfg_branch;
    freq.tag  = 'freq';
    freq.name = 'Regularized band limited covariance';
    freq.val  = {foi,regmethod};
    
    res = freq;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;


ntrials = length(S.trials);
nchans=length(S.channels);
%% now identify frequency bands of interest

Nbands=size(S.foi,1);
Nwoi=numel(S.samples);

if length(unique(length([S.samples{:}])))~=1,
    error('all windows must be of equal length');
end;
nsamples=length(S.samples{1}); %% use length of first window to set up DCT (as all windows fixed at same length)
windowduration=(nsamples/D.fsample);
dctfreq = (0:nsamples-1)/2/windowduration;           % DCT frequencies (Hz)
dctT      = spm_dctmtx(nsamples,nsamples);

freqstr=[];
allfreqind=[];
for fband=1:Nbands, %% allows one to break up spectrum and ignore some frequencies
    freqrange=S.foi(fband,:);
    j      = find( (dctfreq >= freqrange(1)) & (dctfreq<=freqrange(2)) );
    featureind{fband}=j;
    allfreqind=sort(unique([allfreqind j]));
    freqstr=[freqstr sprintf('%3.1f-%3.1f,',dctfreq(min(j)),dctfreq(max(j)))];
end; % for fband=1:Nbands

% Hanning operator (if requested)
%----------------------------------------------------------------------
S.windowing='Hanning'; %% set by default for now
switch S.windowing,
    case 'Hanning',
        W  = repmat(spm_hanning(nsamples)',nchans,1);
    case 'None'
        W  = ones(nchans,nsamples);
end

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end





allY=[];
Tband=dctT(:,allfreqind); % filter to this band
for i = 1:ntrials
    for j = 1:numel(S.samples)
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));
        
        
        Y = detrend(Y', 'constant')';
        Y=Y.*W;
        dctY=Y*Tband; %% frequency representation
        
        %YY = YY+(dctY*dctY');
        allY=[allY dctY];
        
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


YY=allY*allY';

[U,alls] = svd(YY);

cumpower=cumsum(diag(alls))./sum(diag(alls));
nmodes99=min(find(cumpower>0.99));


spm_progress_bar('Clear');

%% now for the regularisation

Mopt=[];
if isfield(S.regmethod,'bayespca'),
    switch S.regmethod.bayespca,
        case 'Minka trunc',
            %%% WIll's code- based on MInka- based on trunctation to optimum model order
 
 
             [M_opt,log_ev,lambda1] = spm_pca_order (allY);
             disp(sprintf('Estimated covariance matrix order %d',M_opt));
             noisevar=mean(lambda1(M_opt+1:max(find(lambda1>0)))); %% some eigenvals can come out negative
             redYY=U(:,1:M_opt)*alls(1:M_opt,1:M_opt)*U(:,1:M_opt)';
             C=redYY; %% rank deficient but with full complement of channels
             Cr=U(:,1:M_opt)'*YY*U(:,1:M_opt); %% compact version of the covariance matrix
         
             
             
             noise_id=eye(size(allY,1))*noisevar; %% noise power
        case 'Minka reg'  %% use Will's bayes pca to get noise level then use this to augment diagonal
            error('not defined yet');
             disp('Will Penny Bayes PCA to get noise and augment diagonal');
             [M_opt,log_ev,lambda1] = spm_pca_order (allY);
             noisevar=mean(lambda1(M_opt+1:max(find(lambda1>0)))); %% some eigenvals can come out negative
             
             opstats.eff_reg=100*noisevar./mean(allsvd);
             disp(sprintf('effective regularisation =%3.2f percent',opstats.eff_reg));
             C=YY+eye(size(YY,1))*noisevar;
             
             noise_id=eye(size(YY,1))*noisevar;
    end; %% switch
else
    %% manual choice of reg value
    
    lambda = (S.lambda/100) * trace(YY)/size(YY,1);
    C= YY + lambda * eye(size(YY));    
end;

Cinv=pinv(C);

features=[];

features.C=C;
features.Cinv=Cinv;


%BF.data.Cinv=Cinv; %% write inverse cov to data structure
%BF.data.C=C; %% write cov to data structure


if ~isempty(M_opt),
    
    
    features.montage.chanind=S.channels;
    features.montage.U=U(:,1:M_opt);
    features.montage.C=Cr;
    features.montage.Cinv=pinv(Cr);
end;
    


res = features;