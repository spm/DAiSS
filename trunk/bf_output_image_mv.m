function res = bf_output_image_mv(BF, S)
% Computes multivariate test on a number of frequency bands
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes, modified from Vladimir Litvak's example code
% $Id$

%--------------------------------------------------------------------------

%% no covariance matrix creation -so need to check bands are within this window

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
    
    
    
    %     design = cfg_const;
    %     design.tag = 'design';
    %     design.name = 'design';
    %     design.help = {'Use default settings for the inversion'};
    %     design.val  = {1};
    
    design = cfg_files;
    design.tag = 'design';
    design.name = 'design matrix';
    design.filter = 'mat';
    design.num = [1 Inf];
    design.help = {'Select the design matrix'};
    
    
    
    woi = cfg_entry;
    woi.tag = 'woi';
    woi.name = 'Time windows of interest';
    woi.strtype = 'r';
    woi.num = [Inf 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time windows'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.val = {[-Inf Inf]};
    foi.help = {'Freq bands'};
    
    
    datafeatures = cfg_menu;
    datafeatures.tag = 'datafeatures';
    datafeatures.name = 'Features';
    datafeatures.labels =get_data_features;
    datafeatures.values =get_data_features;
    
    datafeatures.help = {'Data features of interest'};
    datafeatures.val =datafeatures.values(1);
    
    contrast = cfg_entry;
    contrast.tag = 'contrast';
    contrast.name = 'Time contrast';
    contrast.strtype = 'i';
    contrast.num = [1 Inf];
    contrast.val = {1};
    
    result         = cfg_menu;
    result.tag     = 'result';
    result.name    = 'What to output';
    result.help    = {'Specify output type'};
    result.labels  = {
        'chi square'
        'BIC'
        'r square'
        }';
    result.values  = result.labels;
    result.val = {'chi square'};
    
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
    
    custom = cfg_branch;
    custom.tag = 'custom';
    custom.name = 'Custom';
    custom.help = {'Define custom settings for the inversion'};
    custom.val  = {whatconditions, contrast, woi};
    
    isdesign = cfg_choice;
    isdesign.tag = 'isdesign';
    isdesign.name = 'Design matrix parameters';
    isdesign.help = {'Choose whether to load custom design'};
    isdesign.values = {design, custom};
    isdesign.val = {design};
    
    
    image_mv      = cfg_branch;
    image_mv.tag  = 'image_mv';
    image_mv.name = 'Mv image';
    image_mv.val  = { isdesign, datafeatures, foi,  result, modality};
    
    res = image_mv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end




%keyboard;
D = BF.data.D;


if isfield(S.isdesign,'whatconditions'),
    %% gui specified conditions and contrast
    samples = {};
    for i = 1:size(S.woi, 1)
        samples{i} = D.indsample(S.woi(i, 1)):D.indsample(S.woi(i, 2));
    end
    
    if isfield(S.whatconditions, 'all')
        trials{1} = 1:D.ntrials;
    else
        for i = 1:numel(S.whatconditions.condlabel)
            if isempty(D.indtrial(S.whatconditions.condlabel{i}, 'GOOD'))
                error('No trials matched the selection.');
            end
            trials{i} = D.indtrial(S.whatconditions.condlabel{i}, 'GOOD');
        end
        if isempty(trials)
            error('No trials matched the selection, check the specified condition labels');
        end
    end
    
else %%  conditions and contrast  specified in a file
    if ~exist(cell2mat(S.isdesign.design)),
        error('Cannot load design matrix');
    end;
    disp('loading design matrix');
    a=load(cell2mat(S.isdesign.design));
    X=a.design.X; %% design matrix
    contrast=a.design.contrast;
    ntrials=size(X,1);
    if (size(a.design.Xstartlatencies,1)~=ntrials)||(size(a.design.Xtrials,1)~=ntrials) 
        error('start latencies and Xtrials and X should have a value per row of the design');
    end;
    
    trials=a.design.Xtrials; %% indices of trials to use
    for j=1:ntrials,
        allsamples(j,1)=D.indsample(a.design.Xstartlatencies(j));
        allsamples(j,2)=D.indsample(a.design.Xstartlatencies(j)+a.design.Xwindowduration);
    end;
    
end;



goodchannels = D.indchantype(S.modality, 'GOOD'); %% THIS SHOULD BE DONE IN THE PRE PROC STAGE

W = BF.inverse.W.(S.modality);
if numel(W{1})==length(goodchannels)
    %% no montage used
    U=eye(length(goodchannels));
    chanind=goodchannels;
else
    if numel(W{1})==size(BF.features.(S.modality).montage.U,2),
        U=BF.features.(S.modality).montage.U;
        chanind=BF.features.(S.modality).montage.chanind;
        disp('Picking up data reduction montage');
        
        %% data reduction montage was applied earlier
    else
        error('size of weights does not match number of (reduced) channels');
    end;
end;

Nchans=size(U,2); %% effective number of channels

YY       = {};
nsamples = unique(allsamples(:,2)-allsamples(:,1));

if length(nsamples) > 1
    error('All time windows should be equal lentgh')
end

alltrials = spm_vec(trials);
ntrials   = length(alltrials);


%% now identify frequency bands of interest

Nbands=size(S.foi,1);
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

% Tfull      = dctT(:,allfreqind); %% A filter for all bands (not necessarily continuous)
%% end of freq band section

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials , 'Computing covariance'); drawnow;
if ntrials  > 100, Ibar = floor(linspace(1, ntrials ,100));
else Ibar = 1:ntrials; end

%% load in data and make up simple design matrix
%ncond=numel(samples); %% number of conditions= columns in design matrix
%Nt=ntrials*ncond; %% total number of time windows under consideration
%X=zeros(Nt,ncond);

flatdata=zeros(ntrials*nsamples,Nchans);
%% want flatdata in form Nchans,Nt*Nsamples

%count=0;
YY       = {};
for i = 1:ntrials
    %for j = 1:numel(samples)
     %   count=count+1;
     %   X(count,j)=1;
        Y  = U'*squeeze(D(chanind, allsamples(i,1):allsamples(i,2)-1, alltrials(i)));
        Y  = detrend(Y'); %% detrend and throw away low freq drift
        flatdata((i-1)*nsamples+1:i*nsamples,:) =Y;
        
    %end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

W = BF.inverse.W.(S.modality);
nvert = numel(W);
S.regressout=[]; %% turn off for now
regressout=S.regressout;


%% set up the data features
weights=-1; %% set up flag
Yfull=get_data_features(flatdata,nsamples,ntrials,weights,dctT,S.datafeatures,featureind,regressout); %% set up data structures



spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, 'Scanning grid points'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

outval = nan(1, nvert);


for i = 1:nvert
    
    if ~isnan(W{i})
        
        w    = W{i};
        
        %% returns columns of a matrix with rows as different observations
        [Yfull,vedata]=get_data_features(flatdata,nsamples,ntrials,w,dctT,S.datafeatures,featureind,regressout); %% extract the data features
        
        Yfull=Yfull-repmat(mean(Yfull),size(Yfull,1),1); %% remove dc level from each column/feature
        Yfull=Yfull./repmat(std(Yfull),size(Yfull,1),1); %% normalize features to have unit variance by default
        [chival,BIC,cva] = output_image_mv_cva(X,Yfull,contrast); %% run the multivariate test
        
        
        switch S.result
            case 'chi square'
                resultstr='chisq';
                outval(i) = chival(1);
            case 'r square'
                outval(i) = cva.ccorr.^2;
                resultstr='rsq';
            case 'BIC'
                outval(i)=BIC(1);
                resultstr='BIC';
        end;      
        
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');


image(1).val   = outval;

image(1).label = ['mv' resultstr S.datafeatures freqstr  spm_file(D.fname, 'basename')];


res = image;