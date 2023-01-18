%% Script to train decoders and run sequenceness analysis on compound stimuli

%% Train decoders - compound stimuli during PRE localizer
%test during initial rest and mid rest period

do_decoder_training = 1; %train classifiers (1), or just reuse (0)

addpath(genpath(is.spm_path)) 
whenInMS = is.tss*is.msPerSample;
is.whichTimes = 1:length(is.tss);

iSj = is.nSubj;
start_sub = 1;%specify first subject
n_subjects = 11;%is.nSubj;%specify last subject

GoodChannel = nan(is.nSubj,is.numChan,length(find(strcmp(is.MEGruns{iSj},'loc1')))); %nsubject*nsensors*nruns
num_trial = nan(is.nSubj, length(find(strcmp(is.MEGruns{iSj},'loc1'))));         
num_tpts = nan(is.nSubj, length(find(strcmp(is.MEGruns{iSj},'loc1'))));           

%set rng seed for reproducibility
rng(101)

% Good Channels
for iSj = start_sub:n_subjects
    
    
    all_sess = find(strcmp(is.MEGruns{iSj},'loc1'));
    
    for isession = 1:length(all_sess)
        
        %specify folder and load preprocessed data from respective run
        tempdir = [is.OPTPath is.fnSID{iSj} filesep num2str(all_sess(isession),'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];
        load(fullfile(tempdir, 'opt.mat'));
        
        if isempty (opt.results.spm_files_epoched{1,1})
            SessionMEG = spm_eeg_load(fullfile(tempdir,[opt.results.spm_files_basenames{1,1}]));
        else
            SessionMEG = spm_eeg_load(fullfile(tempdir,[opt.results.spm_files_epoched_basenames{1,1}]));
        end
        
        chan_inds = indchantype(SessionMEG,'meeg','GOOD');
        Channindex  = indchantype(SessionMEG,'meeg','ALL');
        [~,ind] = setdiff(Channindex,chan_inds);
        
        index = ones(1,is.numChan);
        index(ind) = 0;
        GoodChannel(iSj,:,isession) = index;
        clear opt;
        
        num_trial(iSj, isession) = size(SessionMEG, 3);
        num_tpts(iSj, isession) = size(SessionMEG, 2);
        
    end
end

%train on pre localizer, define number of stimuli - differs with respect
%to condition 1 or 2
if mod(iSj,2) == 1
    nstates = 12;
else
    nstates = 14;
end

for iSj = start_sub:n_subjects   
    
    tic
    
    %load behavioral data to get presented stimulus names
    is.all_markov = [];
    
    is.all_stim_onsets_PRE = [];
    
    clear all_stim_onsets_PRE_idx
    clear temp_stim_PRE
    clear temp_stim_PRE_names
    clear unique_stim_PRE_names
    
    %define triggers/stimulus onsets    
    load(['behav/markov_localizer_PRE_LOG_', is.fnSID{iSj},'.mat'])
    
    all_stim_onsets_PRE_idx = find(cell2mat(markov_localizer_PRE.trigger.stimulusnameMEGLog(:,4)) == 32);
    is.all_stim_onsets_PRE = markov_localizer_PRE.trigger.stimulusnameMEGLog(all_stim_onsets_PRE_idx,2:3);
        
    is.all_markov = load(['behav/all_data/markov_LOG_All', is.fnSID{iSj},'.mat']);
    
    
    temp_stim_PRE = is.all_stim_onsets_PRE;
    
    %strip .png and create joint stimulus name
    temp_stim_PRE = strip(temp_stim_PRE,'right','g');
    temp_stim_PRE = strip(temp_stim_PRE,'right','n');
    temp_stim_PRE = strip(temp_stim_PRE,'right','p');
    temp_stim_PRE = strip(temp_stim_PRE,'right','.');
    
    %get stimulus names
    for i = 1:size(temp_stim_PRE,1)
        if strcmp(temp_stim_PRE(i,2), 'None')
            temp_stim_PRE_names(i,1) = temp_stim_PRE(i,1);
        else
            temp_stim_PRE_names(i,1) = strcat(temp_stim_PRE(i,1), ' /', {' '},temp_stim_PRE(i,2));
        end
    end
    
    unique_stim_PRE_names = unique(temp_stim_PRE_names);
    unique_stim_PRE_names = unique_stim_PRE_names(21:end);
    
    
    % creates gnf, which is a cell array of structs, with dimensions of 1 * trainTime * n states * L1p
    gnf = cell(1, length(is.tss), nstates, length(is.ENL1));
    
    %training on compound images during pre localizer
    lciInd = find(ismember(is.MEGruns{iSj},'loc1'));
    
    
    % avoiding hard coding (this assumes that all FL sessions for a subj have same number of trials and identical epoching)
    temp_num_tpts = min(num_tpts(iSj, :));
    temp_num_trials = max(num_trial(iSj, :));
    condition_log = [];
    data = nan(is.numChan, temp_num_tpts, temp_num_trials, length(lciInd));  % [channel, epoch_time_points, trials, sessions]
    trialindex = zeros(temp_num_trials, length(lciInd));
    
    % Getting the data
    for ii = 1:length(lciInd)
        
        % load the epoched data
        dir_temp = [is.OPTPath is.fnSID{iSj} filesep num2str(lciInd(ii),'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];
        
        opt = load(fullfile(dir_temp, 'opt.mat'));
        opt = opt.opt;
        DLCI = spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_epoched_basenames{1,1}]));
        
        % load good channels:
        chan_MEG = indchantype(DLCI,'meeg');
        
        % load good trials
        allconds = unique_stim_PRE_names';
        
        all_trls = sort([DLCI.indtrial(allconds(:))]);
        good_trls = sort([DLCI.indtrial(allconds(:),'good')]);
        good_ind = ismember(all_trls,good_trls);
        good_ind = all_trls(good_ind);
        
        % change this for subject 512, run 2, first 42 trials not recorded
        if sum(is.fnSID{iSj} == '512') == 3 && ii == 2
            trialindex(good_ind+42,ii)=1;
            % get clean data
            cldata = DLCI(chan_MEG,:,good_trls);
            data(:,:,good_ind+42,ii)= cldata;
        else
            trialindex(good_ind,ii)=1;
            % get clean data
            cldata = DLCI(chan_MEG,:,good_trls);
            data(:,:,good_ind,ii)= cldata;
        end
        
        % get labels to use below
        condition_log = [condition_log; DLCI.conditions'];
    end
    
    clear Cleandata
    clear Goodtrialindex
    clear lcidata
    
    % reshape
    Cleandata = reshape(data,[size(data,1),size(data,2),size(data,3)*size(data,4)]); % combine all runs of data;
    Goodtrialindex = reshape(trialindex,[size(trialindex,1)*size(trialindex,2),1]);  % combine all runs of good trials;
    
    clear subsetindex
    clear labStim
    
    subsetindex = Goodtrialindex;
    for k = 1:size(condition_log,1)
        if strlength(condition_log(k)) > 5
            labStim(k,1) = find(strcmp(condition_log(k),allconds)); % converts 'S1' --> '1'
        else
            labStim(k,1) = NaN;
        end
    end
    
    Cleandata_subset = Cleandata(:,:,logical(subsetindex)); % good trials
    
    if sum(is.fnSID{iSj} == '512') == 3
        subsetindex = [subsetindex(1:204); subsetindex(247:end)];
    end
    
    labStim_subset = labStim(logical(subsetindex));
    stimlabel = labStim_subset;
   
    sjchannel = squeeze(GoodChannel(iSj,:,:));
    goodchannindex = nansum(sjchannel') == length(all_sess);
    lcidata = Cleandata_subset(logical(goodchannindex),:,:);
    
    
    if do_decoder_training == 1
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare Dataset(traindata, nulldata) for training classifer
        for iTT = 1:length(is.whichTimes) % which specific time to be trained on (180 ms)
            nullData=[];
            trueData=[];
            % subtract ISI to get a baseline
            for ii=1:size(lcidata,3)
                nullData(:,ii) = squeeze(lcidata(:,is.picOnsetSlice - round(is.nullslice/is.msPerSample),ii)); % nullsclice
                trueData(:,ii) = squeeze(lcidata(:,is.picOnsetSlice + is.tss(is.whichTimes(iTT)),ii));
            end
            
            
           nullData = scaleFunc(nullData');  
           trueData = scaleFunc(trueData');

            
            clear labels
            % Set up the classifer
            for iShuf=1
                
                for iC=1:size(allconds,2)   % train classifiers on null and iT sample of data.
                    labels = [stimlabel == iC; zeros(size(nullData,1),1)];
                    if iShuf==1
                        labels = labels;
                    else
                        labels = labels(randperm(length(labels)));
                    end
                    
                    for iL1=1:length(is.ENL1)  % loop over L1 penalty values: only one value specified apriori
                        l1p = is.ENL1(iL1); l2p = 0; % alpha = L1/(2*L2+L1) ; lambda = 2*L2+L1 --> L2 penalty set to 0 to avoid elastic net
                        [beta, fitInfo] = lassoglm([trueData; nullData], labels, 'binomial', ...
                            'Alpha', l1p / (2*l2p+l1p), 'Lambda', 2*l2p + l1p, 'Standardize', false);
                        gnf{iShuf,is.whichTimes(iTT),iC,iL1}.beta = beta; gnf{iShuf,is.whichTimes(iTT),iC,iL1}.Intercept = fitInfo.Intercept;
                    end
                end
            end
            
        end
        
        save([is.OutputDir 'trained_decoders_compounds/' is.fnSID{iSj}], 'gnf', 'whenInMS', 'is')

        disp(['sj' num2str(iSj) ': No. of good trials = ' num2str(sum(subsetindex))]);
        toc;
        
    end
end
    

%% Decode Resting 
for iRun = 1:2 % 2 rest sessions, before PRE localizer and between prior learning and transfer phase in main experiment
    is.nShuf=1;

    for iSj= start_sub:n_subjects
        
        %load behavioral data to get presented stimulus names
        is.all_markov = [];
        
        all_stim_onsets_PRE_idx = [];
        is.all_stim_onsets_PRE = [];
        
        clear all_stim_onsets_PRE_idx
        clear temp_stim_PRE
        clear temp_stim_PRE_names
        clear unique_stim_PRE_names
        
        %define triggers/stimulus onsets        
        load(['behav/markov_localizer_PRE_LOG_', is.fnSID{iSj},'.mat'])
        load(['behav/markov_localizer_POST_LOG_', is.fnSID{iSj},'.mat'])
        
        all_stim_onsets_PRE_idx = find(cell2mat(markov_localizer_PRE.trigger.stimulusnameMEGLog(:,4)) == 32);
        is.all_stim_onsets_PRE = markov_localizer_PRE.trigger.stimulusnameMEGLog(all_stim_onsets_PRE_idx,2:3);
                
        is.all_markov = load(['behav/all_data/markov_LOG_All', is.fnSID{iSj},'.mat']);
        
        temp_stim_PRE = is.all_stim_onsets_PRE;
           
        unique_stim_PRE_tmp_names = unique(temp_stim_PRE(:,1));
        unique_stim_PRE_tmp_names = unique_stim_PRE_tmp_names(21:end);
        
       
        %strip .png and create joint stimulus name
        temp_stim_PRE = strip(temp_stim_PRE,'right','g');
        temp_stim_PRE = strip(temp_stim_PRE,'right','n');
        temp_stim_PRE = strip(temp_stim_PRE,'right','p');
        temp_stim_PRE = strip(temp_stim_PRE,'right','.');

        %get stimulus names
        for i = 1:size(temp_stim_PRE,1)
            if strcmp(temp_stim_PRE(i,2), 'None')
                temp_stim_PRE_names(i,1) = temp_stim_PRE(i,1);
            else
                temp_stim_PRE_names(i,1) = strcat(temp_stim_PRE(i,1), ' /', {' '},temp_stim_PRE(i,2));
            end
        end

        unique_stim_PRE_names = unique(temp_stim_PRE_names);
        unique_stim_PRE_names = unique_stim_PRE_names(21:end);
        
        
        allconds = unique_stim_PRE_names';
    
        % load classifier
        gnf = load([is.OutputDir 'trained_decoders_compounds/' is.fnSID{iSj} '.mat'], 'gnf') ; % get the gnf variable with the regression models
        gnf = gnf.gnf;

        % Load Resting State Data
        dstInd1 = find(ismember(is.MEGruns{iSj},'iniRest'));
        dstInd2 = find(ismember(is.MEGruns{iSj},'midRest'));
        dstInd = [dstInd1 dstInd2];
        
        dir_temp = [is.OPTPath is.fnSID{iSj} filesep num2str(dstInd(iRun),'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];   
        opt = load(fullfile(dir_temp, 'opt.mat'));
        opt = opt.opt;    
        RST = spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_basenames{1,1}]));

        % Good timepoints/trials
        good_samples = ~all(badsamples(RST,':',':',':'));

        % Select MEG channel:
        chan_meg = indchantype(RST,'meeg');
        rRST = RST(chan_meg,:,:);

        sjchannel=squeeze(GoodChannel(iSj,:,:));
        goodchannindex = nansum(sjchannel') == sum(cellfun(@(x) ~isempty(x), num2cell(all_sess)));

        rRST = rRST(logical(goodchannindex),:,:);

        % get the interval between the start and the end
        evtypes = RST.events; evtypes = {evtypes(:).type}';
        evvals = RST.events; evvals = {evvals(:).value}';
        evvals(all(cellfun(@ischar,evvals),2),:) = {0}; % replace string with 0
        evvals = cell2mat(evvals); %convert cell to double
        evtimes = RST.events; evtimes = [evtimes(:).time]';
        
        if strcmp(is.MEGruns{iSj}{dstInd(iRun)},'iniRest')
            %initial rest
            RestingTriggerVals = is.rest_boundaries(1);  % defined in initialize script           
            EndTriggerVals = is.rest_boundaries(2);  %defined in initialize script
            
        elseif strcmp(is.MEGruns{iSj}{dstInd(iRun)},'midRest')
            %middle rest
            RestingTriggerVals = is.rest_boundaries(3);  % defined in initialize script           
            EndTriggerVals = is.rest_boundaries(4);  %defined in initialize script
        end
        
        RestStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, RestingTriggerVals);    % onset of Resting
        RestStmTimes = evtimes(RestStmInds);
        
        EndStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, EndTriggerVals);  % end of Resting
        EndStmTimes = evtimes(EndStmInds);
        
        
        wholetimeindex = zeros(size(RST,2),1);

        wholetimeindex(ceil(RestStmTimes(end)*(is.sampleRate/is.smoothFact)) : floor(EndStmTimes(end)*(is.sampleRate/is.smoothFact)))=1; % taken from MMN, avoid hard coding

        % set the artifacts to zero
        badwithin_index = wholetimeindex == 1 & double(good_samples') == 0; 
        rRST(:,badwithin_index) = 0;
        data = rRST(:,logical(wholetimeindex));

        data = scaleFunc(data');    

        % apply classifier to the clean data
        if isempty(RestStmTimes) 
            nTr = 1;
        else
            nTr = length(RestStmTimes);
        end

        if do_decoder_training == 1
            
            Rreds = cell(is.nShuf,length(is.whichTimes),nTr,length(is.ENL1));
            
            for whichShuf=1:is.nShuf
                for iTT = 1:length(is.whichTimes)
                    iT=is.whichTimes(iTT);
                    for iTr=1:nTr
                        for iL1=1:length(is.ENL1)
                            for iC=1:size(allconds,2)
                                % use regression models --> use sigmoid
                                % transformed predictions
                                Rreds{whichShuf,iT,iTr,iL1}(:,iC) = 1 ./ (1 + exp(-(data * gnf{whichShuf,iT,iC,iL1}.beta + gnf{whichShuf,iT,iC,iL1}.Intercept)));   % matlab's built-in lassoglm
                            end
                        end
                    end
                end
            end
            
            save([is.OutputDir 'decoded_rest_compounds/' is.fnSID{iSj} '_' num2str(iRun)], 'Rreds', 'is', '-v7.3')% save, as long as only using single subject data

        end
        
        disp(['sj' num2str(iSj) ' finished']);    
        
    end
    
    %% TDLM
    for iSj = start_sub:n_subjects 
        
       %% Prepare transition matrices for TDLM        
        load(['behav/markov_localizer_PRE_LOG_', is.fnSID{iSj},'.mat'])
        load(['behav/markov_localizer_POST_LOG_', is.fnSID{iSj},'.mat'])
        
        all_stim_onsets_PRE_idx = find(cell2mat(markov_localizer_PRE.trigger.stimulusnameMEGLog(:,4)) == 32);
        is.all_stim_onsets_PRE = markov_localizer_PRE.trigger.stimulusnameMEGLog(all_stim_onsets_PRE_idx,2:3);
                
        is.all_markov = load(['behav/all_data/markov_LOG_All', is.fnSID{iSj},'.mat']);
        
        temp_stim_PRE = is.all_stim_onsets_PRE;
        
        unique_stim_PRE_tmp_names = unique(temp_stim_PRE(:,1));
        unique_stim_PRE_tmp_names = unique_stim_PRE_tmp_names(21:end);
        
        %find transition structure that has been used to produce observed sequences
        prior_learning_raw = is.all_markov.all_learning_trials_prior;
        
        all_states = unique_stim_PRE_tmp_names';
        
        %translate to numbers
        for k = 1:size(prior_learning_raw,1)
            schedule_prior(k,1) = find(strcmp(prior_learning_raw(k), all_states));
        end
        
        %make transition matrices from schedule
        state_visit_counts = zeros(max(schedule_prior),max(schedule_prior));
        
        for i = 1:size(schedule_prior,1)-1
            
            current_state = schedule_prior(i);
            next_state = schedule_prior(i+1);
                        
            %keep count of specific transition
            state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
        end
        
        %make TM
        TM_prior = state_visit_counts ./ sum(state_visit_counts,2);
        %check if all rows of TM sum to 1
        if sum(nansum(TM_prior,2) == ones(size(TM_prior,1),1)) ~= size(TM_prior,1)
            disp('TM_prior row does not add to 1!')
        end
      
        uniquePerms = [];
        
        for iter = 1:1000
            unique_perms_tmp = [];

            %create random permutation
            rp = randperm(size(TM_prior,1));
            unique_perms_tmp = TM_prior(rp,rp);
            
            %concatenate current random TM to all random TM
            uniquePerms = cat(3, uniquePerms, unique_perms_tmp);
            
            all_diag_diff(iter) = sum(diag(uniquePerms(:,:,iter)) - ones(size(uniquePerms(:,:,iter),1),1));
        end
        
        if sum(all_diag_diff ~= size(TM_prior,1) * -1) > 0
            disp('random TMs have self transitions!')
        end

        uniquePerms = cat(3, TM_prior, uniquePerms);
        
        is.nShuf = size(uniquePerms,3); % 1st ID is ground truth
       
        sf_All = nan(length(is.lgncy)+1, is.nShuf, is.nSubj, 1, length(is.whichTimes), length(is.ENL1));
        
        Tfwd = TM_prior;

    %% Compute Sequenceness measure, using TDLM

        maxTrials = 1; 

        sf = cell(is.nShuf, maxTrials);
        sb = cell(is.nShuf, maxTrials);
        
        %load decoding data during rest
        S = load([is.OutputDir 'decoded_rest_compounds/' is.fnSID{iSj} '_' num2str(iRun)]);
        
        Rreds = S.Rreds;  % load this subject's preds

        nTr = size(Rreds,3);
        L1l = length(is.ENL1); 
        Ltt = length(is.whichTimes);
        maxLag = length(is.lgncy);

        %define number of states - differs with respect
        %to condition 1 or 2
        nstates = size(TM_prior,1);

        tic
            for iTr = 1:nTr
                tmp = squeeze(shiftdim(Rreds(1,:,iTr,:),-1));
                
                mpty = all(all(cellfun(@isempty, tmp),2),3); tmp(mpty,:,:) = [];  % remove empty trainingTimes
                prSj = permute(cell2mat(shiftdim(tmp,-2)),[1 2 4 3]); % this is samples*states*alpha*lambda*trainingTimes
                prSj = prSj(:,:,:); % alpha and trainingTimes put into a single dimension for the vectorization of sequenceness
                
                nP = size(prSj,3);
                sf_temp = nan(maxLag,is.nShuf,nP);
                sb_temp = nan(maxLag,is.nShuf,nP);
                
                for vec = 1:nP
                    
                    X = squeeze(prSj(:,:,vec));
                    
                    if ~any(~isnan(X(:)))
                        continue
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    nbins = maxLag+1;
                    
                    warning off
                    dm = [toeplitz(X(:,1),[zeros(nbins,1)])]; %create design matrix
                    dm = dm(:,2:end);
                    
                    %concatenate design matrices for all states
                    for kk = 2:nstates
                        temp = toeplitz(X(:,kk),[zeros(nbins,1)]);
                        temp = temp(:,2:end);
                        dm = [dm temp];
                    end
                    
                    warning on
                    
                    Y = X;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    betas = nan(nstates*maxLag, nstates);
                    bins = 10;
                    
                    %run temporally delayed linear models - predicting state
                    %reactivation in next state based on previous state at
                    %different temporal lags
                    for ilag = 1:bins%maxLag
                        zinds = (1:maxLag:nstates*maxLag) + ilag - 1;
                        temp_zinds = (1:bins:nstates*maxLag) + ilag - 1;
   
                        temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y;
                        betas(temp_zinds,:) = temp(1:end-1,:);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %compare transition matrix of beta weights of states with
                    %ground truth transition matrix
                    X = [];
                    betasr = reshape(betas,[maxLag nstates nstates]);
                    betasnbins_statessquared = reshape(betas,[maxLag nstates^2]);
                    
                    for iShuf = 1:is.nShuf
                        
                        if iShuf == 1  % the real transition matrix
                            T1 = Tfwd;
                            T2 = Tfwd';
                        else  % use randomly shuffled transition matrix
                            Tfwd_temp = uniquePerms(:,:,iShuf);
                            T1 = Tfwd_temp; %forward transitions
                            T2 = T1'; % backwards transitions
                        end
                        
                        bbb = pinv([T1(:) T2(:) squash(eye(nstates)) squash(ones(nstates))])*(betasnbins_statessquared');
                        sf_temp(:,iShuf,vec) = bbb(1,:);
                        sb_temp(:,iShuf,vec) = bbb(2,:);
                    end
                end
                
                for iShuf = 1:is.nShuf
                    sf{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l, Ltt);
                    sb{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l, Ltt);
                    
                    sf{iShuf,iTr}(1, 1, 2:end, :, :) = reshape(squeeze(sf_temp(:,iShuf,:)),[maxLag L1l Ltt]);
                    sb{iShuf,iTr}(1, 1, 2:end, :, :) = reshape(squeeze(sb_temp(:,iShuf,:)),[maxLag L1l Ltt]);
                    
                end
                
            end
            
        % RESHAPE
        % aim is -> latency(61) * shuffles(1000) * trials (1) * trainTimes (180ms) * L1 param
        sf2 = permute(cell2mat(sf), [3 1 2 5 4]); % sf2 = latency(61) * shuffles(1000) * trials (1) * trainTimes (180ms) * L1 param
        sf2 = cat(3, sf2, nan(length(is.lgncy)+1, is.nShuf, maxTrials - size(sf2,3), Ltt, L1l));  % [pad with nans for trials after the end]
        sb2 = permute(cell2mat(sb), [3 1 2 5 4]);
        sb2 = cat(3, sb2, nan(length(is.lgncy)+1, is.nShuf, maxTrials - size(sb2,3), Ltt, L1l));  % [pad with nans for trials after the end]

        sfAll(:,:,iSj,:,:,:) = sf2; 
        sbAll(:,:,iSj,:,:,:) = sb2;

        disp(['Sub' num2str(iSj)  ' Sequence ' num2str(iRun) ' Finished' ])  
        toc
    end
        
    save([is.OutputDir 'resting_results/'  'StimAll_rst' num2str(iRun-1)] ,'sfAll','sbAll','is','-v7.3');

end

%% Plot replay effect, get p-value
meg_summarise_and_plot_compounds_prereg


