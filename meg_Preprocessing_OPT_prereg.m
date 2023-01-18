%% Preprocessing loop
%Load data (neural and behavioral files)
%correction of onsets/trigger events to recorded photodiode times
%Crop unused data (not within first and last trial trigger in run)
%Filtering
%Downsampling to 100 Hz and removal of bad segments
%ICA, regressing out ICs correlated with EOG (eyetracker) signals
%Epoching of trials (not in resting blocks) 
%Assignment of trial types to trigger events
%Coregistration

for iSj = 11:is.nSubj 

    is.fnDate{iSj} = is.fnMEG{iSj}(end-7:end-1); %date of raw MEG data
    spL = [is.prereg_path filesep is.OPTPath is.fnSID{iSj} filesep];   % output path 
    spN = [is.rawMEG filesep is.fnSID{iSj} filesep is.fnMEG{iSj}];   % raw data path
    
    is.all_markov = [];
    
    all_stim_onsets_PRE_idx = [];
    is.all_stim_onsets_PRE = [];
    
    all_stim_onsets_POST_idx = [];
    is.all_stim_onsets_POST = [];
    
    %define triggers/stimulus onsets
    load([is.behav_path, filesep, 'markov_localizer_PRE_LOG_', is.fnSID{iSj},'.mat'])
    load([is.behav_path, filesep, 'markov_localizer_POST_LOG_', is.fnSID{iSj},'.mat'])
    load([is.behav_path, filesep, 'markov_LOG_', is.fnSID{iSj},'.mat'])

    all_stim_onsets_PRE_idx = find(cell2mat(markov_localizer_PRE.trigger.stimulusnameMEGLog(:,4)) == 32);
    is.all_stim_onsets_PRE = markov_localizer_PRE.trigger.stimulusnameMEGLog(all_stim_onsets_PRE_idx,2:3);

    all_stim_onsets_POST_idx = find(cell2mat(markov_localizer_POST.trigger.stimulusnameMEGLog(:,4)) == 32);
    is.all_stim_onsets_POST = markov_localizer_POST.trigger.stimulusnameMEGLog(all_stim_onsets_POST_idx,2);

    all_stim_onsets_main_idx = find(cell2mat(markov.trigger.stimulusnameMEGLog(:,4)) == 32);
    is.all_stim_onsets_main = markov.trigger.stimulusnameMEGLog(all_stim_onsets_main_idx,2:3);

    is.all_markov = load([is.behav_path, filesep, 'all_data', filesep, 'markov_LOG_All', is.fnSID{iSj},'.mat']);
    
 
    for iR = [5:6 11] %[1:6 11] %rest periods and localizer I runs
        if ~isempty(is.MEGruns{iSj}{iR})           
            inSt = ['_', num2str(iR, '%02d')];

            localPath = [spL inSt(2:end) '.ds'];
            mkdir(localPath);
            ds_folder = [spN filesep is.fnMEG{iSj} inSt '.ds'];         
                                
            %% Import - CTF resting state using OPT
            S = struct;
            S.fif_file = ds_folder;
            S.spm_file = fullfile(localPath,'spmeeg.mat');
            
            %find other channels (eyetracker + photodiode + trigger channel) in montage
            S.other_channels = {'UADC001','UADC002','UADC003','UADC004', 'UPPT001'};
            
            D = osl_import(S.fif_file, 'outfile', S.spm_file, 'other_channels', S.other_channels); 

            %% add up other channels
            D = D.chantype(find(strcmp(D.chanlabels,'UADC001')),'EOG1');
            D = D.chantype(find(strcmp(D.chanlabels,'UADC002')),'EOG2');
            D = D.chantype(find(strcmp(D.chanlabels,'UADC003')),'EOG3');
            D = D.chantype(find(strcmp(D.chanlabels,'UADC004')),'photodiode');
            D = D.chantype(find(strcmp(D.chanlabels,'UPPT001')),'triggers');

            D.save()

            %% replace trigger times with photodiode times
            if strcmp(is.MEGruns{iSj}{iR},'loc1')
                           
                 load([localPath, filesep, 'spmeeg.mat']);

                 %find trigger onsets
                 trigger_onsets_tmp = [];
                 trigger_onsets = [];

                 trigger_onsets_tmp = find(diff(D.data(279,:) == 32));
                 if sum(diff(trigger_onsets_tmp) > 10 & diff(trigger_onsets_tmp) < 20) == length(trigger_onsets_tmp)/2
                    trigger_onsets(1,:) = trigger_onsets_tmp(1:2:end) + 1;
                    trigger_onsets(2,:) = trigger_onsets_tmp(2:2:end) + 1;
                 else
                     disp('triggers do not add up!')
                     keyboard
                 end
                 
                 %interval of trigger
                 trigger_onsets(3,:) = trigger_onsets(2,:) - trigger_onsets(1,:);
                 
                 %find photodiode change points 
                 diode_onsets_tmp = [];
                 diode_onsets = [];

                 diode_onsets_tmp = find(ischange(D.data(278,:),'Threshold',1));

                 
                 %find diode change point closest in time to stimulus triggers
                 for k = 1:size(diode_onsets_tmp,2)
                     current_idx = [];
                     %only if reasonable, 200 sample points after trigger
                     %onset
                     if numel(find((diode_onsets_tmp(k) - trigger_onsets(1,:)  < 200) ...
                                 & (diode_onsets_tmp(k) - trigger_onsets(1,:)  >= 0))) > 0 
                             
                         current_idx = find((diode_onsets_tmp(k) - trigger_onsets(1,:) < 200) ...
                                          & (diode_onsets_tmp(k) - trigger_onsets(1,:) >= 0));

                         trigger_onsets(4,current_idx) = diode_onsets_tmp(k) - 1;
 
                     end
                    
                 end
                 
                trigger_onsets(5,:) = trigger_onsets(4,:) + trigger_onsets(3,:);

                %save old triggers                
                old_triggers = [];
                old_triggers = D.data(279,:);
                save([localPath, filesep, 'old_triggers.mat'], 'old_triggers')
                
                %delete old triggers and replace
                old_trigger_idx = [];
                old_trigger_idx = find(D.data(279,:) == 32);
                D.data(279,old_trigger_idx) = 0;
                
                for k = 1:size(trigger_onsets,2)
                    diode_intervall = [];
                    %specify diode intervall and replace stimulus trigger
                    %values
                    diode_intervall = trigger_onsets(4,k):trigger_onsets(5,k);
                    D.data(279,diode_intervall) = 32;
                end

                save([localPath, filesep, 'spmeeg.mat'], 'D');
                clear D
                load([localPath, filesep, 'spmeeg.mat']);
               
                trigger_onsets_stimuli = trigger_onsets;

            %align resting state onset/offset triggers
            elseif strcmp(is.MEGruns{iSj}{iR},'iniRest') ...
                || strcmp(is.MEGruns{iSj}{iR},'midRest')
                
                 load([localPath, filesep, 'spmeeg.mat']);

                 %find trigger onsets
                 trigger_start_onset = [];
                 trigger_end_onset = [];
                 trigger_onsets = [];
                 
                 if strcmp(is.MEGruns{iSj}{iR},'iniRest')
                    trigger_numbers = is.rest_boundaries(1:2);
                 else
                     trigger_numbers = is.rest_boundaries(3:4);
                 end
                 
                 trigger_start_onset = find(D.data(279,:) == trigger_numbers(1));
                 trigger_end_onset = find(D.data(279,:) == trigger_numbers(2));

                 trigger_onsets(1,1) = trigger_start_onset(1);
                 trigger_onsets(2,1) = trigger_start_onset(end);
                 trigger_onsets(1,2) = trigger_end_onset(1);
                 trigger_onsets(2,2) = trigger_end_onset(end);
                 
                 trigger_onsets(3,:) = trigger_onsets(2,:) - trigger_onsets(1,:);
                 
                 %find photodiode change points 
                 diode_onsets_tmp = [];
                 diode_onsets = [];

                 diode_onsets_tmp = find(ischange(D.data(278,:),'Threshold',1));
                 
                 %find diode change point closest in time to stimulus triggers
                 for k = 1:size(diode_onsets_tmp,2)
                     current_idx = [];
                     %only if reasonable, 300 sample points after trigger
                     %onset
                     if numel(find((diode_onsets_tmp(k) - trigger_onsets(1,:) < 300) ...
                                 & (diode_onsets_tmp(k) - trigger_onsets(1,:)  >= 0))) > 0 
                             
                         current_idx = find((diode_onsets_tmp(k) - trigger_onsets(1,:) < 300) ...
                                          & (diode_onsets_tmp(k) - trigger_onsets(1,:) >= 0));

                         trigger_onsets(4,current_idx) = diode_onsets_tmp(k) - 1;
                         
                         
                     end
                    
                 end
                 
                 %replace end point in case end point of break was not recorded
                 if trigger_onsets(4,2) == 0
                    trigger_onsets(4,2) = trigger_onsets(1,2) + (trigger_onsets(4,1) - trigger_onsets(1,1));
                 end
                 
                trigger_onsets(5,:) = trigger_onsets(4,:) + trigger_onsets(3,:);

                %save old triggers                
                old_triggers = [];
                old_triggers = D.data(279,:);
                save([localPath, filesep, 'old_triggers.mat'], 'old_triggers')
                
                %delete old triggers and replace
                old_trigger_idx = [];
                old_trigger_idx = find(D.data(279,:) == trigger_numbers(1) | D.data(279,:) == trigger_numbers(2));
                
                D.data(279,old_trigger_idx) = 0;
                
                for k = 1:size(trigger_onsets,2)
                    diode_intervall = [];
                    %specify diode intervall and replace stimulus trigger
                    %values
                    diode_intervall = trigger_onsets(4,k):trigger_onsets(5,k);
                    D.data(279,diode_intervall) = trigger_numbers(k);
                end
                   
                save([localPath, filesep, 'spmeeg.mat'], 'D');
                clear D
                load([localPath, filesep, 'spmeeg.mat']);
                
                if strcmp(is.MEGruns{iSj}{iR},'iniRest')
                    time_points_pauses = [24000 48000];
                else
                    time_points_pauses = [1 12000];
                end
                
                trigger_onsets_pause = trigger_onsets;
            end
            

            %% Crop unused data
            S = struct;
            S.D = fullfile(localPath,'spmeeg.mat');
            S.prefix='';            
            event = ft_read_event(ds_folder);
            
            %replace trigger timepoints with diode timepoints
            if strcmp(is.MEGruns{iSj}{iR},'loc1')
           
                trigger_onsets_to_del = trigger_onsets_stimuli;
                for k = 1:size(event,2)
                    if event(k).value == 32
                       event(k).sample = trigger_onsets_to_del(4,1);
                       
                       trigger_onsets_to_del(:,1) = [];  
                    end
                    
                end
                
            %also with resting state phase onsets/offsets
            elseif strcmp(is.MEGruns{iSj}{iR},'iniRest') ...
                || strcmp(is.MEGruns{iSj}{iR},'midRest')    
                trigger_onsets_to_del = trigger_onsets_pause;

                for k = 1:size(event,2)
                    if event(k).value == trigger_numbers(1)
                       event(k).sample = trigger_onsets_to_del(4,1);
                       
                       trigger_onsets_to_del(:,1) = [];  
                       
                    elseif event(k).value == trigger_numbers(2)
                       event(k).sample = trigger_onsets_to_del(4,1);
                       
                       trigger_onsets_to_del(:,1) = [];  
                       
                    end
                end
                
            end
            sample = [event(find(strcmp('frontpanel trigger', {event.type}))).sample]';
            S.timewin =  [0,round(sample(end)/is.sampleRate*1000)];

            D=spm_eeg_crop_YL(S);
            D.save()
            
            %% OPT Preprocessing Pipeline
            
            %% Phase 1 - Filter
            opt=[];
            
            spm_files{1}=fullfile(localPath,'spmeeg.mat'); 
            structural_files{1}=[]; % leave empty if no .nii structural file available
            
            opt.spm_files=spm_files;
            opt.datatype='ctf';
            
            % HIGHPASS
            opt.highpass.cutoff = is.highpass;
            opt.dirname=fullfile(localPath,['highpass_',num2str(0.5)]); 
            opt.highpass.do=1;
            
            % Notch filter settings
            opt.mains.do=1;
            
            % DOWNSAMPLING
            opt.downsample.do=0;
            
            % IDENTIFYING BAD SEGMENTS 
            opt.bad_segments.do=0;
            
            % Set to 0 for now
            opt.africa.do=0;
            opt.epoch.do=0;
            opt.outliers.do=0;
            opt.coreg.do=0;
            
            %%%%%%%%%%%%%%%%%%%%%
            opt = osl_run_opt(opt);
            
            %% Phase 2 - downsample + bad segment            
            opt2=[];
            
            opt2.spm_files = opt.results.spm_files;
            opt2.datatype='ctf';

            % optional inputs
            opt2.dirname=opt.dirname; % directory opt settings and results will be stored    

            % DOWNSAMPLING
            opt2.downsample.do=1;
            opt2.downsample.freq=100;
            
            % IDENTIFYING BAD SEGMENTS 
            opt2.bad_segments.do=1;
            
            % Set to 0 for now
            opt2.africa.do=0;
            opt2.epoch.do=0;
            opt2.outliers.do=0;
            opt2.coreg.do=0;
            opt2.highpass.do=0;
            opt2.mains.do=0;
            
            %%%%%%%%%%%%%%%%%%%%%
            opt2 = osl_run_opt(opt2);         
  
            %% phase 3 - ICA 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            rng(101) %set seed for reproducibility
            
            opt3 = [];
            % required inputs
            opt3.spm_files = opt2.results.spm_files;
            opt3.datatype='ctf';

            % optional inputs
            opt3.dirname=opt2.dirname; % directory opt settings and results will be stored

            opt3.africa.do=1;
            opt3.africa.ident.artefact_chans = {'EOG1','EOG2','EOG3'};

            opt3.africa.precompute_topos = false;
            opt3.africa.ident.mains_kurt_thresh = 0.5;
            opt3.africa.ident.do_kurt = true;
            opt3.africa.ident.do_cardiac = true; 
            opt3.africa.ident.do_plots = true;
            opt3.africa.ident.do_mains = false;
             
            opt3.africa.ident.do_ident = 'auto';

            % turn the remaining options off
            opt3.downsample.do=0;
            opt3.highpass.do=0;
            opt3.bad_segments.do=0;
            opt3.epoch.do=0;
            opt3.outliers.do=0;
            opt3.coreg.do=0;
            
            %%%%%%%%%%%%%%%%%%%%%
            opt3 = osl_run_opt(opt3);
                           
            %% Phase 4 - Epoch + outliter + coreg
            opt4 = [];
            % required inputs
            opt4.spm_files = opt3.results.spm_files;

            opt4.datatype='ctf';
            
            % optional inputs
            opt4.dirname=opt3.dirname; % directory opt settings and results will be stored
            
            %% Epoching settings - different for localizer, rest and main task 
  
            %% functional localizer
            if strcmp(is.MEGruns{iSj}{iR},'loc1')
                
                opt4.epoch.do=1;
                opt4.epoch.time_range = is.epoch_win_s; 

                opt4.epoch.trialdef(1).conditionlabel = 'Stimulus Onset';
                opt4.epoch.trialdef(1).eventtype = 'frontpanel trigger';
                opt4.epoch.trialdef(1).eventvalue = is.sfl_triggers(1); 
 
                opt4.bad_segments.do=0;
                opt4.outliers.do=1;
                
            %% rest periods, no epoching
                
            elseif strcmp(is.MEGruns{iSj}{iR},'iniRest') ...
                    || strcmp(is.MEGruns{iSj}{iR},'midRest')  
                opt4.epoch.do=0;
                opt4.outliers.do=0;
                opt4.bad_segments.do=1;
                
            end

            %% coreg for subsequent source analysis
            opt4.coreg.do=1;
            opt4.coreg.use_rhino=0;
            opt4.coreg.useheadshape=0;
            % turn the remaining options off
            opt4.downsample.do=0;
            opt4.africa.todo.ica=0;
            opt4.africa.todo.ident=0;
            opt4.africa.todo.remove=0;
            
            %%%%%%%%%%%%%%%%%%%%%
            opt4=osl_run_opt(opt4);
            
            %% Display Results
            opt4 = osl_load_opt(opt4.dirname);
            close all;
            fclose('all'); 
            
        end
        
        %% Assign stimulus labels to onset times
        % assign stimuli/trials to trigger events
        if strcmp(is.MEGruns{iSj}{iR},'loc1')
            
            clear D
            clear temp_stim_PRE
            clear temp_stim_POST
            
            load([opt4.dirname, filesep, 'ReABdffspmeeg']);
            
            
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
            
            
            %find specific trials for this run
            if sum(is.fnSID{iSj} == '512') == 3
                %only last 162 trials in run 2 were recorded due to
                %technical error (42 trials less)
                %change starting and end points of runs
                start_points = [1, 205, 205+162, 367+204, 571+204];
                end_points = [204, 204+162, 366+204, 570+204, 774+204];
                
                %cut down stimulus labels
                temp_stim_PRE_names = [temp_stim_PRE_names(1:204); temp_stim_PRE_names(247:end)];
                all_run_idx = find(strcmp(is.MEGruns{iSj},'loc1'));
                
            else
                trial_length = size(temp_stim_PRE_names,1) / 5;
                start_points = 1:trial_length:size(temp_stim_PRE_names,1);
                end_points = trial_length:trial_length:size(temp_stim_PRE_names,1);
                all_run_idx = find(strcmp(is.MEGruns{iSj},'loc1'));
            end
            
            temp_stim_PRE_names = temp_stim_PRE_names(start_points(find(all_run_idx == iR)):end_points(find(all_run_idx == iR)));
            
            %replace stimulus names in trial structure (MEG data)
            for i = 1:size(D.trials,2)
                D.trials(i).label = temp_stim_PRE_names(i);
            end
            
            save([opt4.dirname, filesep, 'ReABdffspmeeg'], 'D');
        end
        
    end %end run loop
    
end %end subject loop


            
            
            
