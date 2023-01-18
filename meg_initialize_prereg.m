%% Initialize paths, parameters etc

% sets paths
% populates the 'is' parameter structure

is = struct;
is.prereg_path = pwd;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% PATHS 
%--------------------------------------------------------------------------
% toolboxes and analysis scripts 
cd('/Users/lennart/Documents/MATLAB/osl/osl-core');
OSLDIR = getenv('OSLDIR');
addpath(OSLDIR);
osl_startup(OSLDIR);
is.spm_path = '/Users/lennart/Documents/MATLAB/osl/spm12';

cd(is.prereg_path)
is.rawMEG =  '/Users/lennart/Desktop/Data_Transfer_HHU_Mac/PostDoc/Projects/Markov/task/data/data/meg_data/meg_study/neural/OneDrive - University College London/raw';                     % raw MEG from scanner
is.OPTPath =  'preproc/';            % savedir for preprocessed data
is.OutputDir = 'neural/';     % savedir for TDLM
is.behav_path = 'behav';                     % behavioral data files

if ~isfolder([is.OutputDir 'trained_decoders_compounds'])
    mkdir([is.OutputDir 'trained_decoders_compounds']);
    mkdir([is.OutputDir 'decoded_rest_compounds']);
    mkdir([is.OutputDir 'resting_results']);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% PRE-PROCESSING PARAMS
%--------------------------------------------------------------------------
is.numChan = 273;           % 2 of 275 MEG channels down
is.sampleRate = 1200;       % sample rate of raw MEG data 
is.smoothFact = 12;         % to get 100Hz final sample rate
is.msPerSample = 1000 * is.smoothFact / is.sampleRate; % 10 ms, ms per sample of preprocessed MEG data (100Hz)
is.highpass = 0.5;          % high pass filter Hz

% functional localizer and main task epoching
%Relevant triggers
is.sfl_triggers = [32 122];

%epoching windows
is.epoch_win_s = [-0.2 .8]; %[-0.7 1] before;  % epoch win w.r.t. pic onset trigger
is.epoch_win_bins = 1000*is.epoch_win_s/is.msPerSample;
is.picOnsetSlice = -(is.epoch_win_bins(1)) + 1;         % time bin of pic onset in epoched MEG data

% rest session start and end triggers
is.rest_boundaries = [65 122 222 237];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% DECODING AND TDLM ANALYSIS PARAMS
%--------------------------------------------------------------------------
% Decoder training params
is.ENL1 = 0.008;                % L1 parameter = 0.008
is.nullslice = 200; %specify null slice time point, 200 ms before stimulus onset

% decoder train time
is.tss = 18; % bins post-pic onset for decoder training

is.lgncy = 1:60;                % lags in units of samples

is.Ltimes = is.lgncy*is.msPerSample/1000;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SUBJECTS
%--------------------------------------------------------------------------
meg_SpecifySubjects_prereg
