
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file)

function [config, mergeindex] = dtx_setparams_patients_lgi1(config)

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-Patients_LGI1\';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\Patients_LGI1\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis);
mergeindex = {[1, 2]};%, [3, 4, 5]};

%% Patient 1
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%LFP.NAME : toujours commencer par right
%VERIF ALL ELECTRODES. PAS RAJOUTER CELLES EN PLUS
config{1}.os                        = os;
config{1}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{1}.prefix                    = 'pat_LGI1_008-EEG_129-';
config{1}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'Crise_Start','Crise_End'; 'Crise_End','SlowWave'};   % start and end Muse marker
%config{1}.foldername                = 'pat_LGI1_008';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{1}.datasavedir               = datasavedir;         % where to write data
config{1}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
config{1}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
config{1}.labels.emg                = {'EMG1+','EMG2+'};
config{1}.directorylist{1}          = {'EEG_129'}; %dir = eeg file with all the electrodess

%config{1}.preproc_eog %TO DO

config{1}.align.name                = {'SlowWave_R'};%{'SlowWave_R','SlowWave_R'};
config{1}.align.channel             = {'C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.align.flip                = {'no'};
config{1}.align.abs                 = {'no'};
config{1}.align.method              = {'min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{1}.align.filter              = {'lp'};
config{1}.align.freq                = {5};          % lowpass filter freq to smooth peak detection (Hz)
config{1}.align.hilbert             = {'no'};
config{1}.align.thresh              = [0];
config{1}.align.toiplot{1}          = [-1,  1];     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toiactive{1}        = [-0.5, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toibaseline{1}      = [-1, -0.5];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.reref       = 'yes';
config{1}.align.rerefmethod = 'avg';
config{1}.align.refchannel  = config{1}.labels.macro';
config{1}.align.notch       = 'yes';


config{1}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'}; %alwa
config{1}.LFP.emg                   = {'EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 0;
config{1}.LFP.resamplefs            = 256;
config{1}.LFP.baseline              = 'no';
config{1}.LFP.baselinewindow{1}     = [-2, -1];
config{1}.LFP.slidestep             = 0.01;
config{1}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{1}.LFP.reref                 = 'yes';
config{1}.LFP.rerefmethod           = 'avg';
config{1}.LFP.refchannel            = config{1}.labels.macro';
config{1}.LFP.bsfilter              = 'yes';
config{1}.LFP.bsfreq                = [49 51];
config{1}.LFP.hasemg                = 'yes';
config{1}.EMG.hpfilter              = 'yes';
config{1}.EMG.hpfreq                = 10;
config{1}.EMG.bsfilter              = 'yes';
config{1}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{1}.epoch.toi{1}              = [-10, 10];
config{1}.epoch.toi{2}              = [-10, 10];
config{1}.epoch.toi{3}              = [1, -2];
config{1}.epoch.pad{1}              = 10;
config{1}.epoch.pad{2}              = 10;
config{1}.epoch.pad{3}              = 0.5;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{2}.os                        = os;
config{2}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{2}.prefix                    = 'pat_LGI1_008-EEG_131-';
config{2}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{2}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{2}.datasavedir               = datasavedir;         % where to write data
config{2}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
config{2}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{2}.labels.emg                = {'EMG1+','EMG2+'};
config{2}.directorylist{1}          = {'EEG_131'}; %dir = eeg file with all the electrodess

%config{2}.preproc_eog %TO DO

config{2}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{2}.align.channel             = {'C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.align.flip                = {'no','no'};
config{2}.align.abs                 = {'no','no'};
config{2}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{2}.align.filter              = {'lp','lp'};
config{2}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{2}.align.hilbert             = {'no','no'};
config{2}.align.thresh              = [0, 0];
config{2}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.reref               = 'yes';
config{2}.align.rerefmethod         = 'avg';
config{2}.align.refchannel          = config{2}.labels.macro';
config{2}.align.notch               = 'yes';


config{2}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{2}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{2}.LFP.hpfilter              = 'no';
config{2}.LFP.hpfreq                = 0;
config{2}.LFP.resamplefs            = 256;
config{2}.LFP.baseline              = 'no';
config{2}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{2}.LFP.slidestep             = 0.01;
config{2}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{2}.LFP.reref                 = 'yes';
config{2}.LFP.rerefmethod           = 'avg';
config{2}.LFP.refchannel            = config{2}.labels.macro';
config{2}.LFP.bsfilter              = 'yes';
config{2}.LFP.bsfreq                = [49 51];
config{2}.LFP.hasemg                = 'yes';
config{2}.EMG.hpfilter              = 'yes';
config{2}.EMG.hpfreq                = 10;
config{2}.EMG.bsfilter              = 'yes';
config{2}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{2}.epoch.toi{1}              = [-10, 10];
config{2}.epoch.toi{2}              = [-10, 10];
config{2}.epoch.toi{3}              = [1, -2];
config{2}.epoch.pad{1}              = 10;
config{2}.epoch.pad{2}              = 10;
config{2}.epoch.pad{3}              = 0.5;
end


