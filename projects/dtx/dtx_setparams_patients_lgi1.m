
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
mergeindex = {[1, 2], [3, 4, 5], [6, 7], [9, 10, 11, 12, 13, 14]};

%% Patient 1
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%ADD EOG CHANNEL IF ANY
%LFP.NAME : toujours commencer par right
%VERIF ALL ELECTRODES. PAS RAJOUTER CELLES EN PLUS
%DEFINIR CONGIG_COMMON.ALIGN .LFP .EMG .EPOCH
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
config{2}.epoch.toi{3}              = [-10, 10];
config{2}.epoch.pad{1}              = 10;
config{2}.epoch.pad{2}              = 10;
config{2}.epoch.pad{3}              = 10;

%% Patient 3
config{3}.os                        = os;
config{3}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{3}.prefix                    = 'pat_LGI1_001-EEG_1-';
config{3}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{3}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{3}.datasavedir               = datasavedir;         % where to write data
config{3}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{3}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{3}.labels.emg                = {'EMG1+','EMG2+'};
config{3}.directorylist{1}          = {'EEG_1'}; %dir = eeg file with all the electrodess

%config{3}.preproc_eog %TO DO

config{3}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{3}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{3}.align.flip                = {'no','no'};
config{3}.align.abs                 = {'no','no'};
config{3}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{3}.align.filter              = {'lp','lp'};
config{3}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{3}.align.hilbert             = {'no','no'};
config{3}.align.thresh              = [0, 0];
config{3}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{3}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{3}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{3}.align.reref               = 'yes';
config{3}.align.rerefmethod         = 'avg';
config{3}.align.refchannel          = config{3}.labels.macro';
config{3}.align.notch               = 'yes';


config{3}.LFP.name                  = {'SlowWave_R','SlowWave_L'};
config{3}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{3}.LFP.hpfilter              = 'no';
config{3}.LFP.hpfreq                = 0;
config{3}.LFP.resamplefs            = 256;
config{3}.LFP.baseline              = 'no';
config{3}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{3}.LFP.slidestep             = 0.01;
config{3}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{3}.LFP.reref                 = 'yes';
config{3}.LFP.rerefmethod           = 'avg';
config{3}.LFP.refchannel            = config{3}.labels.macro';
config{3}.LFP.bsfilter              = 'yes';
config{3}.LFP.bsfreq                = [49 51];
config{3}.LFP.hasemg                = 'yes';
config{3}.EMG.hpfilter              = 'yes';
config{3}.EMG.hpfreq                = 10;
config{3}.EMG.bsfilter              = 'yes';
config{3}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{3}.epoch.toi{1}              = [-10, 10];
config{3}.epoch.toi{2}              = [-10, 10];
config{3}.epoch.toi{3}              = [-10, 10];
config{3}.epoch.pad{1}              = 10;
config{3}.epoch.pad{2}              = 10;
config{3}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{4}.os                        = os;
config{4}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{4}.prefix                    = 'pat_LGI1_001-EEG_14-';
config{4}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{4}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{4}.datasavedir               = datasavedir;         % where to write data
config{4}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{4}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{4}.labels.emg                = {'EMG1+','EMG2+'};
config{4}.directorylist{1}          = {'EEG_14'}; %dir = eeg file with all the electrodess

%config{4}.preproc_eog %TO DO

config{4}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{4}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{4}.align.flip                = {'no','no'};
config{4}.align.abs                 = {'no','no'};
config{4}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{4}.align.filter              = {'lp','lp'};
config{4}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{4}.align.hilbert             = {'no','no'};
config{4}.align.thresh              = [0, 0];
config{4}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{4}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{4}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{4}.align.reref               = 'yes';
config{4}.align.rerefmethod         = 'avg';
config{4}.align.refchannel          = config{4}.labels.macro';
config{4}.align.notch               = 'yes';


config{4}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{4}.LFP.emg                   = {'EMG1+','EMG2+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{4}.LFP.hpfilter              = 'no';
config{4}.LFP.hpfreq                = 0;
config{4}.LFP.resamplefs            = 256;
config{4}.LFP.baseline              = 'no';
config{4}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{4}.LFP.slidestep             = 0.01;
config{4}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{4}.LFP.reref                 = 'yes';
config{4}.LFP.rerefmethod           = 'avg';
config{4}.LFP.refchannel            = config{4}.labels.macro';
config{4}.LFP.bsfilter              = 'yes';
config{4}.LFP.bsfreq                = [49 51];
config{4}.LFP.hasemg                = 'yes';
config{4}.EMG.hpfilter              = 'yes';
config{4}.EMG.hpfreq                = 10;
config{4}.EMG.bsfilter              = 'yes';
config{4}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{4}.epoch.toi{1}              = [-10, 10];
config{4}.epoch.toi{2}              = [-10, 10];
config{4}.epoch.toi{3}              = [-10, 10];
config{4}.epoch.pad{1}              = 10;
config{4}.epoch.pad{2}              = 10;
config{4}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{5}.os                        = os;
config{5}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{5}.prefix                    = 'pat_LGI1_001-EEG_10-';
config{5}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{5}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{5}.datasavedir               = datasavedir;         % where to write data
config{5}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{5}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{5}.labels.emg                = {'EMG1+','EMG2+'};
config{5}.directorylist{1}          = {'EEG_10'}; %dir = eeg file with all the electrodess

%config{5}.preproc_eog %TO DO

config{5}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{5}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{5}.align.flip                = {'no','no'};
config{5}.align.abs                 = {'no','no'};
config{5}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{5}.align.filter              = {'lp','lp'};
config{5}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{5}.align.hilbert             = {'no','no'};
config{5}.align.thresh              = [0, 0];
config{5}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{5}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{5}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{5}.align.reref               = 'yes';
config{5}.align.rerefmethod         = 'avg';
config{5}.align.refchannel          = config{5}.labels.macro';
config{5}.align.notch               = 'yes';


config{5}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{5}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{5}.LFP.hpfilter              = 'no';
config{5}.LFP.hpfreq                = 0;
config{5}.LFP.resamplefs            = 256;
config{5}.LFP.baseline              = 'no';
config{5}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{5}.LFP.slidestep             = 0.01;
config{5}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{5}.LFP.reref                 = 'yes';
config{5}.LFP.rerefmethod           = 'avg';
config{5}.LFP.refchannel            = config{5}.labels.macro';
config{5}.LFP.bsfilter              = 'yes';
config{5}.LFP.bsfreq                = [49 51];
config{5}.LFP.hasemg                = 'yes';
config{5}.EMG.hpfilter              = 'yes';
config{5}.EMG.hpfreq                = 10;
config{5}.EMG.bsfilter              = 'yes';
config{5}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{5}.epoch.toi{1}              = [-10, 10];
config{5}.epoch.toi{2}              = [-10, 10];
config{5}.epoch.toi{3}              = [-10, 10];
config{5}.epoch.pad{1}              = 10;
config{5}.epoch.pad{2}              = 10;
config{5}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{6}.os                        = os;
config{6}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{6}.prefix                    = 'pat_LGI1_007-EEG_113-';
config{6}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{6}.datasavedir               = datasavedir;         % where to write data
config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
config{6}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{6}.labels.emg                = {'EMG1+','EMG2+'};
config{6}.directorylist{1}          = {'EEG_113'}; %dir = eeg file with all the electrodess

%config{6}.preproc_eog %TO DO

config{6}.align.name                = {'SlowWave_R'};%,'SlowWave_L'};
config{6}.align.channel             = {'C4'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{6}.align.flip                = {'no','no'};
config{6}.align.abs                 = {'no','no'};
config{6}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{6}.align.filter              = {'lp','lp'};
config{6}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{6}.align.hilbert             = {'no','no'};
config{6}.align.thresh              = [0, 0];
config{6}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{6}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{6}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{6}.align.reref               = 'yes';
config{6}.align.rerefmethod         = 'avg';
config{6}.align.refchannel          = config{6}.labels.macro';
config{6}.align.notch               = 'yes';


config{6}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'};%,'SlowWave_L'};
config{6}.LFP.emg                   = {'EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{6}.LFP.hpfilter              = 'no';
config{6}.LFP.hpfreq                = 0;
config{6}.LFP.resamplefs            = 256;
config{6}.LFP.baseline              = 'no';
config{6}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{6}.LFP.slidestep             = 0.01;
config{6}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{6}.LFP.reref                 = 'yes';
config{6}.LFP.rerefmethod           = 'avg';
config{6}.LFP.refchannel            = config{6}.labels.macro';
config{6}.LFP.bsfilter              = 'yes';
config{6}.LFP.bsfreq                = [49 51];
config{6}.LFP.hasemg                = 'yes';
config{6}.EMG.hpfilter              = 'yes';
config{6}.EMG.hpfreq                = 10;
config{6}.EMG.bsfilter              = 'yes';
config{6}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{6}.epoch.toi{1}              = [-10, 10];
config{6}.epoch.toi{2}              = [-10, 10];
config{6}.epoch.toi{3}              = [-10, 10];
config{6}.epoch.pad{1}              = 10;
config{6}.epoch.pad{2}              = 10;
config{6}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{7}.os                        = os;
config{7}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{7}.prefix                    = 'pat_LGI1_007-EEG_123-';
config{7}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{7}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{7}.datasavedir               = datasavedir;         % where to write data
config{7}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
config{7}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{7}.labels.emg                = {'EMG1+','EMG2+'};
config{7}.directorylist{1}          = {'EEG_123'}; %dir = eeg file with all the electrodess

%config{7}.preproc_eog %TO DO

config{7}.align.name                = {'SlowWave_R'};%,'SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{7}.align.channel             = {'C4'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{7}.align.flip                = {'no','no'};
config{7}.align.abs                 = {'no','no'};
config{7}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{7}.align.filter              = {'lp','lp'};
config{7}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{7}.align.hilbert             = {'no','no'};
config{7}.align.thresh              = [0, 0];
config{7}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{7}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{7}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{7}.align.reref               = 'yes';
config{7}.align.rerefmethod         = 'avg';
config{7}.align.refchannel          = config{7}.labels.macro';
config{7}.align.notch               = 'yes';


config{7}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'};%,'SlowWave_L'};
config{7}.LFP.emg                   = {'EMG1+'};%,'no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{7}.LFP.hpfilter              = 'no';
config{7}.LFP.hpfreq                = 0;
config{7}.LFP.resamplefs            = 256;
config{7}.LFP.baseline              = 'no';
config{7}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{7}.LFP.slidestep             = 0.01;
config{7}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{7}.LFP.reref                 = 'yes';
config{7}.LFP.rerefmethod           = 'avg';
config{7}.LFP.refchannel            = config{7}.labels.macro';
config{7}.LFP.bsfilter              = 'yes';
config{7}.LFP.bsfreq                = [49 51];
config{7}.LFP.hasemg                = 'yes';
config{7}.EMG.hpfilter              = 'yes';
config{7}.EMG.hpfreq                = 10;
config{7}.EMG.bsfilter              = 'yes';
config{7}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{7}.epoch.toi{1}              = [-10, 10];
config{7}.epoch.toi{2}              = [-10, 10];
config{7}.epoch.toi{3}              = [-10, 10];
config{7}.epoch.pad{1}              = 10;
config{7}.epoch.pad{2}              = 10;
config{7}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{8}.os                        = os;
config{8}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{8}.prefix                    = 'pat_LGI1_009-EEG_156-';
config{8}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{8}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_009');
config{8}.datasavedir               = datasavedir;         % where to write data
config{8}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_009');       % where to print images
config{8}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{8}.labels.emg                = {'EMG1+','EMG2+'};
config{8}.directorylist{1}          = {'EEG_156'}; %dir = eeg file with all the electrodess

%config{8}.preproc_eog %TO DO

config{8}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{8}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{8}.align.flip                = {'no','no'};
config{8}.align.abs                 = {'no','no'};
config{8}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{8}.align.filter              = {'lp','lp'};
config{8}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{8}.align.hilbert             = {'no','no'};
config{8}.align.thresh              = [0, 0];
config{8}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{8}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{8}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{8}.align.reref               = 'yes';
config{8}.align.rerefmethod         = 'avg';
config{8}.align.refchannel          = config{8}.labels.macro';
config{8}.align.notch               = 'yes';


config{8}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{8}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{8}.LFP.hpfilter              = 'no';
config{8}.LFP.hpfreq                = 0;
config{8}.LFP.resamplefs            = 256;
config{8}.LFP.baseline              = 'no';
config{8}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{8}.LFP.slidestep             = 0.01;
config{8}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{8}.LFP.reref                 = 'yes';
config{8}.LFP.rerefmethod           = 'avg';
config{8}.LFP.refchannel            = config{8}.labels.macro';
config{8}.LFP.bsfilter              = 'yes';
config{8}.LFP.bsfreq                = [49 51];
config{8}.LFP.hasemg                = 'yes';
config{8}.EMG.hpfilter              = 'yes';
config{8}.EMG.hpfreq                = 10;
config{8}.EMG.bsfilter              = 'yes';
config{8}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{8}.epoch.toi{1}              = [-10, 10];
config{8}.epoch.toi{2}              = [-10, 10];
config{8}.epoch.toi{3}              = [-10, 10];
config{8}.epoch.pad{1}              = 10;
config{8}.epoch.pad{2}              = 10;
config{8}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{9}.os                        = os;
config{9}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{9}.prefix                    = 'pat_LGI1_010-EEG_172-';
config{9}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{9}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{9}.datasavedir               = datasavedir;         % where to write data
config{9}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{9}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{9}.labels.emg                = {'EMG1+','EMG2+'};
config{9}.directorylist{1}          = {'EEG_172'}; %dir = eeg file with all the electrodess

%config{9}.preproc_eog %TO DO

config{9}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{9}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{9}.align.flip                = {'no','no'};
config{9}.align.abs                 = {'no','no'};
config{9}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{9}.align.filter              = {'lp','lp'};
config{9}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{9}.align.hilbert             = {'no','no'};
config{9}.align.thresh              = [0, 0];
config{9}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{9}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{9}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{9}.align.reref               = 'yes';
config{9}.align.rerefmethod         = 'avg';
config{9}.align.refchannel          = config{9}.labels.macro';
config{9}.align.notch               = 'yes';


config{9}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{9}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{9}.LFP.hpfilter              = 'no';
config{9}.LFP.hpfreq                = 0;
config{9}.LFP.resamplefs            = 256;
config{9}.LFP.baseline              = 'no';
config{9}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{9}.LFP.slidestep             = 0.01;
config{9}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{9}.LFP.reref                 = 'yes';
config{9}.LFP.rerefmethod           = 'avg';
config{9}.LFP.refchannel            = config{9}.labels.macro';
config{9}.LFP.bsfilter              = 'yes';
config{9}.LFP.bsfreq                = [49 51];
config{9}.LFP.hasemg                = 'yes';
config{9}.EMG.hpfilter              = 'yes';
config{9}.EMG.hpfreq                = 10;
config{9}.EMG.bsfilter              = 'yes';
config{9}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{9}.epoch.toi{1}              = [-10, 10];
config{9}.epoch.toi{2}              = [-10, 10];
config{9}.epoch.toi{3}              = [-10, 10];
config{9}.epoch.pad{1}              = 10;
config{9}.epoch.pad{2}              = 10;
config{9}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{10}.os                        = os;
config{10}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{10}.prefix                    = 'pat_LGI1_010-EEG_203-';
config{10}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{10}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{10}.datasavedir               = datasavedir;         % where to write data
config{10}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{10}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{10}.labels.emg                = {'EMG1+','EMG2+'};
config{10}.directorylist{1}          = {'EEG_203'}; %dir = eeg file with all the electrodess

%config{10}.preproc_eog %TO DO

config{10}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{10}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{10}.align.flip                = {'no','no'};
config{10}.align.abs                 = {'no','no'};
config{10}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{10}.align.filter              = {'lp','lp'};
config{10}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{10}.align.hilbert             = {'no','no'};
config{10}.align.thresh              = [0, 0];
config{10}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{10}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{10}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{10}.align.reref               = 'yes';
config{10}.align.rerefmethod         = 'avg';
config{10}.align.refchannel          = config{10}.labels.macro';
config{10}.align.notch               = 'yes';


config{10}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{10}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{10}.LFP.hpfilter              = 'no';
config{10}.LFP.hpfreq                = 0;
config{10}.LFP.resamplefs            = 256;
config{10}.LFP.baseline              = 'no';
config{10}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{10}.LFP.slidestep             = 0.01;
config{10}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{10}.LFP.reref                 = 'yes';
config{10}.LFP.rerefmethod           = 'avg';
config{10}.LFP.refchannel            = config{10}.labels.macro';
config{10}.LFP.bsfilter              = 'yes';
config{10}.LFP.bsfreq                = [49 51];
config{10}.LFP.hasemg                = 'yes';
config{10}.EMG.hpfilter              = 'yes';
config{10}.EMG.hpfreq                = 10;
config{10}.EMG.bsfilter              = 'yes';
config{10}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{10}.epoch.toi{1}              = [-10, 10];
config{10}.epoch.toi{2}              = [-10, 10];
config{10}.epoch.toi{3}              = [-10, 10];
config{10}.epoch.pad{1}              = 10;
config{10}.epoch.pad{2}              = 10;
config{10}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{11}.os                        = os;
config{11}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{11}.prefix                    = 'pat_LGI1_010-EEG_174-';
config{11}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{11}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{11}.datasavedir               = datasavedir;         % where to write data
config{11}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{11}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{11}.labels.emg                = {'EMG1+','EMG2+'};
config{11}.directorylist{1}          = {'EEG_174'}; %dir = eeg file with all the electrodess

%config{11}.preproc_eog %TO DO

config{11}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{11}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{11}.align.flip                = {'no','no'};
config{11}.align.abs                 = {'no','no'};
config{11}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{11}.align.filter              = {'lp','lp'};
config{11}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{11}.align.hilbert             = {'no','no'};
config{11}.align.thresh              = [0, 0];
config{11}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{11}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{11}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{11}.align.reref               = 'yes';
config{11}.align.rerefmethod         = 'avg';
config{11}.align.refchannel          = config{11}.labels.macro';
config{11}.align.notch               = 'yes';


config{11}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{11}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{11}.LFP.hpfilter              = 'no';
config{11}.LFP.hpfreq                = 0;
config{11}.LFP.resamplefs            = 256;
config{11}.LFP.baseline              = 'no';
config{11}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{11}.LFP.slidestep             = 0.01;
config{11}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{11}.LFP.reref                 = 'yes';
config{11}.LFP.rerefmethod           = 'avg';
config{11}.LFP.refchannel            = config{11}.labels.macro';
config{11}.LFP.bsfilter              = 'yes';
config{11}.LFP.bsfreq                = [49 51];
config{11}.LFP.hasemg                = 'yes';
config{11}.EMG.hpfilter              = 'yes';
config{11}.EMG.hpfreq                = 10;
config{11}.EMG.bsfilter              = 'yes';
config{11}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{11}.epoch.toi{1}              = [-10, 10];
config{11}.epoch.toi{2}              = [-10, 10];
config{11}.epoch.toi{3}              = [-10, 10];
config{11}.epoch.pad{1}              = 10;
config{11}.epoch.pad{2}              = 10;
config{11}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{12}.os                        = os;
config{12}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{12}.prefix                    = 'pat_LGI1_010-EEG_177-';
config{12}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{12}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{12}.datasavedir               = datasavedir;         % where to write data
config{12}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{12}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{12}.labels.emg                = {'EMG1+','EMG2+'};
config{12}.directorylist{1}          = {'EEG_177'}; %dir = eeg file with all the electrodess

%config{12}.preproc_eog %TO DO

config{12}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{12}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{12}.align.flip                = {'no','no'};
config{12}.align.abs                 = {'no','no'};
config{12}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{12}.align.filter              = {'lp','lp'};
config{12}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{12}.align.hilbert             = {'no','no'};
config{12}.align.thresh              = [0, 0];
config{12}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{12}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{12}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{12}.align.reref               = 'yes';
config{12}.align.rerefmethod         = 'avg';
config{12}.align.refchannel          = config{12}.labels.macro';
config{12}.align.notch               = 'yes';


config{12}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{12}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{12}.LFP.hpfilter              = 'no';
config{12}.LFP.hpfreq                = 0;
config{12}.LFP.resamplefs            = 256;
config{12}.LFP.baseline              = 'no';
config{12}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{12}.LFP.slidestep             = 0.01;
config{12}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{12}.LFP.reref                 = 'yes';
config{12}.LFP.rerefmethod           = 'avg';
config{12}.LFP.refchannel            = config{12}.labels.macro';
config{12}.LFP.bsfilter              = 'yes';
config{12}.LFP.bsfreq                = [49 51];
config{12}.LFP.hasemg                = 'yes';
config{12}.EMG.hpfilter              = 'yes';
config{12}.EMG.hpfreq                = 10;
config{12}.EMG.bsfilter              = 'yes';
config{12}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{12}.epoch.toi{1}              = [-10, 10];
config{12}.epoch.toi{2}              = [-10, 10];
config{12}.epoch.toi{3}              = [-10, 10];
config{12}.epoch.pad{1}              = 10;
config{12}.epoch.pad{2}              = 10;
config{12}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{13}.os                        = os;
config{13}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{13}.prefix                    = 'pat_LGI1_010-EEG_180-';
config{13}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{13}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{13}.datasavedir               = datasavedir;         % where to write data
config{13}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{13}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{13}.labels.emg                = {'EMG1+','EMG2+'};
config{13}.directorylist{1}          = {'EEG_180'}; %dir = eeg file with all the electrodess

%config{13}.preproc_eog %TO DO

config{13}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{13}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{13}.align.flip                = {'no','no'};
config{13}.align.abs                 = {'no','no'};
config{13}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{13}.align.filter              = {'lp','lp'};
config{13}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{13}.align.hilbert             = {'no','no'};
config{13}.align.thresh              = [0, 0];
config{13}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{13}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{13}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{13}.align.reref               = 'yes';
config{13}.align.rerefmethod         = 'avg';
config{13}.align.refchannel          = config{13}.labels.macro';
config{13}.align.notch               = 'yes';


config{13}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{13}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{13}.LFP.hpfilter              = 'no';
config{13}.LFP.hpfreq                = 0;
config{13}.LFP.resamplefs            = 256;
config{13}.LFP.baseline              = 'no';
config{13}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{13}.LFP.slidestep             = 0.01;
config{13}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{13}.LFP.reref                 = 'yes';
config{13}.LFP.rerefmethod           = 'avg';
config{13}.LFP.refchannel            = config{13}.labels.macro';
config{13}.LFP.bsfilter              = 'yes';
config{13}.LFP.bsfreq                = [49 51];
config{13}.LFP.hasemg                = 'yes';
config{13}.EMG.hpfilter              = 'yes';
config{13}.EMG.hpfreq                = 10;
config{13}.EMG.bsfilter              = 'yes';
config{13}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{13}.epoch.toi{1}              = [-10, 10];
config{13}.epoch.toi{2}              = [-10, 10];
config{13}.epoch.toi{3}              = [-10, 10];
config{13}.epoch.pad{1}              = 10;
config{13}.epoch.pad{2}              = 10;
config{13}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{14}.os                        = os;
config{14}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{14}.prefix                    = 'pat_LGI1_010-EEG_199-';
config{14}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{14}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{14}.datasavedir               = datasavedir;         % where to write data
config{14}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{14}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{14}.labels.emg                = {'EMG1+','EMG2+'};
config{14}.directorylist{1}          = {'EEG_199'}; %dir = eeg file with all the electrodess

%config{14}.preproc_eog %TO DO

config{14}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{14}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{14}.align.flip                = {'no','no'};
config{14}.align.abs                 = {'no','no'};
config{14}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{14}.align.filter              = {'lp','lp'};
config{14}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{14}.align.hilbert             = {'no','no'};
config{14}.align.thresh              = [0, 0];
config{14}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{14}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{14}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{14}.align.reref               = 'yes';
config{14}.align.rerefmethod         = 'avg';
config{14}.align.refchannel          = config{14}.labels.macro';
config{14}.align.notch               = 'yes';


config{14}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{14}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{14}.LFP.hpfilter              = 'no';
config{14}.LFP.hpfreq                = 0;
config{14}.LFP.resamplefs            = 256;
config{14}.LFP.baseline              = 'no';
config{14}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{14}.LFP.slidestep             = 0.01;
config{14}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{14}.LFP.reref                 = 'yes';
config{14}.LFP.rerefmethod           = 'avg';
config{14}.LFP.refchannel            = config{14}.labels.macro';
config{14}.LFP.bsfilter              = 'yes';
config{14}.LFP.bsfreq                = [49 51];
config{14}.LFP.hasemg                = 'yes';
config{14}.EMG.hpfilter              = 'yes';
config{14}.EMG.hpfreq                = 10;
config{14}.EMG.bsfilter              = 'yes';
config{14}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{14}.epoch.toi{1}              = [-10, 10];
config{14}.epoch.toi{2}              = [-10, 10];
config{14}.epoch.toi{3}              = [-10, 10];
config{14}.epoch.pad{1}              = 10;
config{14}.epoch.pad{2}              = 10;
config{14}.epoch.pad{3}              = 10;

%% Patient 2
%Same patient, other EEG.
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%FILTRE HIGH PASS POUR LES EMG : séparer preproc EEG et EMG et appendata
%LISTE DIRECTORY = DOSSIER. + LISTE FILE
%ADD EOG CHANNEL IF ANY
%ATTENTION same number of EMG in all different patients EEG
config{15}.os                        = os;
config{15}.name                      = {'SlowWave_R','Seizure','InterIctal'};
config{15}.prefix                    = 'pat_LGI1_012-EEG_213-';
config{15}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{15}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
config{15}.datasavedir               = datasavedir;         % where to write data
config{15}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
config{15}.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
%config{15}.labels.emg                = {'EMG1+','EMG2+'};
config{15}.directorylist{1}          = {'EEG_213'}; %dir = eeg file with all the electrodess

%config{15}.preproc_eog %TO DO

config{15}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{15}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{15}.align.flip                = {'no','no'};
config{15}.align.abs                 = {'no','no'};
config{15}.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{15}.align.filter              = {'lp','lp'};
config{15}.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
config{15}.align.hilbert             = {'no','no'};
config{15}.align.thresh              = [0, 0];
config{15}.align.toiplot             = {[-1,  1],[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{15}.align.toiactive           = {[-0.5, 0.5],[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{15}.align.toibaseline         = {[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{15}.align.reref               = 'yes';
config{15}.align.rerefmethod         = 'avg';
config{15}.align.refchannel          = config{15}.labels.macro';
config{15}.align.notch               = 'yes';


config{15}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{15}.LFP.emg                   = {'no','EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{15}.LFP.hpfilter              = 'no';
config{15}.LFP.hpfreq                = 0;
config{15}.LFP.resamplefs            = 256;
config{15}.LFP.baseline              = 'no';
config{15}.LFP.baselinewindow        = {[-2, -1],[-2, -1]};
config{15}.LFP.slidestep             = 0.01;
config{15}.LFP.electrodeToPlot       = [2 1 3 10 11];
config{15}.LFP.reref                 = 'yes';
config{15}.LFP.rerefmethod           = 'avg';
config{15}.LFP.refchannel            = config{15}.labels.macro';
config{15}.LFP.bsfilter              = 'yes';
config{15}.LFP.bsfreq                = [49 51];
config{15}.LFP.hasemg                = 'yes';
config{15}.EMG.hpfilter              = 'yes';
config{15}.EMG.hpfreq                = 10;
config{15}.EMG.bsfilter              = 'yes';
config{15}.EMG.bsfreq                = [49 51];


%1, 2 and 3 associated with config.muse.startend
config{15}.epoch.toi{1}              = [-10, 10];
config{15}.epoch.toi{2}              = [-10, 10];
config{15}.epoch.toi{3}              = [-10, 10];
config{15}.epoch.pad{1}              = 10;
config{15}.epoch.pad{2}              = 10;
config{15}.epoch.pad{3}              = 10;
end


