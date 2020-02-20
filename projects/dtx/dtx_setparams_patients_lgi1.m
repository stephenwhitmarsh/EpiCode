
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



%% Congig common for all patients

datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis);
mergeindex = {[1, 2], [3, 4, 5], [6, 7],[8], [9, 10, 11, 12, 13, 14],[15]};

configcommon.os                        = os;
configcommon.name                      = {'SlowWave_R','SlowWave_L','InterIctal'};
configcommon.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
configcommon.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};

%configcommon.preproc_eog %TO DO

configcommon.align.flip                = {'no','no'};
configcommon.align.abs                 = {'no','no'};
configcommon.align.method              = {'min','min'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp','lp'};
configcommon.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no','no'};
configcommon.align.thresh              = [0, 0];
configcommon.align.toiplot             = {[-1,  1], [-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.toiactive           = {[-0.5, 0.5], [-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline         = {[-1, -0.5], [-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.reref               = 'yes';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = configcommon.labels.macro';
configcommon.align.notch               = 'yes';

configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 0;
configcommon.LFP.resamplefs            = 256;
configcommon.LFP.baseline              = 'no';
configcommon.LFP.baselinewindow        = {[-2, -1], [-2, -1]};
configcommon.LFP.slidestep             = 0.01;
%configcommon.LFP.electrodeToPlot       = [2 1 3 10 11];
configcommon.LFP.reref                 = 'yes';
configcommon.LFP.rerefmethod           = 'avg';
configcommon.LFP.refchannel            = configcommon.labels.macro';
configcommon.LFP.bsfilter              = 'yes';
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.hasemg                = 'yes';
configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];

%1, 2 and 3 associated with config.muse.startend
configcommon.epoch.toi{1}              = [-10, 10];
configcommon.epoch.toi{2}              = [-10, 10];
configcommon.epoch.toi{3}              = [1, -2];
configcommon.epoch.pad{1}              = 10;
configcommon.epoch.pad{2}              = 10;
configcommon.epoch.pad{3}              = 0.5;



%% Patient 1
%REVOIR CHANNEL ALIGNEMENT DE CE PREMIER ENREGISTREMENT
%ADD EOG CHANNEL IF ANY
%LFP.NAME : toujours commencer par right
%VERIF ALL ELECTRODES. PAS RAJOUTER CELLES EN PLUS
%DEFINIR CONGIG_COMMON.ALIGN .LFP .EMG .EPOCH
config{1}                           = configcommon;
config{1}.prefix                    = 'pat_LGI1_008-EEG_129-';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{1}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
config{1}.directorylist{1}          = {'EEG_129'}; %dir = eeg file with all the electrodess
config{1}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker. For defining trials
config{1}.align.name                = {'SlowWave_R'};%{'SlowWave_R','SlowWave_R'};
config{1}.align.channel             = {'C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'}; %alwa
config{1}.LFP.emg                   = {'EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{1}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

%% Patient 2
%Same patient, other EEG.
config{2}                           = configcommon;
config{2}.prefix                    = 'pat_LGI1_008-EEG_131-';
config{2}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{2}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
config{2}.directorylist{1}          = {'EEG_131'}; %dir = eeg file with all the electrodess
config{2}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{2}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{2}.align.channel             = {'C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{2}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{2}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 3
config{3}                           = configcommon;
config{3}.prefix                    = 'pat_LGI1_001-EEG_1-';
config{3}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{3}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{3}.directorylist{1}          = {'EEG_1'}; %dir = eeg file with all the electrodess
config{3}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{3}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{3}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{3}.LFP.name                  = {'SlowWave_R','SlowWave_L'};
config{3}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{3}.continous                 = false; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 4
config{4}                           = configcommon;
config{4}.prefix                    = 'pat_LGI1_001-EEG_14-';
config{4}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{4}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{4}.directorylist{1}          = {'EEG_14'}; %dir = eeg file with all the electrodess
config{4}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{4}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{4}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{4}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{4}.LFP.emg                   = {'EMG1+','EMG2+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{4}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 5
config{5}                           = configcommon;
config{5}.prefix                    = 'pat_LGI1_001-EEG_10-';
config{5}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{5}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{5}.directorylist{1}          = {'EEG_10'}; %dir = eeg file with all the electrodess
config{5}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{5}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{5}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{5}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{5}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{5}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 6
config{6}                           = configcommon;
config{6}.prefix                    = 'pat_LGI1_007-EEG_113-';
config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
config{6}.directorylist{1}          = {'EEG_113'}; %dir = eeg file with all the electrodess
config{6}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{6}.align.name                = {'SlowWave_R'};%,'SlowWave_L'};
config{6}.align.channel             = {'C4'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{6}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'};%,'SlowWave_L'};
config{6}.LFP.emg                   = {'EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{6}.continous                 = false; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 7
config{7}                           = configcommon;
config{7}.prefix                    = 'pat_LGI1_007-EEG_123-';
config{7}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{7}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
config{7}.directorylist{1}          = {'EEG_123'}; %dir = eeg file with all the electrodess
config{7}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{7}.align.name                = {'SlowWave_R'};%,'SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{7}.align.channel             = {'C4'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{7}.LFP.name                  = {'SlowWave_R'};%,'SlowWave_L'};%,'SlowWave_L'};
config{7}.LFP.emg                   = {'EMG1+'};%,'no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{7}.continous                 = false; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 8
config{8}                           = configcommon;
config{8}.prefix                    = 'pat_LGI1_009-EEG_156-';
config{8}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_009');
config{8}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_009');       % where to print images
config{8}.directorylist{1}          = {'EEG_156'}; %dir = eeg file with all the electrodess
config{8}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{8}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{8}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{8}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{8}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{8}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 9
config{9}                           = configcommon;
config{9}.prefix                    = 'pat_LGI1_010-EEG_172-';
config{9}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{9}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{9}.directorylist{1}          = {'EEG_172'}; %dir = eeg file with all the electrodess
config{9}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{9}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{9}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{9}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{9}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{9}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 10
config{10}                           = configcommon;
config{10}.prefix                    = 'pat_LGI1_010-EEG_203-';
config{10}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{10}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{10}.directorylist{1}          = {'EEG_203'}; %dir = eeg file with all the electrodess
config{10}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{10}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{10}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{10}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{10}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{10}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 11
config{11}                           = configcommon;
config{11}.prefix                    = 'pat_LGI1_010-EEG_174-';
config{11}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{11}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{11}.directorylist{1}          = {'EEG_174'}; %dir = eeg file with all the electrodess
config{11}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{11}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{11}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{11}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{11}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{11}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 12
config{12}                           = configcommon;
config{12}.prefix                    = 'pat_LGI1_010-EEG_177-';
config{12}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{12}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{12}.directorylist{1}          = {'EEG_177'}; %dir = eeg file with all the electrodess
config{12}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{12}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{12}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{12}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{12}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{12}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 13
config{13}                           = configcommon;
config{13}.prefix                    = 'pat_LGI1_010-EEG_180-';
config{13}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{13}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{13}.directorylist{1}          = {'EEG_180'}; %dir = eeg file with all the electrodess
config{13}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{13}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{13}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{13}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{13}.LFP.emg                   = {'EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{13}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 14
config{14}                           = configcommon;
config{14}.prefix                    = 'pat_LGI1_010-EEG_199-';
config{14}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{14}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{14}.directorylist{1}          = {'EEG_199'}; %dir = eeg file with all the electrodess
config{14}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{14}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{14}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{14}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{14}.LFP.emg                   = {'no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{14}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 15
config{15}                           = configcommon;
config{15}.prefix                    = 'pat_LGI1_012-EEG_228-';
config{15}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
config{15}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
config{15}.directorylist{1}          = {'EEG_228'}; %dir = eeg file with all the electrodess
config{15}.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{15}.align.name                = {'SlowWave_R','SlowWave_L'};%{'SlowWave_R','SlowWave_R'};
config{15}.align.channel             = {'F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{15}.LFP.name                  = {'SlowWave_R','SlowWave_L'};%,'SlowWave_L'};
config{15}.LFP.emg                   = {'no','EMG1+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{15}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

end


