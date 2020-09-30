
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file)

function [config] = dtx_setparams_patients_lgi1(config)

disp('setting parameters for LGI1 patients');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-Patients_LGI1';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/Patients_LGI1';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-Patients_LGI1\';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\Patients_LGI1\';
    os                  = 'windows';
else
    error('Platform not supported')
end



%% Congig common for all patients

datasavedir = fullfile(rootpath_analysis, 'data', 'spike');
imagesavedir = fullfile(rootpath_analysis, 'image_lfp');
%mergeindex = {[1, 2], [3, 4, 5], [6, 7],[8], [9, 10, 11, 12, 13, 14],[15]};

configcommon.os                        = os;
% SlowWave_R_begin, SlowWave_L_begin
% SlowWave_R_peak, SlowWave_L_peak
% configcommon.name                      = {'SlowWave_R_begin','SlowWave_L_begin','SlowWave_R_EMGalign','SlowWave_L_EMGalign','SlowWavealign_EMG_R','SlowWavealign_EMG_L'};
configcommon.name                      = {'SlowWave_R','SlowWave_L'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
configcommon.LFP.channel               = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
% configcommon.merge                     = true; %if "true" : the last "ipart" is a merge of all the previous parts.Usefull for good naming of the analysis files

configcommon.muse.backupdir                   = fullfile(rootpath_analysis, 'Musemarkers_backup');
configcommon.muse.startend(1,1:2)             = {'SlowWave_R','SlowWave_R'};   % start and end Muse marker. For defining trials
configcommon.muse.startend(2,1:2)             = {'SlowWave_L','SlowWave_L'};   % start and end Muse marker. For defining trials
% configcommon.muse.startend(3,1:2)             = {'SlowWave_R_EMG','SlowWave_R_EMG'};   % start and end Muse marker. For defining trials
% configcommon.muse.startend(4,1:2)             = {'SlowWave_L_EMG','SlowWave_L_EMG'};   % start and end Muse marker. For defining trials
% configcommon.muse.startend(5,1:2)             = {'SlowWave_R','SlowWave_R'};   % start and end Muse marker. For defining trials
% configcommon.muse.startend(6,1:2)             = {'SlowWave_L','SlowWave_L'};   % start and end Muse marker. For defining trials

%Index are associated with config.muse.startend
configcommon.epoch.toi{1}              = [-10, 10];
configcommon.epoch.toi{2}              = [-10, 10];
% configcommon.epoch.toi{3}              = [-10, 10];
% configcommon.epoch.toi{4}              = [-10, 10];
% configcommon.epoch.toi{5}              = [-10, 10];
% configcommon.epoch.toi{6}              = [-10, 10];
configcommon.epoch.pad{1}              = 10;
configcommon.epoch.pad{2}              = 10;
% configcommon.epoch.pad{3}              = 10;
% configcommon.epoch.pad{4}              = 10;
% configcommon.epoch.pad{5}              = 10;
% configcommon.epoch.pad{6}              = 10;

% ALIGN XCORR : 
% configcommon.align.name                = {'SlowWave_R', 'SlowWave_L'};                               % Name of markers/patterns to align
% configcommon.align.channel             = [];    %Setup for each patient separately
% configcommon.align.demean              = 'yes';
% configcommon.align.baselinewindow      = [-2.5 -1.5];
% configcommon.align.reref               = 'no';
% configcommon.align.refmethod           = 'bipolar';
% configcommon.align.latency             = [-1.5, 1];  

% OLD ALIGN
%for peak detection
configcommon.align.name                = {'SlowWave_R','SlowWave_L'};
% configcommon.align.flip                = {'no','no'};
configcommon.align.abs                 = {'no','no','no','no'};
configcommon.align.method              = {'nearestmin','nearestmin','nearestmin','nearestmin'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp','lp','lp','lp'};
configcommon.align.freq                = {5,5,5,5};          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no','no','no','no'};
configcommon.align.thresh.value        = [1, 1, 1, 1];
configcommon.align.thresh.method       = {'trial','trial','trial','trial'};%'medianbl','both';
configcommon.align.toiplot             = {[-1.5,  1],[-1.5,  1],[-1.5,  1],[-1.5,  1]}; 
configcommon.align.toiactive           = {[-1, 0.5], [-1, 0.5], [-1, 0.5], [-1, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline         = {[-11, -1], [-11, -1], [-11, -1], [-11, -1], [-11, -1], [-11, -1]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.reref               = 'yes';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = configcommon.labels.macro';
configcommon.align.notch               = 'yes';
configcommon.align.maxtimeshift        = 0.8;
configcommon.align.flip                = 'yes';
configcommon.align.demean              = 'yes';
% %for detection of begin of event
% configcommon.align.begin.name          = {'SlowWave_R_begin', 'SlowWave_L_begin'};
% configcommon.align.begin.maxtimeshift  = 0.8;
% configcommon.align.begin.thresh        = 0.2; % percent of peak
% configcommon.alignEMG.name             = {'SlowWave_R_EMGalign','SlowWave_L_EMGalign'}; %name of analysis
% %configcommon.alignEMG.channel         %set separately for each patient
% configcommon.alignEMG.toiplot          = {[-1,  2],[-1,  2]}; 
% configcommon.alignEMG.toiactive        = {[-0.5, 0.5], [-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% configcommon.alignEMG.toibaseline      = {[-1, -0.5], [-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% configcommon.alignEMG.maxtimeshift     = 0.25; %if abs(timeshift) is superior, trial is considered as misdetected and remived 
% configcommon.alignEMG.hpfreq           = 10;

configcommon.LFP.name                  = {'SlowWave_R','SlowWave_L'};
configcommon.LFP.write                 = true;
% configcommon.LFP.emgmarker             = {'no','no','SlowWave_R_EMG','SlowWave_L_EMG','SlowWave_R_EMG','SlowWave_L_EMG'};
configcommon.LFP.motorcortex           = {'C4','C3'};
% configcommon.LFP.TFR.doTFR             = false;

configcommon.LFP.flip                  = true;
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 0;
configcommon.LFP.resamplefs            = 256;
configcommon.LFP.baseline              = 'no';%No baseline when reading LFP. Baseline applied for each trial after trial definition
% configcommon.LFP.baselinewindow        = {[-2, -1], [-2, -1]};
configcommon.LFP.reref                 = 'yes';
configcommon.LFP.rerefmethod           = 'avg';
configcommon.LFP.refchannel            = configcommon.LFP.channel;
configcommon.LFP.bsfilter              = 'no';
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.lpfilter              = 'yes';
configcommon.LFP.lpfreq                = 10;
configcommon.LFP.lpfilttype            = 'fir';

configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];
configcommon.EMG.bsfiltord             = 3;
configcommon.EMG.reref                 = 'no';
configcommon.EMG.rerefmethod           = [];
configcommon.EMG.refchannel            = [];
% configcommon.EMG.envmethod             = 'rms';
% configcommon.EMG.envparam              = 50;
% configcommon.EMG.toi                   = [-5 5];

configcommon.remove.plotdata           = 'yes';
configcommon.remove.keepindexes        = 'yes';

configcommon.topoplot.marker_list      = [1 2];
configcommon.topoplot.toi_topoplot     = [-0.125 0.125];
configcommon.topoplot.toi_multiplot    = [-2 2];

configcommon.morpho.channame           = []; %set for each patient
configcommon.morpho.toiplot            = [-2 2];
configcommon.morpho.measurehalfwidth    = 'yes';
configcommon.morpho.measureamplitude    = 'yes';
configcommon.morpho.blmethod           = 'bl';
configcommon.morpho.toiac              = [-1 1];
configcommon.morpho.toibl              = [-2 -1];

% configcommon.stats.toibaseline          = {[-11, -1], [-11, -1], [-11, -1], [-11, -1], [-11, -1], [-11, -1]};%{[-5 -2],[-5 -2],[-5 -2],[-5 -2],[-5 -2],[-5 -2]};
% configcommon.stats.alpha                = 0.05;
% configcommon.stats.toiploteeg           = [-10 10];
% configcommon.stats.toilatency           = [-0.3 0.3]; %for selecting electrodes
% configcommon.topoplot.timestep          = 0.1;
% configcommon.topoplot.toi               = [-0.5 0.5];

%% Patient 1
config{1}                           = configcommon;
config{1}.prefix                    = 'pat_LGI1_001-';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{1}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
% config{1}.directorylist{1}          = {'EEG_14'}; %dir = eeg file with all the electrodess
% config{1}.directorylist{2}          = {'EEG_10'}; %dir = eeg file with all the electrodess
config{1}.directorylist{1}          = {'EEG_14','EEG_10'}; %dir = eeg file with all the electrodess
config{1}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.alignEMG.channel          = {'EMG1+','EMG2+'}; 
config{1}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% config{1}.LFP.emg                   = {'no','no','EMG1+','EMG2+','EMG1+','EMG2+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{1}.EMG.channel               = {'EMG1+','EMG2+'};
config{1}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 2
config{2}                           = configcommon;
config{2}.prefix                    = 'pat_LGI1_008-';
config{2}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{2}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
% config{2}.directorylist{1}          = {'EEG_129'}; %dir = eeg file with all the electrodess
% config{2}.directorylist{2}          = {'EEG_131'};
config{2}.directorylist{1}          = {'EEG_129', 'EEG_131'};
config{2}.align.channel             = {'C4','Fp1','C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.alignEMG.channel          = {'EMG1+','no'}; 
config{2}.morpho.channame           = {'C4','Fp1','C4','Fp1','C4','Fp1','C4','Fp1'};
% config{2}.LFP.emg                   = {'no','no','EMG1+','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{2}.EMG.channel               = {'EMG1+',[]};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{2}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

% config{2}.align.toibaseline         = {[-1.5, -0.5], [-1.5, -0.5], [-1.5, -0.5], [-1.5, -0.5], [-1.5, -0.5], [-1.5, -0.5]}; 
% config{2}.align.begin.thresh        = 0.4; % percent of peak



%% Patient 3
config{3}                           = configcommon;
config{3}.prefix                    = 'pat_LGI1_007-';
config{3}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{3}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
% config{3}.directorylist{1}          = {'EEG_113'}; %dir = eeg file with all the electrodess
% config{3}.directorylist{2}          = {'EEG_123'}; %dir = eeg file with all the electrodess
config{3}.directorylist{1}          = {'EEG_113','EEG_123'}; %dir = eeg file with all the electrodess
config{3}.align.channel             = {'C4','no','C4','no'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{3}.alignEMG.channel          = {'EMG1+','no'}; 
config{3}.morpho.channame           = {'C4','no','C4','no','C4','no','C4','no'};
% config{3}.LFP.emg                   = {'no','no','EMG1+','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{3}.EMG.channel               = {'EMG1+',[]};
config{3}.continous                 = false; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 4
config{4}                           = configcommon;
config{4}.prefix                    = 'pat_LGI1_009-';
config{4}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_009');
config{4}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_009');       % where to print images
config{4}.directorylist{1}          = {'EEG_156'}; %dir = eeg file with all the electrodess
config{4}.align.channel             = {'F4','F3','F4','F3'}; % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{4}.alignEMG.channel          = {'EMG1+','no'}; 
config{4}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% config{4}.LFP.emg                   = {'no','no','EMG1+','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{4}.EMG.channel               = {'EMG1+',[]};
config{4}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 5
config{5}                           = configcommon;
config{5}.prefix                    = 'pat_LGI1_010-';
config{5}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{5}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
% config{5}.directorylist{1}          = {'EEG_172'}; %dir = eeg file with all the electrodess
% config{5}.directorylist{2}          = {'EEG_203'}; %dir = eeg file with all the electrodess
% config{5}.directorylist{3}          = {'EEG_174'}; %dir = eeg file with all the electrodess
% config{5}.directorylist{4}          = {'EEG_177'}; %dir = eeg file with all the electrodess
% config{5}.directorylist{5}          = {'EEG_180'}; %dir = eeg file with all the electrodess
% config{5}.directorylist{6}          = {'EEG_199'}; %dir = eeg file with all the electrodess
config{5}.directorylist{1}          = {'EEG_172','EEG_203','EEG_174','EEG_177','EEG_180','EEG_199'}; %dir = eeg file with all the electrodess
config{5}.align.channel             = {'F4','F3','F4','F3'};       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{5}.alignEMG.channel          = {'EMG1+','no'}; 
config{5}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% config{5}.LFP.emg                   = {'no','no','EMG1+','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{5}.EMG.channel               = {'EMG1+',[]};
config{5}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

%% Patient 6
config{6}                           = configcommon;
config{6}.prefix                    = 'pat_LGI1_012-';
config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
config{6}.directorylist{1}          = {'EEG_228'}; %dir = eeg file with all the electrodess
config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{6}.alignEMG.channel          = {'no','no'}; 
config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{6}.EMG.channel               = {[],[]};
config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 9
config{7}                           = configcommon;
config{7}.prefix                    = 'pat_LGI1_015-';
config{7}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_015');
config{7}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_015');       % where to print images
config{7}.directorylist{1}          = {'EEG_408682','EEG_408822','EEG_408823','EEG_408825','EEG_408826','EEG_408830',...
    'EEG_408834','EEG_408836','EEG_408839'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures
%% Patient 10
config{8}                           = configcommon;
config{8}.prefix                    = 'pat_LGI1_016-';
config{8}.rawdir                     = fullfile(rootpath_data,'pat_LGI1_016');
config{8}.imagesavedir               = fullfile(imagesavedir,'pat_LGI1_016');       % where to print images
config{8}.directorylist{1}           = {'EEG_1786'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures
%% Patient 11
config{9}                           = configcommon;
config{9}.prefix                    = 'pat_LGI1_017-';
config{9}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_017');
config{9}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_017');       % where to print images
config{9}.directorylist{1}          = {'2010.12.15_1019','2010.12.30_1107','2010.12.31_1324','2011.05.27_1328'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures
%% Patient 12
config{10}                           = configcommon;
config{10}.prefix                    = 'pat_LGI1_018-';
config{10}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_018');
config{10}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_018');       % where to print images
config{10}.directorylist{1}          = {'2007.12.12_1127','2007.12.13_1128','2007.12.20_1150','2007.12.28_1107','2007.12.31_1023','2008.01.29_0940'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures
%% Patient 13
config{11}                           = configcommon;
config{11}.prefix                    = 'pat_LGI1_019-';
config{11}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_019');
config{11}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_019');       % where to print images
config{11}.directorylist{1}          = {'2009.03.26_0937','2009.11.03_0936','2009.12.08_1114','2009.12.08_1114_2','2009.12.14_1337',...
    '2009.12.17_1458','2010.01.04_1426','2010.01.06_1335','2010.03.18_1431','2010.06.17_1501','2011.04.29_1354','2011.11.02_1455','2012.03.21_1023'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

%% Patient 12
config{12}                           = configcommon;
config{12}.prefix                    = 'pat_LGI1_002-';
config{12}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_002');
config{12}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_002');       % where to print images
config{12}.directorylist{1}          = {'EEG_34','EEG_37','EEG_39'}; %dir = eeg file with all the electrodess
% config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % config{6}.alignEMG.channel          = {'no','no'}; 
% config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% config{6}.EMG.channel               = {[],[]};
% config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


% 
% 
% %% Patient 7
% config{7}                           = configcommon;
% config{7}.prefix                    = 'pat_LGI1_013-';
% % config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
% % config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
% % config{6}.directorylist{1}          = {'EEG_228'}; %dir = eeg file with all the electrodess
% % config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % % config{6}.alignEMG.channel          = {'no','no'}; 
% % config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% % config{6}.EMG.channel               = {[],[]};
% % config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures
% %% Patient 8
% config{8}                           = configcommon;
% config{8}.prefix                    = 'pat_LGI1_014-';
% % config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
% % config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
% % config{6}.directorylist{1}          = {'EEG_228'}; %dir = eeg file with all the electrodess
% % config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% % % config{6}.alignEMG.channel          = {'no','no'}; 
% % config{6}.morpho.channame           = {'F4','F3','F4','F3','F4','F3','F4','F3'};
% % % config{6}.LFP.emg                   = {'no','no','no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
% % config{6}.EMG.channel               = {[],[]};
% % config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

end


