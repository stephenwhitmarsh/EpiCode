
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file)

function [config] = dtx_setparams_patients_lgi1(config)

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
%mergeindex = {[1, 2], [3, 4, 5], [6, 7],[8], [9, 10, 11, 12, 13, 14],[15]};

configcommon.os                        = os;
configcommon.name                      = {'SlowWave_R','SlowWave_L','SlowWave_R_EMG','SlowWave_L_EMG'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'};
configcommon.muse.startend             = {'SlowWave_R','SlowWave_R'; 'SlowWave_L','SlowWave_L'; 'SlowWave_R_EMG','SlowWave_R_EMG'; 'SlowWave_L_EMG','SlowWave_L_EMG'};   % start and end Muse marker. For defining trials
configcommon.merge                     = true; %if "true" : the last "ipart" is a merge of all the previous parts.Usefull for good naming of the analysis files

configcommon.align.name                = {'SlowWave_R','SlowWave_L','SlowWave_R_EMG','SlowWave_L_EMG'};
configcommon.align.flip                = {'no','no'};
configcommon.align.abs                 = {'no','no'};
configcommon.align.method              = {'nearestmin','nearestmin','nearestmin','nearestmin'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp','lp','lp','lp'};
configcommon.align.freq                = {5,5,5,5};          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no','no','no','no'};
configcommon.align.thresh              = [0, 0, 0, 0];
configcommon.align.toiplot             = {[-1,  1], [-1,  1], [-1,  1], [-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.toiactive           = {[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline         = {[-1, -0.5], [-1, -0.5], [-1, -0.5], [-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.reref               = 'yes';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = configcommon.labels.macro';
configcommon.align.notch               = 'yes';

configcommon.LFP.name                  = {'SlowWave_R','SlowWave_L','SlowWave_R_EMG','SlowWave_L_EMG'};
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 0;
configcommon.LFP.resamplefs            = 256;
configcommon.LFP.baseline              = 'no';
configcommon.LFP.baselinewindow        = {[-2, -1], [-2, -1], [-2, -1], [-2, -1]};
configcommon.LFP.slidestep             = 0.01;
configcommon.LFP.reref                 = 'yes';
configcommon.LFP.rerefmethod           = 'avg';
configcommon.LFP.refchannel            = configcommon.labels.macro';
configcommon.LFP.bsfilter              = 'no';
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.lpfilter              = 'yes';
configcommon.LFP.lpfreq                = 30;
configcommon.LFP.lpfilttype            = 'fir';
configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];

%1, 2 and 3 associated with config.muse.startend
configcommon.epoch.toi{1}              = [-10, 10];
configcommon.epoch.toi{2}              = [-10, 10];
configcommon.epoch.toi{3}              = [-10, 10];
configcommon.epoch.toi{4}              = [-10, 10];
configcommon.epoch.pad{1}              = 10;
configcommon.epoch.pad{2}              = 10;
configcommon.epoch.pad{3}              = 10;
configcommon.epoch.pad{4}              = 10;


%% Patient 1
config{1}                           = configcommon;
config{1}.prefix                    = 'pat_LGI1_008-';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{1}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');       % where to print images
config{1}.directorylist{1}          = {'EEG_129'}; %dir = eeg file with all the electrodess
config{1}.directorylist{2}          = {'EEG_131'};
config{1}.directorylist{3}          = {'EEG_129', 'EEG_131'};
config{1}.align.channel             = {'C4','Fp1','C4','Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.LFP.emg                   = {'no','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{1}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 2
config{2}                           = configcommon;
config{2}.prefix                    = 'pat_LGI1_001-';
config{2}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{2}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       % where to print images
config{2}.directorylist{1}          = {'EEG_14'}; %dir = eeg file with all the electrodess
config{2}.directorylist{2}          = {'EEG_10'}; %dir = eeg file with all the electrodess
config{2}.directorylist{3}          = {'EEG_14','EEG_10'}; %dir = eeg file with all the electrodess
config{2}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.LFP.emg                   = {'no','no','EMG1+','EMG2+'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{2}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 3
config{3}                           = configcommon;
config{3}.prefix                    = 'pat_LGI1_007-';
config{3}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{3}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');       % where to print images
config{3}.directorylist{1}          = {'EEG_113'}; %dir = eeg file with all the electrodess
config{3}.directorylist{2}          = {'EEG_123'}; %dir = eeg file with all the electrodess
config{3}.directorylist{3}          = {'EEG_113','EEG_123'}; %dir = eeg file with all the electrodess
config{3}.align.channel             = {'C4','no','C4','no'};%,'Fp1'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{3}.LFP.emg                   = {'no','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{3}.continous                 = false; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 4
config{4}                           = configcommon;
config{4}.prefix                    = 'pat_LGI1_009-';
config{4}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_009');
config{4}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_009');       % where to print images
config{4}.directorylist{1}          = {'EEG_156'}; %dir = eeg file with all the electrodess
config{4}.align.channel             = {'F4','F3','F4','F3'}; % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{4}.LFP.emg                   = {'no','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{4}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures


%% Patient 5
config{5}                           = configcommon;
config{5}.prefix                    = 'pat_LGI1_010-';
config{5}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{5}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       % where to print images
config{5}.directorylist{1}          = {'EEG_172'}; %dir = eeg file with all the electrodess
config{5}.directorylist{2}          = {'EEG_203'}; %dir = eeg file with all the electrodess
config{5}.directorylist{3}          = {'EEG_174'}; %dir = eeg file with all the electrodess
config{5}.directorylist{4}          = {'EEG_177'}; %dir = eeg file with all the electrodess
config{5}.directorylist{5}          = {'EEG_180'}; %dir = eeg file with all the electrodess
config{5}.directorylist{6}          = {'EEG_199'}; %dir = eeg file with all the electrodess
config{5}.directorylist{7}          = {'EEG_172','EEG_203','EEG_174','EEG_177','EEG_180','EEG_199'}; %dir = eeg file with all the electrodess
config{5}.align.channel             = {'F4','F3','F4','F3'};       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{5}.LFP.emg                   = {'no','no','EMG1+','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{5}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

%% Patient 6
config{6}                           = configcommon;
config{6}.prefix                    = 'pat_LGI1_012-';
config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       % where to print images
config{6}.directorylist{1}          = {'EEG_228'}; %dir = eeg file with all the electrodess
config{6}.align.channel             = {'F4','F3','F4','F3'};%{'C4','C3');       % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{6}.LFP.emg                   = {'no','no','no','no'};%same index as associated EEG. 'no' if no EMG associated to this seizure side 
config{6}.continous                 = true; %sometimes EEG data are cut and clinicians onlys keep seizures

end


