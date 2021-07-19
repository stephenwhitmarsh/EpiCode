
%% Setting parameters DTX project Paul Baudin

function [config] = dtx_patients_lgi1_setparams(config)

disp('setting parameters for LGI1 patients');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-Patients_LGI1';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/Patients_LGI1';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-Patients_LGI1\';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\Patients_LGI1\';
else
    error('Platform not supported')
end


%% Congig common for all patients :

%where to save data and images
datasavedir  = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis, 'image');

%general infos
configcommon.name                      = {'SlowWave_R','SlowWave_L','SlowWave_R_EMG_align','SlowWave_L_EMG_align'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.labels.macro              = {'Fp2','F4','C4','P4','O2','F8','T4','T6','Fpz','Fz','Cz','Pz','Oz','Fp1','F3',...
    'C3','P3','O1','F7','T3','T5'}; %the order is important. Plot first on top.
configcommon.LFP.channel               = configcommon.labels.macro; 
configcommon.muse.backupdir               	= fullfile(rootpath_analysis, 'Musemarkers_backup');

%cut into trials according to Muse markers with readLFP.m
configcommon.muse.startmarker.SlowWave_R    = 'SlowWave_R';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_R      = 'SlowWave_R';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_R           = [-10, 10];
configcommon.epoch.pad.SlowWave_R           = 0;

configcommon.muse.startmarker.SlowWave_L    = 'SlowWave_L';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_L      = 'SlowWave_L';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_L           = [-10, 10];
configcommon.epoch.pad.SlowWave_L           = 0;

configcommon.muse.startmarker.SlowWave_R_EMG_align    = 'SlowWave_R_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_R_EMG_align      = 'SlowWave_R_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_R_EMG_align           = [-10, 10];
configcommon.epoch.pad.SlowWave_R_EMG_align           = 0;

configcommon.muse.startmarker.SlowWave_L_EMG_align    = 'SlowWave_L_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_L_EMG_align      = 'SlowWave_L_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_L_EMG_align           = [-10, 10];
configcommon.epoch.pad.SlowWave_L_EMG_align           = 0;

%preprocessing of the LFP and EMG in readLFP.m, and som einformations required for
%plots of the LFP and EMG
configcommon.LFP.name                  = {'SlowWave_R','SlowWave_L','SlowWave_R_EMG_align','SlowWave_L_EMG_align'};
configcommon.LFP.write                 = true;
configcommon.LFP.motorcortex.SlowWave_L  = 'C3';
configcommon.LFP.motorcortex.SlowWave_R  = 'C4';
configcommon.SlowWave_R.channel = {'Fp2','F4','C4','F8','Fz','Cz'};
configcommon.SlowWave_L.channel = {'Fp1','F3','C3','F7','Fz','Cz'};

configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 0;
configcommon.LFP.resamplefs            = 256;
configcommon.LFP.baseline              = 'no';%No baseline when reading LFP. Baseline applied for each trial after trial definition
configcommon.LFP.reref                 = 'yes';
configcommon.LFP.rerefmethod           = 'avg';
configcommon.LFP.refchannel            = configcommon.LFP.channel;
configcommon.LFP.bsfilter              = 'no';%no need because of the low pass filter
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.lpfilter              = 'yes';
configcommon.LFP.lpfreq                = 20;
configcommon.LFP.lpfilttype            = 'fir';

configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];
configcommon.EMG.bsfiltord             = 3;
configcommon.EMG.reref                 = 'no';
configcommon.EMG.rerefmethod           = [];
configcommon.EMG.refchannel            = [];
configcommon.EMG.envmethod             = 'rms';
configcommon.EMG.envparam              = 50;

configcommon.minbadtime.SlowWave_R_EMG_align   = 0;
configcommon.minbadtime.SlowWave_R             = 0;
configcommon.minbadtime.SlowWave_L_EMG_align   = 0;
configcommon.minbadtime.SlowWave_L             = 0;

%aligment of each trial to the peak of the TDW, with alignMuseMarkerPeaks.m
configcommon.align.name                = {'SlowWave_R','SlowWave_L'};
configcommon.align.reref               = 'yes';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = configcommon.labels.macro';
configcommon.align.notch               = 'yes';

configcommon.align.method.SlowWave_R              = 'nearestmin';      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter.SlowWave_R              = 'lp';
configcommon.align.freq.SlowWave_R                = 2;          
configcommon.align.demean.SlowWave_R              = 'yes';
configcommon.align.thresh.value.SlowWave_R        = 1;
configcommon.align.thresh.method.SlowWave_R       = 'trial';%'medianbl','both';
configcommon.align.maxtimeshift.SlowWave_R        = 0.8;
configcommon.align.toiplot.SlowWave_R             = [-1.5,  1]; 
configcommon.align.toiactive.SlowWave_R           = [-1, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline.SlowWave_R         = [-10, -1];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

configcommon.align.method.SlowWave_L              = 'nearestmin';      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter.SlowWave_L              = 'lp';
configcommon.align.freq.SlowWave_L                = 2;          
configcommon.align.demean.SlowWave_L              = 'yes';
configcommon.align.thresh.value.SlowWave_L        = 1;
configcommon.align.thresh.method.SlowWave_L       = 'trial';%'medianbl','both';
configcommon.align.maxtimeshift.SlowWave_L        = 0.8;
configcommon.align.toiplot.SlowWave_L             = [-1.5,  1]; 
configcommon.align.toiactive.SlowWave_L           = [-1, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline.SlowWave_L         = [-10, -1];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

%for amplitude and halfwidth detection, with plot_morpho.m
configcommon.morpho.channame           = []; %set for each patient, use align.channel
configcommon.morpho.negpeak            = 'yes';
configcommon.morpho.toiplot            = [-2 2];
configcommon.morpho.measurehalfwidth   = 'yes';
configcommon.morpho.measureamplitude   = 'yes';
configcommon.morpho.blmethod           = 'bl';
configcommon.morpho.toibl              = [-2 -1];

%% config specific for each patient :

%% Patient 1
%general infos
config{1}                           = configcommon;
config{1}.prefix                    = 'pat_LGI1_001-';
config{1}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_001');
config{1}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_001');       
config{1}.directorylist{1}          = {'EEG_14','EEG_10'}; %dir = eeg file with all the electrodess

%electrode used for peak alignment and morphology measurements
config{1}.align.channel.SlowWave_R  = 'F4';         
config{1}.align.channel.SlowWave_L  = 'F3';  

%name of the EMG electrodes
config{1}.EMG.SlowWave_R            = 'EMG1+';
config{1}.EMG.SlowWave_L            = 'EMG2+';
config{1}.EMG.SlowWave_R_EMG_align  = 'EMG1+';
config{1}.EMG.SlowWave_L_EMG_align  = 'EMG2+';

%time used for computing TDW morpho
config{1}.morpho.toi.SlowWave_R          = [-0.5 0.3];
config{1}.morpho.toi.SlowWave_L          = [-0.4 0.5];
config{1}.morpho.toibl                   = [-1 -0.5];

%electrodes where the TDW is visualy detected (only for plot)
config{1}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{1}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{1}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7'};
config{1}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3'};

%% Patient 2
config{2}                           = configcommon;
config{2}.prefix                    = 'pat_LGI1_002-';
config{2}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_002');
config{2}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_002');       
config{2}.directorylist{1}          = {'EEG_34','EEG_37','EEG_39'}; 

config{2}.align.channel.SlowWave_R   = 'F4';
config{2}.align.channel.SlowWave_L   = 'F3';

config{2}.EMG.SlowWave_R            = [];
config{2}.EMG.SlowWave_L            = 'EMG1+';
config{2}.EMG.SlowWave_R_EMG_align  = [];
config{2}.EMG.SlowWave_L_EMG_align  = 'EMG1+';

config{2}.morpho.toi.SlowWave_R          = [-0.4 0.8];
config{2}.morpho.toi.SlowWave_L          = [-0.5 0.4];

config{2}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{2}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{2}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{2}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};

%% Patient 3
config{3}                           = configcommon;
config{3}.prefix                    = 'pat_LGI1_007-';
config{3}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_007');
config{3}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_007');      
config{3}.directorylist{1}          = {'EEG_113','EEG_123'}; 
config{3}.continous                 = false; %For this patient, EEG data are cut and clinicians onlys keep seizures

config{3}.align.channel.SlowWave_R  = 'C4';
config{3}.align.channel.SlowWave_L  = [];
config{3}.align.freq.SlowWave_R     = 5;

config{3}.EMG.SlowWave_R            = 'EMG1+';
config{3}.EMG.SlowWave_L            = [];
config{3}.EMG.SlowWave_R_EMG_align  = 'EMG1+';
config{3}.EMG.SlowWave_L_EMG_align  = [];

config{3}.morpho.toi.SlowWave_R          = [-0.3 0];
config{3}.morpho.toi.SlowWave_L          = [];

config{3}.morpho.channel.LFP.SlowWave_R  = {'Fz','F4','F8','Cz','C4'};
config{3}.morpho.channel.CSD.SlowWave_R  = {'Fz','F8','Cz','C4'};
config{3}.morpho.channel.LFP.SlowWave_L  = {};
config{3}.morpho.channel.CSD.SlowWave_L  = {};

%% Patient 4
config{4}                           = configcommon;
config{4}.prefix                    = 'pat_LGI1_008-';
config{4}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_008');
config{4}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_008');      
config{4}.directorylist{1}          = {'EEG_129', 'EEG_131'};

config{4}.align.channel.SlowWave_R  = 'F4';
config{4}.align.channel.SlowWave_L  = 'F3';

config{4}.EMG.SlowWave_R            = 'EMG1+';
config{4}.EMG.SlowWave_L            = [];
config{4}.EMG.SlowWave_R_EMG_align  = 'EMG1+';
config{4}.EMG.SlowWave_L_EMG_align  = [];

config{4}.morpho.toi.SlowWave_R          = [-0.6 0.4];
config{4}.morpho.toi.SlowWave_L          = [-0.5 0.5];

config{4}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{4}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','Cz','C4'};
config{4}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{4}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};

%% Patient 5
config{5}                           = configcommon;
config{5}.prefix                    = 'pat_LGI1_009-';
config{5}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_009');
config{5}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_009');      
config{5}.directorylist{1}          = {'EEG_156'}; 

config{5}.align.channel.SlowWave_R  = 'F4';
config{5}.align.channel.SlowWave_L  = 'F3';

config{5}.EMG.SlowWave_R            = 'EMG1+';
config{5}.EMG.SlowWave_L            = [];
config{5}.EMG.SlowWave_R_EMG_align  = 'EMG1+';
config{5}.EMG.SlowWave_L_EMG_align  = [];

config{5}.morpho.toi.SlowWave_R          = [-0.45 0.6];
config{5}.morpho.toi.SlowWave_L          = [-0.55 0.5];

config{5}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{5}.morpho.channel.CSD.SlowWave_R  = {'Fp2','F4','F8',};
config{5}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{5}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};

%% Patient 6
config{6}                           = configcommon;
config{6}.prefix                    = 'pat_LGI1_010-';
config{6}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_010');
config{6}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_010');       
config{6}.directorylist{1}          = {'EEG_172','EEG_203','EEG_174','EEG_177','EEG_180','EEG_199'}; 

config{6}.align.channel.SlowWave_R  = 'F4';
config{6}.align.channel.SlowWave_L  = 'F3';

config{6}.EMG.SlowWave_R            = 'EMG1+';
config{6}.EMG.SlowWave_L            = [];
config{6}.EMG.SlowWave_R_EMG_align  = 'EMG1+';
config{6}.EMG.SlowWave_L_EMG_align  = [];

config{6}.morpho.toi.SlowWave_R          = [-0.7 1];
config{6}.morpho.toi.SlowWave_L          = [-0.7 0.7];

config{6}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{6}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{6}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{6}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};

%% Patient 7
config{7}                           = configcommon;
config{7}.prefix                    = 'pat_LGI1_012-';
config{7}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_012');
config{7}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_012');       
config{7}.directorylist{1}          = {'EEG_228'}; 

config{7}.align.channel.SlowWave_R  = 'F4';
config{7}.align.channel.SlowWave_L  = 'F3';

config{7}.EMG.SlowWave_R            = [];
config{7}.EMG.SlowWave_L            = [];
config{7}.EMG.SlowWave_R_EMG_align  = [];
config{7}.EMG.SlowWave_L_EMG_align  = [];

config{7}.morpho.toi.SlowWave_R          = [-0.6 0.6];
config{7}.morpho.toi.SlowWave_L          = [-0.7 0.6];

config{7}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{7}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{7}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','C3'};
config{7}.morpho.channel.CSD.SlowWave_L  = {'Fp1','F3','F7','C3'};

%% Patient 8
config{8}                           = configcommon;
config{8}.prefix                    = 'pat_LGI1_016-';
config{8}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_016');
config{8}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_016');      
config{8}.directorylist{1}          = {'EEG_1786'}; 

config{8}.align.channel.SlowWave_R  = [];
config{8}.align.channel.SlowWave_L  = 'Cz';

config{8}.EMG.SlowWave_R            = [];
config{8}.EMG.SlowWave_L            = 'EMG1+';
config{8}.EMG.SlowWave_R_EMG_align  = [];
config{8}.EMG.SlowWave_L_EMG_align  = 'EMG1+';

config{8}.morpho.toi.SlowWave_R          = [];
config{8}.morpho.toi.SlowWave_L          = [-0.5 0.3];

config{8}.morpho.channel.LFP.SlowWave_R  = {};
config{8}.morpho.channel.CSD.SlowWave_R  = {};
config{8}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{8}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};

%% Patient 9
config{9}                           = configcommon;
config{9}.prefix                    = 'pat_LGI1_017-';
config{9}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_017');
config{9}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_017');       
config{9}.directorylist{1}          = {'2010.12.15_1019','2010.12.30_1107','2010.12.31_1324','2011.05.27_1328'};

config{9}.align.channel.SlowWave_R  = 'F4';
config{9}.align.channel.SlowWave_L  = 'F3';

config{9}.EMG.SlowWave_R            = [];
config{9}.EMG.SlowWave_L            = 'EMG Chin';
config{9}.EMG.SlowWave_R_EMG_align  = [];
config{9}.EMG.SlowWave_L_EMG_align  = 'EMG Chin';

config{9}.morpho.toi.SlowWave_R          = [-0.5 0.5];
config{9}.morpho.toi.SlowWave_L          = [-0.5 0.5];

config{9}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{9}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{9}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{9}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3','Cz','C3'};

%% Patient 10
config{10}                           = configcommon;
config{10}.prefix                    = 'pat_LGI1_018-';
config{10}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_018');
config{10}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_018');       
config{10}.directorylist{1}          = {'2007.12.12_1127','2007.12.13_1128','2007.12.20_1150','2007.12.28_1107','2007.12.31_1023','2008.01.29_0940'}; 

config{10}.align.channel.SlowWave_R  = 'F4';
config{10}.align.channel.SlowWave_L  = 'F3';

config{10}.EMG.SlowWave_R            = [];
config{10}.EMG.SlowWave_L            = [];
config{10}.EMG.SlowWave_R_EMG_align  = [];
config{10}.EMG.SlowWave_L_EMG_align  = [];

config{10}.morpho.toi.SlowWave_R          = [-0.4 0.4];
config{10}.morpho.toi.SlowWave_L          = [-0.2 0.3];

config{10}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','Cz','C4'};
config{10}.morpho.channel.CSD.SlowWave_R  = {'Fp2','Fz','F4'};
config{10}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{10}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','F3'};

%% Patient 11
config{11}                           = configcommon;
config{11}.prefix                    = 'pat_LGI1_019-';
config{11}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_019');
config{11}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_019');      
config{11}.directorylist{1}          = {'2009.03.26_0937','2009.11.03_0936','2009.12.08_1114','2009.12.08_1114_2','2009.12.14_1337',...
    '2009.12.17_1458','2010.01.04_1426','2010.01.06_1335','2010.03.18_1431','2010.06.17_1501','2011.04.29_1354','2011.11.02_1455','2012.03.21_1023'}; 

config{11}.align.channel.SlowWave_R = 'F4';
config{11}.align.channel.SlowWave_L = 'F3';

config{11}.EMG.SlowWave_R            = [];
config{11}.EMG.SlowWave_L            = 'EMG Chin';
config{11}.EMG.SlowWave_R_EMG_align  = [];
config{11}.EMG.SlowWave_L_EMG_align  = 'EMG Chin';

config{11}.morpho.toi.SlowWave_R          = [-0.7 0.7];
config{11}.morpho.toi.SlowWave_L          = [-0.5 0.5];

config{11}.morpho.channel.LFP.SlowWave_R  = {'Fp2','Fz','F4','F8','C4'};
config{11}.morpho.channel.CSD.SlowWave_R  = {'Fp2','F4','F8','C4'};
config{11}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{11}.morpho.channel.CSD.SlowWave_L  = {'Fp1','F3','C3'};

%% Patient 12
config{12}                           = configcommon;
config{12}.prefix                    = 'pat_LGI1_020-';
config{12}.rawdir                    = fullfile(rootpath_data,'pat_LGI1_020');
config{12}.imagesavedir              = fullfile(imagesavedir,'pat_LGI1_020');      
config{12}.directorylist{1}          = {'2011.09.14.15-24','2011.09.15.13-46','2011.09.19.10-46','2011.09.21.10-11','2011.09.23.11-26'}; 

config{12}.align.channel.SlowWave_R  = [];
config{12}.align.channel.SlowWave_L  = 'Cz';

config{12}.EMG.SlowWave_R            = [];
config{12}.EMG.SlowWave_L            = 'EMG1+';
config{12}.EMG.SlowWave_R_EMG_align  = [];
config{12}.EMG.SlowWave_L_EMG_align  = 'EMG1+';

config{12}.morpho.toi.SlowWave_R          = [];
config{12}.morpho.toi.SlowWave_L          = [-0.9 0.2];

config{12}.morpho.channel.LFP.SlowWave_R  = {};
config{12}.morpho.channel.CSD.SlowWave_R  = {};
config{12}.morpho.channel.LFP.SlowWave_L  = {'Fp1','Fz','F3','F7','Cz','C3'};
config{12}.morpho.channel.CSD.SlowWave_L  = {'Fp1','Fz','Cz','C3'};

end


