function [config] = dtx_eeganesth_setparams(config)
% Setting parameters DTX project
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

disp('setting parameters for eeg of anesthetized rats');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-EEG_ANESTH/';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-EEG_ANESTH\';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-raw\';    
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data'); %removed of config{i} so we can more easily modify it
imagesavedir = fullfile(rootpath_analysis);

%% Congig common for all ratsdatasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis);

configcommon.os                        = os;
configcommon.name                      = {'SlowWave','SlowWave_begin','Crise_End','SlowWave_not_aligned','SlowWave_EMG_begin'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.type                      = 'dtx';

configcommon.missingdata = datetime.empty;%set for rats 4 and 7 for whom data is missing (one because unplug, the other because ampli was shut down on the first night of recording)

configcommon.muse.startmarker.SlowWave      = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave        = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave             = [-10, 5];
configcommon.epoch.pad.SlowWave             = 0;

configcommon.muse.startmarker.SlowWave_begin      = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_begin        = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_begin             = [-10, 5];
configcommon.epoch.pad.SlowWave_begin             = 0;

configcommon.muse.startmarker.Crise_End      = 'Crise_End';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.Crise_End        = 'Crise_End';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.Crise_End             = [0, 5];
configcommon.epoch.pad.Crise_End             = 5;

configcommon.muse.startmarker.SlowWave_not_aligned      = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_not_aligned        = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_not_aligned             = [-10, 5];
configcommon.epoch.pad.SlowWave_not_aligned             = 0;

%not used but kept to be compatible with awake scripts
configcommon.muse.startmarker.SlowWave_EMG_begin      = 'SlowWave_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_EMG_begin        = 'SlowWave_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_EMG_begin             = [-10, 5];
configcommon.epoch.pad.SlowWave_EMG_begin             = 0;

configcommon.seizuretimings.marker_start = 'Crise_Start';
configcommon.seizuretimings.marker_end   = 'Crise_End';
configcommon.seizuretimings.analysis_start = 'Analysis_Start';
configcommon.seizuretimings.analysis_end   = 'Analysis_End';
configcommon.seizuretimings.winsize        = 3600;%s
configcommon.seizuretimings.winstep        = 1200;%s

configcommon.LFP.name                  = configcommon.name;
configcommon.labels.macro              = {};%set for each rat
configcommon.LFP.channel               = {};%set for each rat %do not put the emg channels here
configcommon.LFP.electrodetoplot       = {};%set for each rat %do not put the emg channels here
configcommon.LFP.motorcortex           = {};%set for each rat %do not put the emg channels here

configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.resamplefs            = 512;
configcommon.LFP.reref                 = 'no';
configcommon.LFP.bsfilter              = 'yes';
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.lpfilter              = 'no';
configcommon.LFP.keepcfg               = 'no';
configcommon.LFP.write                 = false; %do not save after readLFP but after having detected artefacts, flipped, and corrected baseline

configcommon.LFP.baseline                          = 'yes';
configcommon.LFP.baselinewindow.SlowWave           = [-2 -1];
configcommon.LFP.baselinewindow.SlowWave_begin     = [-2 -1];
configcommon.LFP.baselinewindow.SlowWave_EMG_begin = [-2 -1];
configcommon.LFP.baselinewindow.Crise_End          = [3 5];
configcommon.LFP.baselinewindow.SlowWave_not_aligned = [-2 -1];

configcommon.align.name                = {'SlowWave','SlowWave_begin'};
configcommon.align.reref               = 'no';
configcommon.align.channel.SlowWave             = [];%Set for each rat
configcommon.align.method.SlowWave              = [];%SET FOR EACH RAT 'nearestmin';      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter.SlowWave              = 'bp';
configcommon.align.freq.SlowWave                = [1, 2];          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.demean.SlowWave              = 'yes';
configcommon.align.thresh.value.SlowWave        = 0;
configcommon.align.thresh.method.SlowWave       = 'trial';%'medianbl','both';
configcommon.align.maxtimeshift.SlowWave        = 0.5;
configcommon.align.toiplot.SlowWave             = [-1,  1]; 
configcommon.align.toiactive.SlowWave           = [-0.5, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline.SlowWave         = [-10, -1];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];


%SlowWave_begin : only when plotting seizures
configcommon.align.reref               = 'no';
configcommon.align.channel.SlowWave_begin             = [];%Set for each rat
configcommon.align.method.SlowWave_begin              = [];%SET FOR EACH RAT 'nearestmin';      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter.SlowWave_begin              = 'bp';
configcommon.align.freq.SlowWave_begin                = [1, 2];          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.demean.SlowWave_begin              = 'yes';
configcommon.align.thresh.value.SlowWave_begin        = 0;
configcommon.align.thresh.method.SlowWave_begin       = 'trial';%'medianbl','both';
configcommon.align.maxtimeshift.SlowWave_begin        = 0.5;
configcommon.align.toiplot.SlowWave_begin             = [-1,  1]; 
configcommon.align.toiactive.SlowWave_begin           = [-0.5, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline.SlowWave_begin         = [-10, -1];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.findbegin.SlowWave_begin          = 'yes';
configcommon.align.beginthresh.SlowWave_begin        = 0.1; % percent of peak %FIXME à vérifier

configcommon.morpho.channame           = []; %set for each patient, use align.channel
configcommon.morpho.negpeak            = 'yes';
configcommon.morpho.toiplot            = [-2 2];
configcommon.morpho.measurehalfwidth   = 'yes';
configcommon.morpho.measureamplitude   = 'yes';
configcommon.morpho.blmethod           = 'bl';
configcommon.morpho.toiac              = []; %set for each rat
configcommon.morpho.toibl              = []; %set for each rat
configcommon.morpho.raw_facealpha      = 1;
configcommon.morpho.plotavg            = 'no';
configcommon.morpho.winsize            = 3600;
configcommon.morpho.winstep            = 1200;

configcommon.TFR.name                  = {'SlowWave_begin', 'Crise_End'};
configcommon.TFR.foi                   = 3:0.1:200;
configcommon.TFR.tapsmofrq             = 5.02;
configcommon.TFR.timewinsize           = 0.35;
configcommon.TFR.timewinstep           = 0.01;
configcommon.TFR.toi.SlowWave      =[-2 2];
configcommon.TFR.toi.SlowWave_begin=[-2 2];
configcommon.TFR.toi.Crise_End     =[0 5];
configcommon.TFR.toi.Baseline      =[-7 -2];

%Intra : 
%Avant DTX22 : hp filter analogue 1Hz => manips non gardées pour la morpho
%de l'OL. 
% Avant DTX31 : ancien respirateur, pbs respiratoires récurrents.

%Probes : 
% garder manips 2, 4, 5, 6, 7
% Voir les résultats probe avec ancien respirateur et si il faut les garder
% ou pas

% plot morpho sur tous et voir où placer la fin du filtre 1Hz


%% Rodent 1
%ATTENTION : SUR MUSE, MAUVAIS NOMS DE CHANNELS : CORRIGES AVEC
%DTX_CORRECTDTX2NAME
config{1}                           = configcommon;
config{1}.domorpho                  = 'no';%no because filtered at 1Hz
config{1}.prefix                    = 'DTX2-';
config{1}.rawdir                    = fullfile(rootpath_data, 'DTX2-M1-10uM', '2019_03_01_DTX-2');
config{1}.imagesavedir              = fullfile(imagesavedir, 'DTX2');       % where to print images
config{1}.directorylist{1}          =  {'2019-03-01_12-14',...
    '2019-03-01_12-33',...
    '2019-03-01_14-14',...
    '2019-03-01_16-14',...
    '2019-03-01_18-14',...
    '2019-03-01_20-14'};

config{1}.labels.macro              = {'ECoGS1'}; % /!\ Mistake during acquisition : ECoGS1 = ECoGM1G, ECoGM1 = ECoGM1D. Corrected in the code after having loaded the data (need acquisitin name for laoding)
config{1}.injectiontime             = datetime('01-Mar-2019 12:33:00');
config{1}.LFP.channel               = config{1}.labels.macro;%ECoG S1 et ECoGM1 renommés respectivement ECoGM1G et ECoGM1D, après l'étape de readLFP

config{1}.align.channel.SlowWave = 'ECoGS1';
config{1}.align.channel.SlowWave_begin = 'ECoGS1';
config{1}.align.method.SlowWave     = 'nearestmax';
config{1}.align.method.SlowWave_begin = 'nearestmax';
config{1}.LFP.flip                  = 'yes';

config{1}.LFP.electrodetoplot       = {'ECoGS1'};
config{1}.LFP.motorcortex           = {'ECoGS1'};
config{1}.plotseizure.h             = 5000;
config{1}.morpho.toiac              = [-0.5 0.5]; 
config{1}.morpho.toibl              = [-2 -0.5];

%% Rodent 2
config{2}                           = configcommon;
config{2}.domorpho                  = 'no';%no because filtered at 1Hz
config{2}.prefix                    = 'DTX4-';
config{2}.rawdir                    = fullfile(rootpath_data, 'DTX4-M1-10uM', '2019_03_08_DTX-4');
config{2}.imagesavedir              = fullfile(imagesavedir,'DTX4');       % where to print images
config{2}.directorylist{1}          =  {'2019-03-08_12-28',...
    '2019-03-08_14-05',...
    '2019-03-08_16-05',...
    '2019-03-08_18-05',...
    '2019-03-08_20-05'};

config{2}.labels.macro              = {'ECoGM1G'};
config{2}.injectiontime             = datetime('08-Mar-2019 13:21:00');
config{2}.LFP.channel               = config{2}.labels.macro;

config{2}.align.channel.SlowWave = 'ECoGM1G';
config{2}.align.channel.SlowWave_begin = 'ECoGM1G';
config{2}.align.method.SlowWave     = 'nearestmax';
config{2}.align.method.SlowWave_begin = 'nearestmax';
config{2}.LFP.flip                  = 'yes';

config{2}.LFP.electrodetoplot       = {'ECoGM1G'};
config{2}.LFP.motorcortex           = {'ECoGM1G'};
config{2}.plotseizure.h             = 5000;
config{2}.morpho.toiac              = [-0.5 0.5]; 
config{2}.morpho.toibl              = [-2 -0.5];

%% Rodent 3
config{3}                           = configcommon;
config{3}.domorpho                  = 'no';%no because filtered at 1Hz
config{3}.prefix                    = 'DTX5-';
config{3}.rawdir                    = fullfile(rootpath_data, 'DTX5-M1-10uM', '2019_03_19_DTX-5');
config{3}.imagesavedir              = fullfile(imagesavedir, 'DTX5');       % where to print images
config{3}.directorylist{1}          =  {'2019-03-19_13-18',...
    '2019-03-19_14-30',...
    '2019-03-19_16-30',...
    '2019-03-19_18-30',...
    '2019-03-19_20-30',...
    '2019-03-19_22-30',...
    '2019-03-20_00-30'};

config{3}.labels.macro              = {'ECoGM1G'};
config{3}.injectiontime             = datetime('19-Mar-2019 13:55:00');
config{3}.LFP.channel               = config{3}.labels.macro;

config{3}.align.channel.SlowWave = 'ECoGM1G';
config{3}.align.channel.SlowWave_begin = 'ECoGM1G';
config{3}.align.method.SlowWave     = 'nearestmax';
config{3}.align.method.SlowWave_begin = 'nearestmax';
config{3}.LFP.flip                  = 'yes';

config{3}.LFP.electrodetoplot       = {'ECoGM1G'};
config{3}.LFP.motorcortex           = {'ECoGM1G'};
config{3}.plotseizure.h             = 5000;
config{3}.morpho.toiac              = [-0.5 0.5]; 
config{3}.morpho.toibl              = [-2 -0.5];

%% Rodent 4
config{4}                           = configcommon;
config{4}.domorpho                  = 'no';%no because filtered at 1Hz
config{4}.prefix                    = 'DTX6-';
config{4}.rawdir                    = fullfile(rootpath_data, 'DTX6-M1-10uM', '2019_03_21_DTX-6');
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX6');       % where to print images
config{4}.directorylist{1}          =  {'2019-03-21_14-12', '2019-03-21_16-12','2019-03-21_18-12', '2019-03-21_20-12'};

config{4}.labels.macro              = {'ECoGM1G'};
config{4}.injectiontime             = datetime('21-Mar-2019 15:08:00');
config{4}.LFP.channel               = config{4}.labels.macro;

config{4}.align.channel.SlowWave = 'ECoGM1G';
config{4}.align.channel.SlowWave_begin = 'ECoGM1G';
config{4}.align.method.SlowWave     = 'nearestmax';
config{4}.align.method.SlowWave_begin = 'nearestmax';
config{4}.LFP.flip                  = 'yes';

config{4}.LFP.electrodetoplot       = {'ECoGM1G'};
config{4}.LFP.motorcortex           = {'ECoGM1G'};
config{4}.plotseizure.h             = 5000;
config{4}.morpho.toiac              = [-0.5 0.5]; 
config{4}.morpho.toibl              = [-2 -0.5];

%% Rodent 5
config{5}                           = configcommon;
config{5}.domorpho                  = 'no';%no because filtered at 1Hz
config{5}.prefix                    = 'DTX7-';
config{5}.rawdir                    = fullfile(rootpath_data, 'DTX7-M1-10uM', '2019_03_22_DTX-7');
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX7');       % where to print images
config{5}.directorylist{1}          =  {'2019-03-22_12-31',...
    '2019-03-22_14-31',...
    '2019-03-22_16-31',...
    '2019-03-22_18-31'};

config{5}.labels.macro              = {'ECoGM1G'};
config{5}.injectiontime             = datetime('22-Mar-2019 13:40:07');
config{5}.LFP.channel               = config{5}.labels.macro;

config{5}.align.channel.SlowWave = 'ECoGM1G';
config{5}.align.channel.SlowWave_begin = 'ECoGM1G';
config{5}.align.method.SlowWave     = 'nearestmax';
config{5}.align.method.SlowWave_begin = 'nearestmax';
config{5}.LFP.flip                  = 'yes';

config{5}.LFP.electrodetoplot       = {'ECoGM1G'};
config{5}.LFP.motorcortex           = {'ECoGM1G'};
config{5}.plotseizure.h             = 5000;
config{5}.morpho.toiac              = [-0.5 0.5]; 
config{5}.morpho.toibl              = [-2 -0.5];

%% Rodent 6
config{6}                           = configcommon;
config{6}.domorpho                  = 'yes';%no because filtered at 1Hz
config{6}.prefix                    = 'DTX31-';
config{6}.rawdir                    = fullfile(rootpath_data, 'DTX31_2019_07_13-INTRA', 'Brainvision');
config{6}.imagesavedir              = fullfile(imagesavedir,'DTX31');       % where to print images
config{6}.directorylist{1}          = []; %automatically set at the end of this script

config{6}.labels.macro              = {'EEG M1G'};
config{6}.injectiontime             = datetime('13-Jul-2019 15:48:53');
config{6}.LFP.channel               = config{6}.labels.macro;

config{6}.align.channel.SlowWave = 'EEG M1G';
config{6}.align.channel.SlowWave_begin = 'EEG M1G';
config{6}.align.method.SlowWave     = 'nearestmin';
config{6}.align.method.SlowWave_begin = 'nearestmin';
config{6}.LFP.flip                  = 'no';

config{6}.LFP.electrodetoplot       = {'EEG M1G'};
config{6}.LFP.motorcortex           = {'EEG M1G'};
config{6}.plotseizure.h             = 5000;
config{6}.morpho.toiac              = [-0.5 0.3]; 
config{6}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{7}                           = configcommon;
config{7}.domorpho                  = 'yes';%no because filtered at 1Hz
config{7}.prefix                    = 'DTX32-';
config{7}.rawdir                    = fullfile(rootpath_data, 'DTX32_2019_08_01-INTRA', 'Brainvision');
config{7}.imagesavedir              = fullfile(imagesavedir,'DTX32');       % where to print images
config{7}.directorylist{1}          = []; %automatically set at the end of this script

config{7}.labels.macro              = {'ECoG-M1G'};
config{7}.injectiontime             = datetime('01-Aug-2019 13:10:06');
config{7}.LFP.channel               = config{7}.labels.macro;

config{7}.align.channel.SlowWave = 'ECoG-M1G';
config{7}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{7}.align.method.SlowWave     = 'nearestmin';
config{7}.align.method.SlowWave_begin = 'nearestmin';
config{7}.LFP.flip                  = 'no';

config{7}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{7}.LFP.motorcortex           = {'ECoG-M1G'};
config{7}.plotseizure.h             = 5000;
config{7}.morpho.toiac              = [-0.5 0.5]; 
config{7}.morpho.toibl              = [-2 -0.5];

%% Rodent 8
config{8}                           = configcommon;
config{8}.domorpho                  = 'yes';%no because filtered at 1Hz
config{8}.prefix                    = 'DTX33-';
config{8}.rawdir                    = fullfile(rootpath_data, 'DTX33_2019_08_06-INTRA', 'Brainvision');
config{8}.imagesavedir              = fullfile(imagesavedir,'DTX33');       % where to print images
config{8}.directorylist{1}          = []; %automatically set at the end of this script

config{8}.labels.macro              = {'ECoG-M1G'};
config{8}.injectiontime             = datetime('06-Aug-2019 13:19:00');
config{8}.LFP.channel               = config{8}.labels.macro;

config{8}.align.channel.SlowWave = 'ECoG-M1G';
config{8}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{8}.align.method.SlowWave     = 'nearestmin';
config{8}.align.method.SlowWave_begin = 'nearestmin';
config{8}.LFP.flip                  = 'no';

config{8}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{8}.LFP.motorcortex           = {'ECoG-M1G'};
config{8}.plotseizure.h             = 5000;
config{8}.morpho.toiac              = [-0.5 0.5]; 
config{8}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{9}                           = configcommon;
config{9}.domorpho                  = 'yes';%no because filtered at 1Hz
config{9}.prefix                    = 'DTX34-';
config{9}.rawdir                    = fullfile(rootpath_data, 'DTX34_2019_08_12-INTRA', 'Brainvision');
config{9}.imagesavedir              = fullfile(imagesavedir,'DTX34');       % where to print images
config{9}.directorylist{1}          = []; %automatically set at the end of this script

config{9}.labels.macro              = {'ECoG-M1G'};
config{9}.injectiontime             = datetime('12-Aug-2019 12:28:2');
config{9}.LFP.channel               = config{9}.labels.macro;

config{9}.align.channel.SlowWave = 'ECoG-M1G';
config{9}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{9}.align.method.SlowWave     = 'nearestmin';
config{9}.align.method.SlowWave_begin = 'nearestmin';
config{9}.LFP.flip                  = 'no';

config{9}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{9}.LFP.motorcortex           = {'ECoG-M1G'};
config{9}.plotseizure.h             = 5000;
config{9}.morpho.toiac              = [-0.5 0.5]; 
config{9}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{10}                           = configcommon;
config{10}.domorpho                  = 'yes';%no because filtered at 1Hz
config{10}.prefix                    = 'DTX35-';
config{10}.rawdir                    = fullfile(rootpath_data, 'DTX35_2019_08_14-INTRA', 'Brainvision');
config{10}.imagesavedir              = fullfile(imagesavedir,'DTX35');       % where to print images
config{10}.directorylist{1}          = []; %automatically set at the end of this script

config{10}.labels.macro              = {'ECoG-M1G'};
config{10}.injectiontime             = datetime('14-Aug-2019 12:43:37');
config{10}.LFP.channel               = config{10}.labels.macro;

config{10}.align.channel.SlowWave = 'ECoG-M1G';
config{10}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{10}.align.method.SlowWave     = 'nearestmin';
config{10}.align.method.SlowWave_begin = 'nearestmin';
config{10}.LFP.flip                  = 'no';

config{10}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{10}.LFP.motorcortex           = {'ECoG-M1G'};
config{10}.plotseizure.h             = 5000;
config{10}.morpho.toiac              = [-0.5 0.5]; 
config{10}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{11}                           = configcommon;
config{11}.domorpho                  = 'yes';%no because filtered at 1Hz
config{11}.prefix                    = 'DTX36-';
config{11}.rawdir                    = fullfile(rootpath_data, 'DTX36_2019_08_21-INTRA', 'Brainvision');
config{11}.imagesavedir              = fullfile(imagesavedir,'DTX36');       % where to print images
config{11}.directorylist{1}          = []; %automatically set at the end of this script

config{11}.labels.macro              = {'ECoG-M1G'};
config{11}.injectiontime             = datetime('21-Aug-2019 12:15:16');
config{11}.LFP.channel               = config{11}.labels.macro;

config{11}.align.channel.SlowWave = 'ECoG-M1G';
config{11}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{11}.align.method.SlowWave     = 'nearestmin';
config{11}.align.method.SlowWave_begin = 'nearestmin';
config{11}.LFP.flip                  = 'no';

config{11}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{11}.LFP.motorcortex           = {'ECoG-M1G'};
config{11}.plotseizure.h             = 5000;
config{11}.morpho.toiac              = [-0.5 0.5]; 
config{11}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
%fichier _002 data manquantes entre 6300 et 7100
config{12}                           = configcommon;
config{12}.domorpho                  = 'yes';%no because filtered at 1Hz
config{12}.prefix                    = 'DTX38-';
config{12}.rawdir                    = fullfile(rootpath_data, 'DTX38_2019_09_02-INTRA', 'Brainvision');
config{12}.imagesavedir              = fullfile(imagesavedir,'DTX38');       % where to print images
config{12}.directorylist{1}          = []; %automatically set at the end of this script

config{12}.labels.macro              = {'ECoG-M1G'};
config{12}.injectiontime             = datetime('02-Sep-2019 12:56:18');
config{12}.LFP.channel               = config{12}.labels.macro;

config{12}.align.channel.SlowWave = 'ECoG-M1G';
config{12}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{12}.align.method.SlowWave     = 'nearestmin';
config{12}.align.method.SlowWave_begin = 'nearestmin';
config{12}.LFP.flip                  = 'no';

config{12}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{12}.LFP.motorcortex           = {'ECoG-M1G'};
config{12}.plotseizure.h             = 5000;
config{12}.morpho.toiac              = [-0.5 0.5]; 
config{12}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{13}                           = configcommon;
config{13}.domorpho                  = 'yes';%no because filtered at 1Hz
config{13}.prefix                    = 'DTX49-';
config{13}.rawdir                    = fullfile(rootpath_data, 'DTX49_2020_08_10-INTRA-DTX', 'Brainvision');
config{13}.imagesavedir              = fullfile(imagesavedir,'DTX49');       % where to print images
config{13}.directorylist{1}          = []; %automatically set at the end of this script

config{13}.labels.macro              = {'ECoG-M1G'};
config{13}.injectiontime             = datetime('10-Aug-2020 12:30:36');
config{13}.LFP.channel               = config{13}.labels.macro;

config{13}.align.channel.SlowWave = 'ECoG-M1G';
config{13}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{13}.align.method.SlowWave     = 'nearestmin';
config{13}.align.method.SlowWave_begin = 'nearestmin';
config{13}.LFP.flip                  = 'no';

config{13}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{13}.LFP.motorcortex           = {'ECoG-M1G'};
config{13}.plotseizure.h             = 5000;
config{13}.morpho.toiac              = [-0.5 0.5]; 
config{13}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
config{14}                           = configcommon;
config{14}.domorpho                  = 'yes';%no because filtered at 1Hz
config{14}.prefix                    = 'DTX50-';
config{14}.rawdir                    = fullfile(rootpath_data, 'DTX50_2020_08_18-INTRA-DTX', 'Brainvision');
config{14}.imagesavedir              = fullfile(imagesavedir,'DTX50');       % where to print images
config{14}.directorylist{1}          = []; %automatically set at the end of this script

config{14}.labels.macro              = {'ECoG-M1G'};
config{14}.injectiontime             = datetime('18-Aug-2020 13:18:26');
config{14}.LFP.channel               = config{14}.labels.macro;

config{14}.align.channel.SlowWave = 'ECoG-M1G';
config{14}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{14}.align.method.SlowWave     = 'nearestmin';
config{14}.align.method.SlowWave_begin = 'nearestmin';
config{14}.LFP.flip                  = 'no';

config{14}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{14}.LFP.motorcortex           = {'ECoG-M1G'};
config{14}.plotseizure.h             = 5000;
config{14}.morpho.toiac              = [-0.5 0.5]; 
config{14}.morpho.toibl              = [-2 -0.5];

%% Rodent 7
%fichier très artefacté. Pour la morpho de l'onde lente, rétrécir la
%période pour détecter si artefacts
config{15}                           = configcommon;
config{15}.domorpho                  = 'yes';%no because filtered at 1Hz
config{15}.prefix                    = 'DTX52-';
config{15}.rawdir                    = fullfile(rootpath_data, 'DTX52-2020_08_18-INTRA-DTX', 'Brainvision');
config{15}.imagesavedir              = fullfile(imagesavedir,'DTX52');       % where to print images
config{15}.directorylist{1}          = []; %automatically set at the end of this script

config{15}.labels.macro              = {'ECoG-M1G'};
config{15}.injectiontime             = datetime('27-Aug-2020 12:02:18');
config{15}.LFP.channel               = config{15}.labels.macro;

config{15}.align.channel.SlowWave = 'ECoG-M1G';
config{15}.align.channel.SlowWave_begin = 'ECoG-M1G';
config{15}.align.method.SlowWave     = 'nearestmin';
config{15}.align.method.SlowWave_begin = 'nearestmin';
config{15}.LFP.flip                  = 'no';

config{15}.LFP.electrodetoplot       = {'ECoG-M1G'};
config{15}.LFP.motorcortex           = {'ECoG-M1G'};
config{15}.plotseizure.h             = 5000;
config{15}.morpho.toiac              = [-1 1]; 
config{15}.morpho.toibl              = [-2 -1];


%% find data files
for irat = 6:size(config,2)
    %config{irat}.directorylist
    filelist = dir(config{irat}.rawdir);
    i=0;
    for ifile = 1:length(filelist)
        [~,~,file_extension] = fileparts(filelist(ifile).name);
        if strncmp(file_extension,'.eeg',4)
            i=i+1;
            config{irat}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
        end
    end
    clear filelist
end

return 

% %% Rodent 16
% config{16}                           = configcommon;
% config{16}.domorpho                  = 'no';%no because filtered at 1Hz
% config{16}.prefix                    = 'DTX20-';
% config{16}.rawdir                    = fullfile(rootpath_data, 'DTX20-M1-10uM-INTRA', 'Brainvision');
% config{16}.imagesavedir              = fullfile(imagesavedir,'DTX20');       % where to print images
% config{16}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{16}.labels.macro              = {'EEG M1G'};
% config{16}.injectiontime             = datetime('21-May-2019 15:11:20');
% config{16}.LFP.channel               = config{16}.labels.macro;
% 
% %% Rodent 17
% config{17}                           = configcommon;
% config{17}.domorpho                  = 'no';%no because filtered at 1Hz
% config{17}.prefix                    = 'DTX23-';
% config{17}.rawdir                    = fullfile(rootpath_data, 'DTX23-M1-10uM-INTRA', 'Brainvision');
% config{17}.imagesavedir              = fullfile(imagesavedir,'DTX23');       % where to print images
% config{17}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{17}.labels.macro              = {'EEG M1G'};
% config{17}.injectiontime             = datetime('29-May-2019 12:44:56');
% config{17}.LFP.channel               = config{17}.labels.macro;
% 
% %% Rodent 18
% config{18}                           = configcommon;
% config{18}.domorpho                  = 'no';%no because filtered at 1Hz
% config{18}.prefix                    = 'DTX24-';
% config{18}.rawdir                    = fullfile(rootpath_data, 'DTX24-M1-10uM-INTRA', 'Brainvision');
% config{18}.imagesavedir              = fullfile(imagesavedir,'DTX24');       % where to print images
% config{18}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{18}.labels.macro              = {'EEG M1G'};
% config{18}.injectiontime             = datetime('03-Jun-2019 14:44:20');
% config{18}.LFP.channel               = config{18}.labels.macro;
% 
% %% Rodent 19
% config{19}                           = configcommon;
% config{19}.domorpho                  = 'no';%no because filtered at 1Hz
% config{19}.prefix                    = 'DTX26-';
% config{19}.rawdir                    = fullfile(rootpath_data, 'DTX26-M1-INTRA', 'Brainvision');
% config{19}.imagesavedir              = fullfile(imagesavedir,'DTX26');       % where to print images
% config{19}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{19}.labels.macro              = {'EEG M1G'};
% config{19}.injectiontime             = datetime('12-Jun-2019 12:15:18');
% config{19}.LFP.channel               = config{19}.labels.macro;
% 
% %% Rodent 20
% config{20}                           = configcommon;
% config{20}.domorpho                  = 'no';%no because filtered at 1Hz
% config{20}.prefix                    = 'DTX27-';
% config{20}.rawdir                    = fullfile(rootpath_data, 'DTX27-M1-INTRA', 'Brainvision');
% config{20}.imagesavedir              = fullfile(imagesavedir,'DTX27');       % where to print images
% config{20}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{20}.labels.macro              = {'EEG M1G'};
% config{20}.injectiontime             = datetime('13-Jun-2019 12:34:25');
% config{20}.LFP.channel               = config{20}.labels.macro;
% 
% %% Rodent 21
% config{21}                           = configcommon;
% config{21}.domorpho                  = 'no';%no because filtered at 1Hz
% config{21}.prefix                    = 'DTX28-';
% config{21}.rawdir                    = fullfile(rootpath_data, 'DTX28-M1-INTRA', 'Brainvision');
% config{21}.imagesavedir              = fullfile(imagesavedir,'DTX28');       % where to print images
% config{21}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{21}.labels.macro              = {'EEG M1G'};
% config{21}.injectiontime             = datetime('18-Jun-2019 13:03:01');
% config{21}.LFP.channel               = config{21}.labels.macro;
% 
% %% Rodent 22
% config{22}                           = configcommon;
% config{22}.domorpho                  = 'no';%no because filtered at 1Hz
% config{22}.prefix                    = 'DTX29-';
% config{22}.rawdir                    = fullfile(rootpath_data, 'DTX29-M1-INTRA', 'Brainvision');
% config{22}.imagesavedir              = fullfile(imagesavedir,'DTX29');       % where to print images
% config{22}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{22}.labels.macro              = {'EEG M1G'};
% config{22}.injectiontime             = datetime('21-Jun-2019 13:33:13');
% config{22}.LFP.channel               = config{22}.labels.macro;
% 
% %% Rodent 23
% config{23}                           = configcommon;
% config{23}.domorpho                  = 'no';%no because filtered at 1Hz
% config{23}.prefix                    = 'DTX30-';
% config{23}.rawdir                    = fullfile(rootpath_data, 'DTX30_2019_07_07-INTRA', 'Brainvision');
% config{23}.imagesavedir              = fullfile(imagesavedir,'DTX30');       % where to print images
% config{23}.directorylist{1}          = []; %automatically set at the end of this script
% 
% config{23}.labels.macro              = {'EEG M1G'};
% config{23}.injectiontime             = datetime('07-Jul-2019 14:47:56');
% config{23}.LFP.channel               = config{23}.labels.macro;



