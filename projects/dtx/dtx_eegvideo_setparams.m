function [config] = dtx_eegvideo_setparams(config)
% Setting parameters DTX project
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\natsort
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/natsort
end

disp('setting parameters for EEG-video');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-EEG-VIDEO/';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-EEGRodents-Brainvision/';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-EEG-VIDEO';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-EEGRodents-Brainvision\';
    os                  = 'windows';
else
    error('Platform not supported')
end


%% Congig common for all rats

datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis,'image');

configcommon.os                        = os;
configcommon.name                      = {'SlowWave','SlowWave_EMG_begin','Crise_End', 'Seizure'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.continous                 = true;
configcommon.group                     = 'dtx';

configcommon.missingdata = datetime.empty;%set for rats 4 and 7 for whom data is missing (one because unplug, the other because ampli was shut down on the first night of recording)

configcommon.muse.startmarker.SlowWave      = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave        = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave             = [-10, 5];
configcommon.epoch.pad.SlowWave             = 0;

configcommon.muse.startmarker.SlowWave_EMG_begin      = 'SlowWave_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.SlowWave_EMG_begin        = 'SlowWave_EMG__START__';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.SlowWave_EMG_begin             = [-10, 5];
configcommon.epoch.pad.SlowWave_EMG_begin             = 0;

configcommon.muse.startmarker.Crise_End      = 'Crise_End';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.Crise_End        = 'Crise_End';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.Crise_End             = [0, 5];
configcommon.epoch.pad.Crise_End             = 5;

configcommon.muse.startmarker.Seizure      = 'SlowWave';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.Seizure        = 'Crise_End';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.Seizure             = [-10, 10];
configcommon.epoch.pad.Seizure             = 0;

configcommon.minbadtime.SlowWave_EMG_begin  = 0;
configcommon.minbadtime.SlowWave            = 0;
configcommon.minbadtime.SlowWave_begin      = 0;
configcommon.minbadtime.Crise_End           = 0;
configcommon.minbadtime.SlowWave_not_aligned=0;
configcommon.minbadtime.Seizure             = 1;

configcommon.seizuretimings.marker_start = 'Crise_Start';
configcommon.seizuretimings.marker_end   = 'Crise_End';
configcommon.seizuretimings.analysis_start = 'Analysis_Start';
configcommon.seizuretimings.analysis_end   = 'Analysis_End';
configcommon.seizuretimings.winsize        = 3600;%s
configcommon.seizuretimings.winstep        = 1200;%s

configcommon.emgtimings.emg_start      = 'SlowWave_EMG__START__';
configcommon.emgtimings.emg_end        = 'SlowWave_EMG__END__';
configcommon.emgtimings.eeg_start      = 'SlowWave_begin';
configcommon.emgtimings.analysis_start =  'Analysis_Start';
configcommon.emgtimings.analysis_end   = 'Analysis_End';
configcommon.emgtimings.winsize        = 3600;%s
configcommon.emgtimings.winstep        = 1200;%s

configcommon.LFP.name                  = {'SlowWave','SlowWave_EMG_begin','SlowWave_begin','SlowWave_not_aligned', 'Seizure'};
configcommon.labels.macro              = {'M1G','M1D','PtA'};%do not put the emg channels here
configcommon.LFP.channel               = {};%set for each rat %do not put the emg channels here
configcommon.LFP.electrodetoplot       = {'PtA','M1D','M1G'};
configcommon.LFP.motorcortex           = {'M1G'};

configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 1;
configcommon.LFP.hpfiltord             = 4;
configcommon.LFP.hpfilttype            = 'but';
configcommon.LFP.resamplefs            = 512;
configcommon.LFP.reref                 = 'no';
configcommon.LFP.bsfilter              = 'yes';
configcommon.LFP.bsfreq                = [49 51];
configcommon.LFP.lpfilter              = 'no';
configcommon.LFP.lpfreq                = 30;
configcommon.LFP.lpfilttype            = 'fir';
configcommon.LFP.keepcfg               = 'no';
configcommon.LFP.write                 = false; %do not save after readLFP but after having removed artefacts, flipped, and corrected baseline

configcommon.LFP.baseline                          = 'yes';
configcommon.LFP.baselinewindow.SlowWave           = [-2 -1];
configcommon.LFP.baselinewindow.SlowWave_EMG_begin = [-2 -1];
configcommon.LFP.baselinewindow.Crise_End          = [3 5];
configcommon.LFP.baselinewindow.Seizure = [-10 -5];

configcommon.EMG.SlowWave              = {'EMG1', 'EMG2'};%'EMG1';%name of EMG channel associated with marker LFP.name. 'no' if no EMG associated to this seizure side 
configcommon.EMG.SlowWave_EMG_begin    = {'EMG1', 'EMG2'};
configcommon.EMG.Crise_End             = [];
configcommon.EMG.Seizure               = {'EMG1', 'EMG2'};
configcommon.EMG.reref                 = 'no';
configcommon.EMG.rerefmethod           = 'bipolar';
configcommon.EMG.refchannel            = 'EMG2';
configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];
configcommon.EMG.bsfiltord             = 3;
configcommon.EMG.envmethod             = 'rms';
configcommon.EMG.envparam              = 50;
configcommon.EMG.toi                   = [-5 5];

configcommon.align.name                = {'SlowWave','SlowWave_begin'};
configcommon.align.reref               = 'no';
configcommon.align.channel.SlowWave             = 'M1G';
configcommon.align.method.SlowWave              = [];%SET FOR EACH RAT 'nearestmin';      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter.SlowWave              = 'bp';
configcommon.align.freq.SlowWave                = [1, 2];          
configcommon.align.demean.SlowWave              = 'yes';
configcommon.align.thresh.value.SlowWave        = 0;
configcommon.align.thresh.method.SlowWave       = 'trial';%'medianbl','both';
configcommon.align.maxtimeshift.SlowWave        = 0.5;
configcommon.align.toiplot.SlowWave             = [-1,  1]; 
configcommon.align.toiactive.SlowWave           = [-0.5, 0.5];  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline.SlowWave         = [-10, -1];   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

configcommon.morpho.channame           = []; 
configcommon.morpho.negpeak            = 'yes';
configcommon.morpho.toiplot            = [-2 2];
configcommon.morpho.measurehalfwidth   = 'yes';
configcommon.morpho.measureamplitude   = 'yes';
configcommon.morpho.blmethod           = 'bl';
configcommon.morpho.toiac              = []; %set for each rat
configcommon.morpho.toibl              = []; %set for each rat
configcommon.morpho.raw_facealpha      = 0.1;
configcommon.morpho.plotavg            = 'no';

configcommon.TFR.foi                   = 3:2:200;
configcommon.TFR.tapsmofrq             = 5.02;
configcommon.TFR.timewinsize           = 0.35;

%% Rodent 1
config{1}                           = configcommon;
config{1}.prefix                    = 'Rat-2020_02_19-1-';
config{1}.rawdir                    = fullfile(rootpath_data,'Rat-2020_02_19-1');
config{1}.rawlabels.oldnames        = {'47','45','48','46','49'}; %for conversion deltamed to brainvision
config{1}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{1}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_02_19-1');       % where to print images

config{1}.injectiontime             = datetime('19-Feb-2020 10:30:00');
config{1}.align.method.SlowWave     = 'nearestmin';
config{1}.align.method.SlowWave_begin = 'nearestmin';
config{1}.LFP.channel               = {'M1G','M1D','PtA'};%set for each rat %do not put the emg channels here
config{1}.LFP.flip                  = 'no';

config{1}.plotseizure.h             = 5000;
config{1}.morpho.toiac              = [-0.5 1]; 
config{1}.morpho.toibl              = [-2 -0.5];

%% Rodent 2 

config{2}                           = configcommon;
config{2}.prefix                    = 'Rat-2020_02_19-2-';
config{2}.rawdir                    = fullfile(rootpath_data,'Rat-2020_02_19-2');
config{2}.rawlabels.oldnames        = {'57','55','56','59'}; %for conversion deltamed to brainvision
config{2}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{2}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_02_19-2'); % where to print images

config{2}.injectiontime             = datetime('19-Feb-2020 15:20:00');
config{2}.align.method.SlowWave     = 'nearestmin';
config{2}.align.method.SlowWave_begin     = 'nearestmin';
config{2}.LFP.channel               = {'M1G','M1D','PtA'};%set for each rat %do not put the emg channels here
config{2}.LFP.flip                  = 'no';
config{2}.plotseizure.h             = 2500;
config{2}.morpho.toiac              = [-0.4 0.6]; 
config{2}.morpho.toibl              = [-2 -0.4];

%% Rodent 3

config{3}                           = configcommon;
config{3}.prefix                    = 'Rat-2020_06_03-1-';
config{3}.rawdir                    = fullfile(rootpath_data,'Rat-2020_06_03-1');
config{3}.rawlabels.oldnames        = {'52','50','51','54'}; %for conversion deltamed to brainvision
config{3}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{3}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_06_03-1');       % where to print images

config{3}.injectiontime             = datetime('03-Jun-2020 16:15:00');
config{3}.align.method.SlowWave     = 'nearestmin';
config{3}.align.method.SlowWave_begin  = 'nearestmin';
config{3}.LFP.channel               = {'M1G','M1D'};%set for each rat %do not put the emg channels here
config{3}.LFP.flip                  = 'no';
config{3}.plotseizure.h             = 20000;
config{3}.morpho.toiac              = [-0.9 1.2]; 
config{3}.morpho.toibl              = [-2 -0.9];

%% Rodent 4 
% début à 17h00, 5h30 post injection, car changement de position à ce moment
% là. 

config{4}                           = configcommon;
config{4}.prefix                    = 'Rat-2020_06_09-1-';
config{4}.rawdir                    = fullfile(rootpath_data,'Rat-2020_06_09-1');
config{4}.rawlabels.oldnames        = {'52','50','51','54'}; %for conversion deltamed to brainvision
config{4}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{4}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_06_09-1');       % where to print images

config{4}.injectiontime             = datetime('09-Jun-2020 11:30:00');
config{4}.missingdata(1)            = datetime('09-Jun-2020 13:26:44'); %fichier 2 perdu
config{4}.missingdata(2)            = datetime('09-Jun-2020 14:26:47');
config{4}.align.method.SlowWave     = 'nearestmin';
config{4}.align.method.SlowWave_begin     = 'nearestmin';
config{4}.LFP.channel               = {'M1G','M1D'};%set for each rat %do not put the emg channels here
config{4}.LFP.flip                  = 'no';
config{4}.plotseizure.h             = 10000;
config{4}.morpho.toiac              = [-0.4 0.5]; 
config{4}.morpho.toibl              = [-2 -0.5];

%% Rodent 5

% l'électrode M1G lache pendant le fichier 2020_06_30_19-37. On voit quand
% même les crises sur l'électrode M1D
config{5}                           = configcommon;
config{5}.prefix                    = 'Rat-2020_06_29-1-';
config{5}.rawdir                    = fullfile(rootpath_data,'Rat-2020_06_29-1');
config{5}.rawlabels.oldnames        = {'1-S1L','1-S1R','1-M1L','1-NTS'}; %for conversion deltamed to brainvision
config{5}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{5}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_06_29-1');       % where to print images

config{5}.injectiontime             = datetime('29-Jun-2020 11:15:00');
config{5}.align.method.SlowWave     = 'nearestmax';
config{5}.align.method.SlowWave_begin     = 'nearestmax';
config{5}.LFP.channel               = {'M1G'};%set for each rat %do not put the emg channels here
config{5}.LFP.flip                  = 'yes';
config{5}.plotseizure.h             = 2000;
config{5}.morpho.toiac              = [-0.3 0.5];
config{5}.morpho.toiac              = [-0.5 0.5]; 
config{5}.morpho.toibl              = [-2 -0.5];

%% Rodent 

config{6}                           = configcommon;
config{6}.prefix                    = 'Rat-2020_06_29-2-';
config{6}.rawdir                    = fullfile(rootpath_data,'Rat-2020_06_29-2');
config{6}.rawlabels.oldnames        = {'2-S1L','2-S1R','2-M1L','2-NTS'}; %for conversion deltamed to brainvision
config{6}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{6}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_06_29-2');       % where to print images

config{6}.injectiontime             = datetime('29-Jun-2020 14:55:00');
config{6}.EMG.reref                 = 'no';
config{6}.EMG.rerefmethod           = [];
config{6}.EMG.refchannel            = [];
config{6}.LFP.channel               = {'M1G','M1D'};%set for each rat %do not put the emg channels here
config{6}.LFP.flip                  = 'yes';
config{6}.align.method.SlowWave     = 'nearestmax';
config{6}.align.method.SlowWave_begin     = 'nearestmax';
config{6}.plotseizure.h             = 2000;
config{6}.morpho.toiac              = [-1 0.9]; 
config{6}.morpho.toibl              = [-2 -1];;

%% Rodent 7

config{7}                           = configcommon;
config{7}.prefix                    = 'Rat-2020_09_24-1-';
config{7}.rawdir                    = fullfile(rootpath_data,'Rat-2020_09_24-1');
config{7}.rawlabels.oldnames        = {'1-S1L','1-S1R','1-M1L','1-NTS'}; %for conversion deltamed to brainvision
config{7}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{7}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_09_24-1');       % where to print images

config{7}.injectiontime             = datetime('23-Sep-2020 17:39:00');
config{7}.missingdata(1)            = datetime('23-Sep-2020 20:38:36');%ampli coupé dans la nuit
config{7}.missingdata(2)            = datetime('24-Sep-2020 10:11:27');
config{7}.EMG.reref                 = 'no';
config{7}.EMG.rerefmethod           = [];
config{7}.EMG.refchannel            = [];
config{7}.LFP.channel               = {'M1G'};%set for each rat %do not put the emg channels here
config{7}.LFP.flip                  = 'yes';
config{7}.align.method.SlowWave     = 'nearestmax';
config{7}.align.method.SlowWave_begin     = 'nearestmax';
config{7}.plotseizure.h             = 2000;
config{7}.morpho.toiac              = [-0.5 1]; 
config{7}.morpho.toibl              = [-2 -0.5];

%% Rodent 8

config{8}                           = configcommon;
config{8}.prefix                    = 'Rat-2020_10_09-1-';
config{8}.rawdir                    = fullfile(rootpath_data,'Rat-2020_10_09-1');
config{8}.rawlabels.oldnames        = {'2-S1L','2-S1R','2-M1L','2-NTS'};  %for conversion deltamed to brainvision
config{8}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{8}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_10_09-1');       % where to print images

config{8}.injectiontime             = datetime('09-Oct-2020 11:36:00');
config{8}.LFP.channel               = {'M1G','M1D'};%set for each rat %do not put the emg channels here
config{8}.LFP.flip                  = 'yes';
config{8}.EMG.hpfreq                = 50;
config{8}.align.method.SlowWave     = 'nearestmax';
config{8}.align.method.SlowWave_begin     = 'nearestmax';
config{8}.plotseizure.h             = 2000;
config{8}.morpho.toiac              = [-0.5 0.7]; 
config{8}.morpho.toibl              = [-2 -0.5];

%% Rodent 9

config{9}                           = configcommon;
config{9}.prefix                    = 'Rat-2020_09_24-2-';
config{9}.rawdir                    = fullfile(rootpath_data,'Rat-2020_09_24-2');
config{9}.rawlabels.oldnames        = {'2-S1L','2-S1R','2-M1L','2-NTS'};  %for conversion deltamed to brainvision
config{9}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{9}.imagesavedir              = fullfile(imagesavedir,'Rat-2020_09_24-2');       % where to print images

config{9}.injectiontime             = datetime('23-Sep-2020 10:54:00');
config{9}.LFP.flip                  = 'yes';
config{9}.LFP.channel               = {'M1G','M1D'};%set for each rat %do not put the emg channels here
config{9}.align.method.SlowWave     = 'nearestmax';
config{9}.align.method.SlowWave_begin     = 'nearestmax';
config{9}.plotseizure.h             = 2000;
config{9}.morpho.toiac              = [-0.5 0.6]; 
config{9}.morpho.toibl              = [-2 -0.5];

%% find files
for irat = 1:size(config,2)
    filelist = dir(config{irat}.rawdir);
    filelist = natsort({filelist.name});
    i=0;
    for ifile = 1:length(filelist)
        [~,~,file_extension] = fileparts(filelist{ifile});
        if strncmp(file_extension,'.eeg',4)
            i=i+1;
            config{irat}.directorylist{1}{i} =  filelist{ifile}(1:end-4);
        end
    end
    clear filelist
end