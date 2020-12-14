
% Setting parameters DTX project Paul Baudin

function [config] = dtx_spikes_setparams(config)

disp('setting parameters for probe spike analysis');
if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-PROBE\';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data'); %removed of config{i} so we can more easily modify it
imagesavedir = fullfile(rootpath_analysis, 'image_spike');

%% Config common for all experiments
configcommon.os                        = os;
configcommon.datasavedir               = datasavedir;
configcommon.muse.backupdir            = fullfile(rootpath_analysis,'Musemarkers_backup');
configcommon.hyplabels                 = [];

configcommon.unit_table = fullfile(rootpath_analysis,'classification_units.xlsx');
configcommon.commonchans               = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

configcommon.LFP.flip                  = 'yes';
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.resamplefs            = 320; %because sampling rate is 3200Hz
configcommon.LFP.baseline              = 'yes';
configcommon.LFP.write                 = true;
configcommon.LFP.keepcfg           	   = false;

%write spyking circus
configcommon.circus.version            = 'fieldtrip';
configcommon.circus.reref              = 'no';
configcommon.circus.outputdir          = 'SpykingCircus';
configcommon.circus.paramfile          = fullfile(rootpath_analysis,'dtx_spikes_SpykingCircus.params');
configcommon.circus.params.detection.spike_thresh  = '8';
configcommon.circus.params.filtering.cut_off       = '500, auto';
configcommon.circus.params.clustering.max_elts     = '30000';
configcommon.circus.params.clustering.nb_repeats   = '15';
configcommon.circus.params.clustering.merging_method= 'distance';
configcommon.circus.params.clustering.merging_param= '2';
configcommon.circus.params.detection.peaks         = 'positive';
configcommon.circus.params.extracting.max_elts     = '30000';
%configcommon.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
%configcommon.circus.postfix            = [];%'-final';
%configcommon.circus.deadfilesuffix     = 'withseizures';
configcommon.circus.channel = {}; %set for each rat
configcommon.circus.channelname = {}; %set for each rat, analyse each channel independently

configcommon.spike.nrsdfbins           = 200;
configcommon.spike.RPV                 = 0.001; %refractory period violation, in seconds
configcommon.spike.ISIbins             = [0:0.0005:0.050]; %in s
configcommon.spikewin.windowsize       = 60;
configcommon.spikewin.windowoverlap    = 0;

configcommon.stats.alpha = 0.05;

%spike waveform
configcommon.spikewaveform.toi         = [-0.0015 0.0015]; %in s
configcommon.spikewaveform.cutoff      = 500; %high pass filter frequency to apply to raw data
configcommon.spikewaveform.nspikes     = 'all'; %maximum number of spike waveforms to load. Can be 'all'. 

%% Config for dtx experiments
configdtx = configcommon;
configdtx.minbadtime.window = 1;
configdtx.minbadtime.SlowWave = 0;
% configdtx.minbadtime.Seizure  = 0;
configdtx.minbadtime.Interictal = 1;

configdtx.type                      = 'dtx';
configdtx.name                      = {'SlowWave', 'Interictal'};%,'Seizure','Interictal'};%,'SlowWave_Larger','xcorr_10_1', 'xorr_2_1', 'xcorr_1_0', 'xcorr_1_1','Baseline'};
configdtx.spike.events_name         = {'SlowWave'};%,'Seizure'};%, 'Seizure'};
configdtx.spike.baseline_name       = 'Interictal';

configdtx.seizuretimings.marker_start     = "SlowWave";
configdtx.seizuretimings.marker_end       = "SlowWave";
configdtx.seizuretimings.analysis_start = 'Analysis_Start';
configdtx.seizuretimings.analysis_end   = 'Analysis_End';
configdtx.seizuretimings.winsize        = 3600;%s
configdtx.seizuretimings.winstep        = 1200;%s

%read LFP
configdtx.LFP.name                  = {'SlowWave','Interictal'};%_EEG', 'SLowWave_Intra'};

%add bad markers for spiking cricus
configdtx.bad.markerStart           = 'Crise_Start';
configdtx.bad.markerEnd             = 'Crise_End';
configdtx.bad.time_from_begin       = -3;
configdtx.bad.time_from_end         = 3;    

%remove trials which intersect BAD markers
configdtx.rmtrials.plotdata         = 'yes';
configdtx.rmtrials.write            = 'no';
configdtx.rmtrials.electrodetoplot  = [];%set separately for each rat. Only for dtx rats

%SlowWave
configdtx.muse.startmarker.SlowWave      = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.muse.endmarker.SlowWave        = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.epoch.toi.SlowWave             = [-2, 2];
configdtx.epoch.pad.SlowWave             = 2;
configdtx.spike.toi.SlowWave             = [-2, 2];

% ALIGN TO BEGIN
configdtx.align.name                = {'SlowWave'};
configdtx.align.method.SlowWave     = 'max';                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configdtx.align.filter.SlowWave              = 'lp';
configdtx.align.freq.SlowWave                = 10;                                                                                  % lowpass filter freq to smooth peak detection (Hz)
configdtx.align.thresh.value.SlowWave        = 1;
configdtx.align.thresh.method.SlowWave       = 'trial';%,'trial','trial'};%'medianbl','both';
configdtx.align.toiplot.SlowWave             = [-1,  1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configdtx.align.toiactive.SlowWave           = [-0.5, 0.5];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configdtx.align.toibaseline.SlowWave         = [-1, -0.5];
configdtx.align.maxtimeshift.SlowWave        = 0.3;
configdtx.align.demean.SlowWave              = 'yes';
%for detection of begin of event
configdtx.align.findbegin.SlowWave          = 'yes';
configdtx.align.beginthresh.SlowWave        = 0.1; % percent of peak %FIXME � v�rifier

configdtx.LFP.baselinewindow.SlowWave       = [-2, -1];
configdtx.stats.bl.SlowWave            = [-2, -1];
% configdtx.stats.ac.SlowWave            = [-1, 1];
configdtx.spike.resamplefs.SlowWave    = 1000;
configdtx.spike.psthbin.SlowWave       = 1/50; %in s
configdtx.spike.sdftimwin.SlowWave     = [-0.05 0.05];
    
% %Interictal
configdtx.muse.startmarker.Interictal      = 'Crise_End';   % start and end Muse marker. For defining trials
configdtx.muse.endmarker.Interictal        = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.epoch.toi.Interictal             = [2, -1];
configdtx.epoch.pad.Interictal             = 2;
configdtx.spike.toi.Interictal             = [2, -1];
configdtx.LFP.baselinewindow.Interictal    = [2, 12];
configdtx.stats.bl.Interictal              = [];%set automatically. all but the last 60 seconds
% configdtx.stats.ac.Interictal              = [12, Inf];
configdtx.spike.resamplefs.Interictal      = 1;%1/size of bar graph bins
configdtx.spike.psthbin.Interictal         = 1; %in s
configdtx.spike.sdftimwin.Interictal       = [-5 5];


% %Seizure
% configdtx.muse.startmarker.Seizure      = 'SlowWave';   % start and end Muse marker. For defining trials
% configdtx.muse.endmarker.Seizure        = 'Crise_End';   % start and end Muse marker. For defining trials
% configdtx.epoch.toi.Seizure             = [-2, 2];
% configdtx.epoch.pad.Seizure             = 0;
% configdtx.spike.toi.Seizure             = [-2, 2];
% configdtx.LFP.baselinewindow.Seizure    = [-2, -1];
% configdtx.stats.bltoi.Seizure           = [-2, -1];
% configdtx.stats.actoi.Seizure           = [-1, Inf];
% configdtx.spike.resamplefs.Seizure      = 100;%1/size of bar graph bins
% configdtx.spike.psthbin.Seizure            = 1/50; %in s
% 

%compute spike stats
% configdtx.stats.numrandomization    = 1000; %Essayer � 1000
% configdtx.stats.alpha               = 0.025;

% %stats over time
% configdtx.statstime.timewin          = 10;
% configdtx.statstime.slidestep        = 2;
% configdtx.statstime.removeempty      = 'yes';

%plot spike quality only on interictal data
% configdtx.spikequal.label_list = 'Interictal';
% ALIGN XCORR : 
% configdtx.align.name                = {'SlowWave'};                               % Name of markers/patterns to align
% configdtx.align.channel.SlowWave             = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP'};    % Channels to use for alignment
% configdtx.align.demean              = 'yes';
% configdtx.align.baselinewindow      = [-2 -1.5];
% configdtx.align.reref               = 'no';
% configdtx.align.refmethod           = 'bipolar';
% configdtx.align.latency             = [-1, 1];  

% to smooth spikerate. Not used for now. Maybe for correlation LFP-spike
% configdtx.spike.toispikerate.SlowWave     = [-0.1 0.1];            % for plotting spikerate
% configdtx.spike.toispikerate.Seizure     = [-0.5 0.5];            % for plotting spikerate
% configdtx.spike.toispikerate.Interictal     = [-10 10];                % for plotting spikerate

%% Config for control experiments
% no seizure-time-locked trials
% markers : Baseline_Start, Analysis_End, Injection, BAD
configctrl = configcommon;

configctrl.type                      = 'ctrl';
configctrl.name                      = {'Control'};
configctrl.spike.baseline_name        = 'Control';
configctrl.spike.events_name          = {};
configctrl.minbadtime.window        = 1;

configctrl.LFP.name                  = []; %do not load LFP
configctrl.LFP.channel               = []; %do not load LFP
configctrl.align.name                = []; %do not align muse markers
configctrl.align.channel             = []; %do not align muse markers

configctrl.muse.write                = true;
configctrl.muse.startmarker.Control  = 'no';
configctrl.muse.endmarker.Control  = 'no';
% list of onset timing with respect to start-marker (s)
% configctrl.epoch.toi.Control              = [0, 0];
% configctrl.epoch.pad.Control              = 10; %for LFP

configctrl.circus.reref              = 'no';
configctrl.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'dtx', 'SpykingCircus');
configctrl.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
configctrl.circus.postfix            = [];%'-final';

configctrl.bad.markerStart           = [];
configctrl.bad.markerEnd             = [];

configctrl.spike.triallength         = 60; %seconds

% configctrl.spikequal.label_list = {'Control'};

% %remove trials which intersect BAD markers
% configctrl.rmtrials.plotdata         = 'no';
% configctrl.rmtrials.write            = 'no';

%% Rodent 1
config{1}                           = configdtx;
config{1}.prefix                    = 'DTX5-';
config{1}.rawdir                    = fullfile(rootpath_data, 'DTX5-M1-10uM', '2019_03_19_DTX-5');
config{1}.imagesavedir              = fullfile(imagesavedir, 'DTX5');       % where to print images
config{1}.directorylist{1}          =  {'2019-03-19_13-18',...
    '2019-03-19_14-30',...
    '2019-03-19_16-30',...
    '2019-03-19_18-30',...
    '2019-03-19_20-30',...
    '2019-03-19_22-30',...
    '2019-03-20_00-30'};

config{1}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{1}.labels.macro              = {'E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{1}.injectiontime             = datetime('19-Mar-2019 13:55:00');

config{1}.align.channel.SlowWave    = 'E12LFP';                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.LFP.channel               = config{1}.labels.macro;
config{1}.LFP.electrodetoplot.SlowWave = 'E12LFP';
config{1}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{1}.circus.channelname        = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{1}.rmtrials.electrodetoplot.SlowWave  = 'E12LFP';

%% Rodent 2
%ATTENTION : SUR MUSE, MAUVAIS NOMS DE CHANNELS? CORRIGES AVEC
%DTC_CORRECTDTX2NAME
config{2}                           = configdtx;
config{2}.prefix                    = 'DTX2-';
config{2}.rawdir                    = fullfile(rootpath_data, 'DTX2-M1-10uM', '2019_03_01_DTX-2');
config{2}.imagesavedir              = fullfile(imagesavedir, 'DTX2');       % where to print images
config{2}.directorylist{1}          =  {'2019-03-01_12-14',...
    '2019-03-01_12-33',...
    '2019-03-01_14-14',...
    '2019-03-01_16-14',...
    '2019-03-01_18-14',...
    '2019-03-01_20-14'};

config{2}.labels.micro              = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{2}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGS1','ECoGM1','ECoGPtA'}; % /!\ Mistake during acquisition : ECoGS1 = ECoGM1G, ECoGM1 = ECoGM1D. Corrected in the code after having loaded the data (need acquisitin name for laoding)

config{2}.injectiontime             = datetime('01-Mar-2019 12:33:00');

config{2}.align.channel.SlowWave             = 'E13LFP';%{'ECoGS1'};% %ATTENTION erreur nom eeg pendant acquisition
config{2}.LFP.channel               = config{2}.labels.macro;%ECoG S1 et ECoGM1 renomm�s respectivement ECoGM1G et ECoGM1D, apr�s l'�tape de readLFP
config{2}.LFP.electrodetoplot.SlowWave       = 'E13LFP'; %ATTENTION ECoGS1 est renomm� ECoGM1G pendant le script LFP

config{2}.circus.channel            = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{2}.circus.channelname        = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{2}.rmtrials.electrodetoplot.SlowWave  = 'E13LFP';


%% Rodent 3
config{3}                           = configdtx;
config{3}.prefix                    = 'DTX4-';
config{3}.rawdir                    = fullfile(rootpath_data, 'DTX4-M1-10uM', '2019_03_08_DTX-4');
config{3}.imagesavedir              = fullfile(imagesavedir,'DTX4');       % where to print images
config{3}.directorylist{1}          =  {'2019-03-08_12-28',...
    '2019-03-08_14-05',...
    '2019-03-08_16-05',...
    '2019-03-08_18-05',...
    '2019-03-08_20-05'};

config{3}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{3}.labels.macro              = {'E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP'...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{3}.injectiontime             = datetime('08-Mar-2019 13:21:00');

config{3}.align.channel.SlowWave             = 'E13LFP';%{'ECoGM1G'};%
config{3}.LFP.channel               = config{3}.labels.macro;
config{3}.LFP.electrodetoplot.SlowWave       = 'E13LFP';
config{3}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{3}.circus.channelname        = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{3}.rmtrials.electrodetoplot.SlowWave  = 'E13LFP';


%% Rodent 4
config{4}                           = configdtx;
config{4}.prefix                    = 'DTX7-';
config{4}.rawdir                    = fullfile(rootpath_data, 'DTX7-M1-10uM', '2019_03_22_DTX-7');
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX7');       % where to print images
config{4}.directorylist{1}          =  {'2019-03-22_12-31',...
    '2019-03-22_14-31',...
    '2019-03-22_16-31',...
    '2019-03-22_18-31'};

config{4}.labels.micro              = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15'};
config{4}.labels.macro              = {'E06LFP','E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{4}.injectiontime             = datetime('22-Mar-2019 13:40:07');

config{4}.align.channel.SlowWave             = 'E13LFP';
config{4}.LFP.channel               = config{4}.labels.macro;
config{4}.LFP.electrodetoplot.SlowWave       = 'E13LFP';
config{4}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15'};
config{4}.circus.channelname        = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15'};
config{4}.rmtrials.electrodetoplot.SlowWave  = 'E13LFP';


%% Rodent 5
%Remove the first 2 files and keep only the 2 lasts because the probe was moved in
%order to have better unit activity
config{5}                           = configdtx;
config{5}.prefix                    = 'DTX6-';
config{5}.rawdir                    = fullfile(rootpath_data, 'DTX6-M1-10uM', '2019_03_21_DTX-6');
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX6');       % where to print images
config{5}.directorylist{1}          =  {'2019-03-21_18-12', '2019-03-21_20-12'};

config{5}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{5}.labels.macro              = {'E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{5}.injectiontime             = datetime('21-Mar-2019 15:08:00');

config{5}.align.channel.SlowWave             = 'E13LFP';
config{5}.LFP.channel               = config{5}.labels.macro;
config{5}.LFP.electrodetoplot.SlowWave       = 'E13LFP';
config{5}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{5}.circus.channelname        = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{5}.rmtrials.electrodetoplot.SlowWave  = 'E13LFP';


%% Rodent 6
config{6}                           = configctrl;
config{6}.prefix                    = 'DTX43-';
config{6}.rawdir                    = fullfile(rootpath_data, 'DTX43_2020-05-22-PROBE-CONTROLE', '2020_05_22_DTX43-ctrl');
config{6}.imagesavedir              = fullfile(imagesavedir,'DTX43');       % where to print images
config{6}.directorylist{1}          =  {'2020-05-22_14-05',...
    '2020-05-22_16-05',...
    '2020-05-22_18-05',...
    '2020-05-22_20-05',...
    '2020-05-22_22-05'};

config{6}.labels.micro              = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{6}.labels.macro              = {'E06LFP','E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{6}.injectiontime             = datetime('22-May-2020 15:08:02');

config{6}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{6}.circus.channelname        = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};

%% Rodent 7
%oblig� de reref avant de lancer spyking circus, car sinon trop d'artefacts
%+++, et il y a un canal sans spike qui contient uniquement les artefacts.
config{7}                           = configctrl;
config{7}.prefix                    = 'DTX44-';
config{7}.rawdir                    = fullfile(rootpath_data, 'DTX44_2020-05-24-Probe-controle', '2020_05_24_DTX44-ctrl-probe');
config{7}.imagesavedir              = fullfile(imagesavedir,'DTX44');       % where to print images
config{7}.directorylist{1}          =  {'2020-05-24_13-10',...
    '2020-05-24_15-10',...
    '2020-05-24_17-10',...
    '2020-05-24_19-10',...
    '2020-05-24_21-10',...
    '2020-05-24_23-10',...
    '2020-05-25_01-10'};

config{7}.labels.micro              = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{7}.labels.macro              = {'E06LFP','E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};



config{7}.injectiontime             = datetime('24-May-2020 14:14:02');

config{7}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15'};
config{7}.circus.channelname        = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15'};
config{7}.circus.reref              = 'yes';
config{7}.circus.refchan            = 'E16';



end




