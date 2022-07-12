function [config] = dtx_spikes_setparams(config)

disp('setting parameters for probe spike analysis');
if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-spikes/';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/lgi1/DTX-raw/';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-spikes\';
    rootpath_data       = '\\lexport\iss01.charpier\raw\lgi1\DTX-raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data', 'reviewer_answers'); 
imagesavedir = fullfile(rootpath_analysis, 'image', 'reviewer_answers');

%% Config common for all experiments
configcommon.datasavedir               = datasavedir;
configcommon.muse.backupdir            = fullfile(rootpath_analysis,'Musemarkers_backup');

configcommon.unit_table = fullfile(rootpath_analysis,'classification_units.xlsx');
configcommon.commonchans               = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

configcommon.LFP.flip                  = 'yes';
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.resamplefs            = 320;
configcommon.LFP.baseline              = 'yes';
configcommon.LFP.write                 = true;
configcommon.LFP.keepcfg           	   = false;

%write spyking circus
configcommon.circus.version            = 'fieldtrip';
configcommon.circus.reref              = 'no';
configcommon.circus.outputdir          = 'SpykingCircus';
configcommon.circus.paramfile          = fullfile(rootpath_analysis,'dtx_spikes_SpykingCircus.params');
configcommon.circus.params.detection.spike_thresh  = '10';
configcommon.circus.params.filtering.cut_off       = '500, auto';
configcommon.circus.params.clustering.max_elts     = '30000';
configcommon.circus.params.clustering.nb_repeats   = '15';
configcommon.circus.params.clustering.merging_method= 'distance';
configcommon.circus.params.clustering.merging_param= '2';
configcommon.circus.params.detection.peaks         = 'positive';
configcommon.circus.params.extracting.max_elts     = '30000';
configcommon.circus.channel = {}; %set for each rat
configcommon.circus.channelname = {}; %set for each rat, analyse each channel independently

configcommon.spike.RPV                 = 0.002; %refractory period violation, in seconds
configcommon.spike.ISIbins             = [0:0.0005:0.050]; %in s

configcommon.stats.alpha = 0.05;

%spike waveform
configcommon.spikewaveform.toi         = [-0.0015 0.0015]; 
configcommon.spikewaveform.cutoff      = 500; 
configcommon.spikewaveform.nspikes     = 'all'; %maximum number of spike waveforms to load. Can be 'all'. 


%% config dtx
configdtx = configcommon;
configctrl.type                     = 'dtx';
configdtx.minbadtime.window         = 1;
configdtx.minbadtime.SlowWave       = 0;
configdtx.minbadtime.Interictal     = 1;

configdtx.type                      = 'dtx';
configdtx.name                      = {'SlowWave', 'Interictal'};
configdtx.spike.events_name         = {'SlowWave'};
configdtx.spike.baseline_name       = 'Interictal';

configdtx.seizuretimings.marker_start     = "SlowWave";
configdtx.seizuretimings.marker_end       = "SlowWave";
configdtx.seizuretimings.analysis_start = 'Analysis_Start';
configdtx.seizuretimings.analysis_end   = 'Analysis_End';
configdtx.seizuretimings.winsize        = 3600;%s
configdtx.seizuretimings.winstep        = 1200;%s

%read LFP
configdtx.LFP.name                  = {'SlowWave','Interictal'};

%add bad markers for spiking cricus
configdtx.bad.markerStart           = 'Crise_Start';
configdtx.bad.markerEnd             = 'Crise_End';
configdtx.bad.time_from_begin       = -3;
configdtx.bad.time_from_end         = 3;    

%SlowWave
configdtx.muse.startmarker.SlowWave      = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.muse.endmarker.SlowWave        = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.epoch.toi.SlowWave             = [-2, 2];
configdtx.epoch.pad.SlowWave             = 2;
configdtx.spike.toi.SlowWave             = [-2, 2];

configdtx.align.name                         = {'SlowWave'};
configdtx.align.method.SlowWave              = 'max';                                                              
configdtx.align.filter.SlowWave              = 'lp';
configdtx.align.freq.SlowWave                = 10;                                                                               
configdtx.align.thresh.value.SlowWave        = 1;
configdtx.align.thresh.method.SlowWave       = 'trial';%,'trial','trial'};%'medianbl','both';
configdtx.align.toiplot.SlowWave             = [-1,  1];                                           
configdtx.align.toiactive.SlowWave           = [-0.5, 0.5];                                           
configdtx.align.toibaseline.SlowWave         = [-1, -0.5];
configdtx.align.maxtimeshift.SlowWave        = 0.3;
configdtx.align.demean.SlowWave              = 'yes';

configdtx.LFP.baselinewindow.SlowWave  = [-2, -1];
configdtx.stats.bl.SlowWave            = [-2, -1];
configdtx.spike.resamplefs.SlowWave    = 1000;
configdtx.spike.psthbin.SlowWave       = 1/50; %s
configdtx.spike.sdftimwin.SlowWave     = [-0.05 0.05];
configdtx.spike.nrsdfbins.SlowWave     = 200;
    
%Interictal
configdtx.muse.startmarker.Interictal      = 'Crise_End';   % start and end Muse marker. For defining trials
configdtx.muse.endmarker.Interictal        = 'SlowWave';   % start and end Muse marker. For defining trials
configdtx.epoch.toi.Interictal             = [2, -1];
configdtx.epoch.pad.Interictal             = 2;
configdtx.spike.toi.Interictal             = [2, -1];
configdtx.LFP.baselinewindow.Interictal    = [2, 12];
configdtx.stats.bl.Interictal              = [];%set automatically. all but the last 60 seconds
configdtx.spike.resamplefs.Interictal      = 1;
configdtx.spike.psthbin.Interictal         = 1;%in s
configdtx.spike.sdftimwin.Interictal       = [-5 5];
configdtx.spike.nrsdfbins.Interictal       = 200;

configdtx.spikewin.windowsize       = 60;
configdtx.spikewin.windowoverlap    = 0;

%% Config for control experiments
% no seizure-time-locked trials
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
configctrl.circus.reref              = 'no';
configctrl.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'dtx', 'SpykingCircus');
configctrl.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
configctrl.circus.postfix            = [];%'-final';

configctrl.bad.markerStart           = [];
configctrl.bad.markerEnd             = [];

configctrl.spike.triallength         = 240; %seconds
configctrl.spikewin.windowsize       = 240;
configctrl.spikewin.windowoverlap    = 0;


%% Rodent 1
config{1}                           = configcommon;
config{1}.prefix                    = 'DTX5-';
config{1}.rawdir                    = fullfile(rootpath_data, 'DTX5-M1-10uM', '2019_03_19_DTX-5');
config{1}.imagesavedir              = fullfile(imagesavedir, 'DTX5');       
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

config{1}.align.channel.SlowWave    = 'E12LFP';                                                                                                                                                            
config{1}.LFP.channel               = config{1}.labels.macro;
config{1}.LFP.electrodetoplot.SlowWave = 'E12LFP';
config{1}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{1}.circus.channelname        = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};

config{1}.probeoffset = 88;

%% Rodent 2
config{2}                           = configcommon;
config{2}.prefix                    = 'DTX2-';
config{2}.rawdir                    = fullfile(rootpath_data, 'DTX2-M1-10uM', '2019_03_01_DTX-2');
config{2}.imagesavedir              = fullfile(imagesavedir, 'DTX2');       
config{2}.directorylist{1}          =  {'2019-03-01_12-14',...
    '2019-03-01_12-33',...
    '2019-03-01_14-14',...
    '2019-03-01_16-14',...
    '2019-03-01_18-14',...
    '2019-03-01_20-14'};

config{2}.labels.micro              = {'E09','E10','E11','E12'};
config{2}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGS1','ECoGM1','ECoGPtA'}; % /!\ Mistake during acquisition : ECoGS1 = ECoGM1G, ECoGM1 = ECoGM1D. Corrected in the code with dtx_correctDTX2name

config{2}.injectiontime             = datetime('01-Mar-2019 12:33:00');

config{2}.align.channel.SlowWave             = 'E13LFP';
config{2}.LFP.channel               = config{2}.labels.macro;
config{2}.LFP.electrodetoplot.SlowWave       = 'E13LFP'; 

config{2}.circus.channel            = {'E09','E10','E11','E12'};
config{2}.circus.channelname        = {'E09','E10','E11','E12'};

config{2}.probeoffset = 264;


%% Rodent 3
config{3}                           = configcommon;
config{3}.prefix                    = 'DTX4-';
config{3}.rawdir                    = fullfile(rootpath_data, 'DTX4-M1-10uM', '2019_03_08_DTX-4');
config{3}.imagesavedir              = fullfile(imagesavedir,'DTX4');       
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

config{3}.probeoffset = 477;

%% Rodent 4
config{4}                           = configcommon;
config{4}.prefix                    = 'DTX7-';
config{4}.rawdir                    = fullfile(rootpath_data, 'DTX7-M1-10uM', '2019_03_22_DTX-7');
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX7');       
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

config{4}.probeoffset = 176;

%% Rodent 5
config{5}                           = configcommon;
config{5}.prefix                    = 'DTX6-';
config{5}.rawdir                    = fullfile(rootpath_data, 'DTX6-M1-10uM', '2019_03_21_DTX-6');
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX6');       
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

config{5}.probeoffset = 250;

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

config{6}.probeoffset = 315;

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

config{7}.probeoffset = 246;