
%% Setting parameters DTX project Paul Baudin

function [config] = dtx_setparams_probe_spikes(config)

disp('setting parameters for probe spike analysis');
muscale = 50; % multiunit scale

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

datasavedir = fullfile(rootpath_analysis, 'data', 'spike'); %removed of config{i} so we can more easily modify it
imagesavedir = fullfile(rootpath_analysis, 'image_spike');

%% Config for dtx experiments

configdtx.type                      = 'dtx';
configdtx.os                        = os;
configdtx.datasavedir               = datasavedir;
configdtx.name                      = {'SlowWave','Seizure','Interictal','SlowWave_Larger','xcorr_10_1', 'xorr_2_1', 'xcorr_1_0', 'xcorr_1_1'};
configdtx.muse.startend             = {'SlowWave','SlowWave'; 'SlowWave', 'Crise_End';'Crise_End','SlowWave';'SlowWave','SlowWave';'SlowWave','SlowWave';'SlowWave','SlowWave';'SlowWave','SlowWave';'SlowWave','SlowWave'};   % 'SlowWave','SlowWave'; for readLFP function : cut data ...s before SlowWave, and ...s after SlowWave

% list of onset timing with respect to start-marker (s)
configdtx.epoch.toi{1}              = [-2, 2];
configdtx.epoch.toi{2}              = [-2, 2];
configdtx.epoch.toi{3}              = [2, -1];
configdtx.epoch.toi{4}              = [-120, 30];
configdtx.epoch.toi{5}              = [-10, -1];
configdtx.epoch.toi{6}              = [-2, -1];
configdtx.epoch.toi{7}              = [-1, 0];
configdtx.epoch.toi{8}              = [-1, 1];
configdtx.epoch.pad{1}              = 10; %for LFP
configdtx.epoch.pad{2}              = 10;
configdtx.epoch.pad{3}              = 10;
configdtx.epoch.pad{4}              = 10;
configdtx.epoch.pad{5}              = 10;
configdtx.epoch.pad{6}              = 10;
configdtx.epoch.pad{7}              = 10;
configdtx.epoch.pad{8}              = 10;

configdtx.commonchans               = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

configdtx.bad.markerStart              = 'SlowWave';
configdtx.bad.markerEnd                = 'Crise_End';
configdtx.bad.time_from_begin          = -2;
configdtx.bad.time_from_end            = 2;

% OLD ALIGN
configdtx.align.name                = {'SlowWave'};
configdtx.align.flip                = {'no'};
configdtx.align.abs                 = {'no'};
configdtx.align.method              = {'max'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configdtx.align.filter              = {'lp'};
configdtx.align.freq                = {10};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
configdtx.align.hilbert             = {'no'};
configdtx.align.thresh.value        = [1, 1];
configdtx.align.thresh.method       = {'trial'};%,'trial','trial'};%'medianbl','both';
configdtx.align.toiplot             = {[-1,  1], [-1, 1]};                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configdtx.align.toiactive           = {[-0.5, 0.5], [-0.5, 0.5]};                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configdtx.align.toibaseline         = {[-1, -0.5], [-1, -0.5]};
configdtx.align.maxtimeshift        = {0.3};
configdtx.align.demean              = {'yes'};
%for detection of begin of event
configdtx.align.begin.doalign       = {'yes'};
configdtx.align.begin.thresh        = {0.3}; % percent of peak
% % ALIGN XCORR : 
% ATTENTION RETIRER ALIGN.CHANNEL RAT PAR RAT
% configdtx.align.name                = {'SlowWave'};                               % Name of markers/patterns to align
% configdtx.align.channel             = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP'};    % Channels to use for alignment
% configdtx.align.demean              = 'no';
% configdtx.align.reref               = 'no';
% configdtx.align.latency             = {[-1, 1]};  

configdtx.LFP.flip                  = 'false';
configdtx.LFP.name                  = {'SlowWave'};%_EEG', 'SLowWave_Intra'};
configdtx.LFP.hpfilter              = 'no';
configdtx.LFP.resamplefs            = 320; %because sampling rate is 3200Hz
configdtx.LFP.baseline              = 'no';
configdtx.LFP.baselinewindow{1}     = [-2, -1];
configdtx.LFP.baselinewindow{2}     = [-2, -1];
configdtx.LFP.baselinewindow{3}     = [0, 1];
configdtx.LFP.write                 = true;

configdtx.bad.markerStart           = 'Crise_Start';
configdtx.bad.markerEnd             = 'Crise_End';
configdtx.bad.time_from_begin       = -2;
configdtx.bad.time_from_end         = 2;    

configdtx.circus.reref              = 'no';
configdtx.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'dtx', 'SpykingCircus');
configdtx.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
configdtx.circus.postfix            = [];%'-final';
configdtx.circus.deadfilesuffix     = 'withseizures';

configdtx.stats.numrandomization     = 1000; %Essayer à 1000
% configdtx.stats.dostat              = {true, true, false}; %not used
configdtx.stats.bltoi{1}            = [-2, -1];
configdtx.stats.bltoi{2}            = [-2, -1];
configdtx.stats.bltoi{3}            = [2, 6];
configdtx.stats.alpha               = 0.025;
%not used for now (but maybe to better select ISI or waveforms of interest) :
configdtx.stats.actoi{1}            = [-1, 1];
configdtx.stats.actoi{2}            = [-1, Inf];
configdtx.stats.actoi{3}            = [6, Inf];

% to smooth spikerate. Not used for now. Maybe for correlation LFP-spike
% configdtx.spike.toispikerate{1}     = [-0.01 0.01];            % for plotting spikerate
% configdtx.spike.toispikerate{1}     = [-1 1];                  % for plotting spikerate
% configdtx.spike.toispikerate{2}     = [-0.05 0.05];            % for plotting spikerate
% configdtx.spike.toispikerate{3}     = [-10 10];                % for plotting spikerate
configdtx.spike.eventsname          = {'SlowWave','Seizure'};%, 'Seizure'};
configdtx.spike.baselinename        = 'Interictal';
configdtx.spike.resamplefs{1}       = 50; %1/size of bar graph bins
configdtx.spike.resamplefs{2}       = 50;
configdtx.spike.resamplefs{3}       = 0.1;
configdtx.spike.RPV                 = 0.003; %refractory period violation, in seconds
configdtx.spike.ISIbins             = [0:0.003:0.150]; %in s
configdtx.spike.sdflatency{1}   = 'maxperiod';
configdtx.spike.sdflatency{2}   = 'maxperiod';
configdtx.spike.sdflatency{3}   = 'maxperiod';
configdtx.spike.sdftimewin{1}   = [-1/2 * 1/50 1/2 * 1/50]; %1/resamplef, so no effect of the convolution
configdtx.spike.sdftimewin{2}   = [-1/2 * 1/50 1/2 * 1/50];
configdtx.spike.sdftimewin{3}   = [-1/2 * 1/0.1 1/2 * 1/0.1];

configdtx.spikewaveform.toi         = [-0.0015 0.0015]; %in s
configdtx.spikewaveform.cutoff      = 300; %high pass filter frequency to apply to raw data
configdtx.spikewaveform.nspikes     = 200000;%'all'; %maximum number of spike waveforms to load. Can be 'all'. 

%% config for control experiments
% no seizure-time-locked trials
% markers : Baseline_Start, Analysis_End, Injection, BAD
configctrl.type                      = 'ctrl';
configctrl.os                        = os;
configctrl.datasavedir               = datasavedir;
configctrl.name                      = {'Control'};

configctrl.commonchans               = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

% configctrl.muse.startend             = {'Injection','Analysis_End'};  

configctrl.LFP.name                  = []; %do not load LFP
configctrl.LFP.channel                  = []; %do not load LFP
configctrl.align.name                = []; %do not align muse markers
configctrl.align.channel                = []; %do not align muse markers

% list of onset timing with respect to start-marker (s)
configctrl.epoch.toi{1}              = [0, 0];
configctrl.epoch.pad{1}              = 10; %for LFP

configctrl.circus.reref              = 'no';
configctrl.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'dtx', 'SpykingCircus');
configctrl.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
configctrl.circus.postfix             = [];%'-final';

configctrl.spike.triallength         = 600; %seconds
configctrl.spike.baselinename        = 'Control';
configctrl.spike.eventsname          = [];
configctrl.spike.RPV                 = 0.001; %refractory period violation, in seconds
configctrl.spike.ISIbins             = [0:0.003:0.150]; %in s

configctrl.spikewaveform.toi         = [-0.0015 0.0015]; %in s
configctrl.spikewaveform.cutoff      = 300; %high pass filter frequency to apply to raw data
configctrl.spikewaveform.nspikes     = 'all'; %maximum number of spike waveforms to load per unit. Can be 'all'. 


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

config{1}.align.channel             = {'E12LFP'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.LFP.channel               = config{1}.labels.macro;
config{1}.LFP.electrodetoplot       = {'ECoGM1G', 'E12LFP'};
config{1}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


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

config{2}.align.channel             = {'E13LFP'};%{'ECoGS1'};% %ATTENTION erreur nom eeg pendant acquisition
config{2}.LFP.channel               = config{2}.labels.macro;%ECoG S1 et ECoGM1 renommés respectivement ECoGM1G et ECoGM1D, après l'étape de readLFP
config{2}.LFP.electrodetoplot       = {'ECoGM1G', 'E13LFP'}; %ATTENTION ECoGS1 est renommé ECoGM1G pendant le script LFP

config{2}.circus.channel            = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};


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

config{3}.align.channel             = {'E13LFP'};%{'ECoGM1G'};%
config{3}.LFP.channel               = config{3}.labels.macro;
config{3}.LFP.electrodetoplot       = {'ECoGM1G', 'E13LFP'};
config{3}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 4
config{4}                           = configdtx;
config{4}.prefix                    = 'DTX7-';
config{4}.rawdir                    = fullfile(rootpath_data, 'DTX7-M1-10uM', '2019_03_22_DTX-7');
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX7');       % where to print images
config{4}.directorylist{1}          =  {'2019-03-22_12-31',...
    '2019-03-22_14-31',...
    '2019-03-22_16-31',...
    '2019-03-22_18-31'};

config{4}.labels.micro              = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{4}.labels.macro              = {'E06LFP','E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{4}.injectiontime             = datetime('22-Mar-2019 13:40:07');

config{4}.align.channel             = {'E13LFP'};
config{4}.LFP.channel               = config{4}.labels.macro;
config{4}.LFP.electrodetoplot       = {'ECoGM1G', 'E13LFP'};
config{4}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 5
%Remove the first 2 files and keep only the 2 lasts because the probe was moved in
%order to have better unit activity
config{5}                           = configdtx;
config{5}.prefix                    = 'DTX6-';
config{5}.rawdir                    = fullfile(rootpath_data, 'DTX6-M1-10uM', '2019_03_21_DTX-6');
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX6');       % where to print images
config{5}.directorylist{1}          =  {'2019-03-21_18-12', '2019-03-21_20-12'};

config{5}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{5}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{5}.injectiontime             = datetime('21-Mar-2019 15:08:00');

config{5}.align.channel             = {'E13LFP'};
config{5}.LFP.channel               = config{5}.labels.macro;
config{5}.LFP.electrodetoplot       = {'ECoGM1G', 'E13LFP'};
config{5}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


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

%% Rodent 7
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

config{7}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


end




