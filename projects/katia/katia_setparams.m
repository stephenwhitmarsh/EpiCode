function [config] = katia_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss02/charpier/analyses/vn_pet';
    rootpath_data       = '/network/lustre/iss02/epimicro/patients/raw';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\l2export\iss02.charpier\analyses\vn_pet\';
    rootpath_data       = '\\l2export\iss02.epimicro\patients\raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');
tablesavedir =  fullfile(rootpath_analysis,'tables');

%% configuration common for all patients

configcommon.name = {'IED', 'interIED'};
configcommon.datasavedir            = datasavedir;                                    %path where to save data
configcommon.imagesavedir           = imagesavedir;                                   %path where to save images
configcommon.tablesavedir           = tablesavedir;                                   %path where to save images

% for MUSE
configcommon.name                          = {'Hspike'}; % name for general bookkeeping
configcommon.muse.startmarker.Hspike       = "Hspike"; % label name in muse
configcommon.muse.endmarker.Hspike         = "Hspike"; % label name in muse
configcommon.muse.startmarker.template1    = "template1";
configcommon.muse.endmarker.template1      = "template1";

%define epochs IED
configcommon.epoch.toi.Hspike       = [-0.5  1]; %entre -0.5s avant le marker et +1 apres le marker
configcommon.epoch.pad.Hspike       = 1;
configcommon.spike.toi.IED          = [-0.5  1];
configcommon.spike.pad.IED          = 0;
configcommon.min_trial_length.interIED = 0; %seconds

%inter IED
configcommon.epoch.toi.interIED     = [2 -2]; 
configcommon.epoch.pad.interIED     = 0;
configcommon.spike.toi.interIED     = [2 -2];%1s apres la premiere pointe et 1s avant la pointe suivante
configcommon.spike.pad.interIED     = 0;
configcommon.min_trial_length.interIED = 5; %seconds

%LFP parameters
configcommon.LFP.name                      = {'IED', 'interIED'};
configcommon.LFP.resamplefs                = 250;
configcommon.LFP.baseline                  = 'yes';
configcommon.LFP.baselinewindow.IED        = [-0.15, -0.05];
configcommon.LFP.baselinewindow.interIED   = [1 2];
configcommon.LFP.lpfilter                  = 'yes';
configcommon.LFP.lpfreq                    = 100; %Hz
configcommon.LFP.bsfilter                  = 'yes';
configcommon.LFP.bsfreq                    = [49 51];
configcommon.LFP.bsinstabilityfix          = 'reduce';

% align xcorr
configcommon.align.name                = {'IED'};
configcommon.align.latency.IED         = [-0.5 1];
configcommon.align.removeartefacts     = true;
configcommon.align.reref               = 'no';
configcommon.align.lpfilter            = 'yes';
configcommon.align.lpfreq              = 30;
configcommon.align.reject              = 'BAD_cnt';

configcommon.TFR.name                    = {'IED'};
configcommon.TFR.keeptrials              = 'yes';
configcommon.TFR.foi.IED                 = 5:1:100;
configcommon.TFR.t_ftimwin.IED           = 7./configcommon.TFR.foi.IED;
configcommon.TFR.toi.IED                 = -1.5 : 0.01 : 2;       

configcommon.FFT.name                    = {'IED', 'interIED'};
configcommon.FFT.foi.IED                 = 5:100;
configcommon.FFT.foi.interIED            = 1:15;


configcommon.circus.channel                        = {'mHaT2_1', 'mHaT2_3', 'mHaT2_4','mHaT2_6', 'mHaT2_7', 'mHaT2_8'};
configcommon.circus.reref                          = 'no';
configcommon.circus.refchan                        = '';
configcommon.circus.outputdir                      = 'SpykingCircus';
configcommon.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'hspike', 'SpykingCircus.params');
configcommon.circus.params.detection.spike_thresh  = '6';
configcommon.circus.params.filtering.cut_off       = '300, auto';
configcommon.circus.params.filtering.remove_median = 'False';
configcommon.circus.params.clustering.max_elts     = '20000';
configcommon.circus.params.data.stream_mode        = 'mapping-file';
configcommon.circus.params.data.mapping_file       = 'filelist.txt';
configcommon.circus.params.detection.peaks         = 'negative';

%stats parameters
configcommon.spike.name             = {'IED', 'interIED'};
configcommon.spike.ISIbins          = 0:0.001:0.025;
configcommon.spike.psthbin.IED      = 1/100; %s
configcommon.spike.psthbin.interIED = 1; %s

configcommon.stats.bl.IED      = [-0.5 -0.15];
configcommon.stats.bl.interIED = [1 10];
configcommon.stats.alpha       = 0.05;
configcommon.spike.RPV         = 0.002;

%% patient 1
config{1}                               = configcommon;

config{1}.prefix                        = '3046-';                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir                        = fullfile(rootpath_data,'pat_03046_1482','eeg'); %path to patient data
config{1}.directorylist{1}              = {'03046_2021-07-07_15-51'};                     %list of folders to analyse
% config{1}.marker_ied_names = {..., ..., ...}; %name of all different muse
% markers used to anotate IDED on the different microelectrodes

config{1}.verifymarkers.one_off = 'mHaT2_4'; %{'Spike_mHaT2_1','Spike_mHaT2_2'}

%define IED epochs
config{1}.muse.startmarker.IED    = 'mHaT2_4';
config{1}.muse.endmarker.IED      = 'mHaT2_4';

%define interIED epochs
config{1}.muse.startmarker.interIED     = {'mHaT2_4', 'Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}
config{1}.muse.endmarker.interIED       = {'mHaT2_4', 'Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}

config{1}.LFP.channel                   = {'_mHaT2_1','_mHaT2_2', 'mHaT2_8', 'mHmT2_5', 'mHmT2_8'}; %include all electrodes used for spike sorting and all electrodes with IED used for alignment
config{1}.align.channel                 = {'_mHaT2_1','_mHaT2_2'}; % all electrodes with ied
config{1}.alignpeak.channel.IED         = '_mHaT2_1'; % one electrode, with the better ied morphology, used for alignment

config{1}.circus.channel                        = {'mHaT2_1', 'mHaT2_2', 'mHaT2_8', 'mHmT2_5', 'mHmT2_8'};
config{1}.circus.channelname                    = {'mHaT2',   'mHaT2',   'mHaT2',   'mHmT2',   'mHmT2'};

config{1}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'PET', 'SpykingCircus.params');
config{1}.circus.params.detection.spike_thresh  = '7';
config{1}.circus.params.filtering.cut_off       = '300, auto';
config{1}.circus.params.filtering.remove_median = 'False';
config{1}.circus.params.clustering.max_elts     = '20000';
config{1}.circus.params.clustering.nb_repeats   = '6';
config{1}.circus.params.fitting.collect_all     = 'True'; %create one garbage template with unfitted events which crossed the threshold
config{1}.circus.params.data.stream_mode        = 'mapping-file';
config{1}.circus.params.data.mapping_file       = 'filelist.txt';
config{1}.circus.params.detection.peaks         = 'negative';
config{1}.circus.outputdir                      = 'SpykingCircus';
config{1}.circus.common_ground                  = []; %pour rereferencer, mettre le numero de l'electrode moins 1 (car commence a compter a zero)


%% patient 2
config{2}                               = configcommon;

config{2}.prefix                        = '3138-';                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir                        = fullfile(rootpath_data,'pat_03138_1601','eeg'); %path to patient data
config{2}.directorylist{1}              = {'03138_2022-03-17_10-23'};                     %list of folders to analyse
% config{1}.marker_ied_names = {..., ..., ...}; %name of all different muse
% markers used to anotate IDED on the different microelectrodes

config{2}.verifymarkers.one_off = {'spike_mAmT2','Hspike'}; %{'Spike_mHaT2_1','Spike_mHaT2_2'}

%define IED epochs
config{2}.muse.startmarker.IED    = 'spike_mAmT2';
config{2}.muse.endmarker.IED      = 'spike_mAmT2';

%define interIED epochs
config{2}.muse.startmarker.interIED     = {'spike_mAmT2','Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}
config{2}.muse.endmarker.interIED       = {'spike_mAmT2','Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}

config{2}.LFP.channel                   = {'_mAmT2_1','_mAmT2_2','_mAmT2_5', '_mAmT2_6','_mAmT2_7' ,'_mAmT2_8'};  %include all electrodes used for spike sorting and all electrodes with IED used for alignment
config{2}.align.channel                 = {'_mAmT2_5','_mAmT2_7','_mAmT2_8'};  % all electrodes with ied
config{2}.alignpeak.channel.IED         = '_mAmT2_5'; % one electrode, with the better ied morphology, used for alignment

config{2}.circus.channel                = {'_mAmT2_1','_mAmT2_2','_mAmT2_5','_mAmT2_6','_mAmT2_7'};
% config{2}.circus.channel                = {'_mAmT2_5','_mAmT2_7'};
config{2}.circus.channelname            = {'mAmT2','mAmT2','mAmT2','mAmT2','mAmT2'};

config{2}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'PET', 'SpykingCircus.params');
config{2}.circus.params.detection.spike_thresh  = '6';
config{2}.circus.params.filtering.cut_off       = '300, auto';
config{2}.circus.params.filtering.remove_median = 'False';
config{2}.circus.params.clustering.max_elts     = '20000';
config{2}.circus.params.clustering.nb_repeats   = '6';
config{2}.circus.params.fitting.collect_all     = 'True'; %create one garbage template with unfitted events which crossed the threshold
config{2}.circus.params.data.stream_mode        = 'mapping-file';
config{2}.circus.params.data.mapping_file       = 'filelist.txt';
config{2}.circus.params.detection.peaks         = 'negative';
config{2}.circus.outputdir                      = 'SpykingCircus';


%% patient 3
config{3}                               = configcommon;

config{3}.prefix                        = '2718-';                                        %patient name. Must end by "-". namepatient-
config{3}.rawdir                        = fullfile(rootpath_data,'pat_02718_1201','eeg'); %path to patient data
config{3}.directorylist{1}              = {'02718_2019-05-16_15-15'};                     %list of folders to analyse
% config{1}.marker_ied_names = {..., ..., ...}; %name of all different muse
% markers used to anotate IDED on the different microelectrodes

config{3}.verifymarkers.one_off = {'Spike_mHaT1_3','Hspike2'}; %{'Spike_mHaT2_1','Spike_mHaT2_2'}

%define IED epochs
config{3}.muse.startmarker.IED    = 'Spike_mHaT1_3';
config{3}.muse.endmarker.IED      = 'Spike_mHaT1_3';

%define interIED epochs
config{3}.muse.startmarker.interIED     = {'Spike_mHaT1_3','Hspike2'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}
config{3}.muse.endmarker.interIED       = {'Spike_mHaT1_3','Hspike2'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}

config{3}.LFP.channel                   = {'_mHaT1_3','_mHaT1_5','_mHaT1_6','_mHaT1_7'};
config{3}.align.channel                 = {'_mHaT1_3','_mHaT1_5','_mHaT1_6','_mHaT1_7'};
config{3}.alignpeak.channel.IED         = '_mHaT1_3';

config{3}.circus.channel                = {'_mHaT1_7'};
config{3}.circus.channelname            = {'mHaT1'};

config{3}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'PET', 'SpykingCircus.params');
config{3}.circus.params.detection.spike_thresh  = '6';
config{3}.circus.params.filtering.cut_off       = '300, auto';
config{3}.circus.params.filtering.remove_median = 'False';
config{3}.circus.params.clustering.max_elts     = '20000';
config{3}.circus.params.clustering.nb_repeats   = '6';
config{3}.circus.params.fitting.collect_all     = 'True'; %create one garbage template with unfitted events which crossed the threshold
config{3}.circus.params.detection.peaks         = 'negative';
config{3}.circus.params.data.stream_mode        = 'mapping-file';
config{3}.circus.params.data.mapping_file       = 'filelist.txt';
config{3}.circus.outputdir                      = 'SpykingCircus';


%% patient 4
config{4}                               = configcommon;

config{4}.prefix                        = '2711-';                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir                        = fullfile(rootpath_data,'pat_02711_1193','eeg'); %path to patient data
config{4}.directorylist{1}              = {'02711_2019-04-18_15-04'};                     %list of folders to analyse
% config{1}.marker_ied_names = {..., ..., ...}; %name of all different muse
% markers used to anotate IDED on the different microelectrodes

config{4}.verifymarkers.one_off = {'mHaT2','Hspike'}; %{'Spike_mHaT2_1','Spike_mHaT2_2'}

%define IED epochs
config{4}.muse.startmarker.IED    = 'mHaT2';
config{4}.muse.endmarker.IED      = 'mHaT2';

%define interIED epochs
config{4}.muse.startmarker.interIED     = {'mHaT2','Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}
config{4}.muse.endmarker.interIED       = {'mHaT2','Hspike'}; %{'Spike_mHaT2_1', 'name_macro1', 'name_macro2'}

config{4}.LFP.channel                   = {'_mHaT2_1','_mHaT2_3','_mHaT2_4','_mHaT2_5','_mHaT2_6','_mHaT2_7','_mHaT2_8'};
config{4}.align.channel                 = {'_mHaT2_7','_mHaT2_8'};
config{4}.alignpeak.channel.IED         = '_mHaT2_8';

config{4}.circus.channel                = {'_mHaT2_1','_mHaT2_3','_mHaT2_4','_mHaT2_6','_mHaT2_7','_mHaT2_8'};
config{4}.circus.channelname            = {'mHaT2','mHaT2','mHaT2','mHaT2','mHaT2','mHaT2'};

config{4}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'PET', 'SpykingCircus.params');
config{4}.circus.params.detection.spike_thresh  = '6';
config{4}.circus.params.filtering.cut_off       = '300, auto';
config{4}.circus.params.filtering.remove_median = 'False';
config{4}.circus.params.clustering.max_elts     = '20000';
config{4}.circus.params.clustering.nb_repeats   = '6';
config{4}.circus.params.fitting.collect_all     = 'True'; %create one garbage template with unfitted events which crossed the threshold
config{4}.circus.params.detection.peaks         = 'negative';
config{4}.circus.params.data.stream_mode        = 'mapping-file';
config{4}.circus.params.data.mapping_file       = 'filelist.txt';
config{4}.circus.outputdir                      = 'SpykingCircus';
