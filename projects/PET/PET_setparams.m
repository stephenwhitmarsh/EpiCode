function [config] = PET_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/vn_PET';
    rootpath_data       = '/network/lustre/iss02/epimicro/patients/raw';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\vn_PET\';
    rootpath_data       = '\\l2export\iss02.epimicro\patients\raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');

%% patient 1
config{1}                               = [];
config{1}.datasavedir                   = datasavedir;                                    %path where to save data
config{1}.imagesavedir                  = imagesavedir;                                   %path where to save images

config{1}.prefix                        = '3046-';                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir                        = fullfile(rootpath_data,'pat_03046_1482','eeg'); %path to patient data
config{1}.directorylist{1}              = {'03046_2021-07-08_01-51'};                     %list of folders to analyse

config{1}.name                          = {'spike_mHa'};
config{1}.muse.startmarker.spike_mHa    = {'Spike_mHaT2_1'};
config{1}.muse.endmarker.spike_mHa      = {'Spike_mHaT2_1'};

config{1}.LFP.name                      = {'spike_mHa'};
config{1}.LFP.channel                   = {'_mHaT2_1','_mHaT2_2','_mHaT2_3','_mHaT2_4','_mHaT2_5'};
config{1}.LFP.hpfilter                  = 'no';
config{1}.LFP.hpfreq                    = 1;
config{1}.LFP.resamplefs                = 250;
config{1}.LFP.baseline                  = 'yes';
config{1}.LFP.baselinewindow.spike_mHa  = [-0.15, -0.05];

config{1}.circus.channel                        = {'mHaT2_1', 'mHaT2_3', 'mHaT2_4','mHaT2_6', 'mHaT2_7', 'mHaT2_8'};
config{1}.circus.reref                          = 'no';
config{1}.circus.refchan                        = '';
config{1}.circus.outputdir                      = 'SpykingCircus';
config{1}.circus.maxchan                        = 'phy';

config{1}.circus.paramfile                      = fullfile(rootpath_analysis, 'EpiCode', 'projects', 'PET', 'SpykingCircus.params');
config{1}.circus.params.detection.spike_thresh  = '6';
config{1}.circus.params.filtering.cut_off       = '300, auto';
config{1}.circus.params.filtering.remove_median = 'False';
config{1}.circus.params.clustering.max_elts     = '20000';
config{1}.circus.params.data.stream_mode        = 'mapping-file';
config{1}.circus.params.data.mapping_file       = 'filelist.txt';
config{1}.circus.params.detection.peaks         = 'negative';

%% patient 2
config{2}                               = config{1};

config{2}.prefix                        = '3046-';                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir                        = fullfile(rootpath_data,'pat_03046_1482','eeg'); %path to patient data
config{2}.directorylist{1}              = {'03046_2021-07-08_01-51'};                     %list of folders to analyse

config{2}.name                          = {'spike_mHa'};
config{2}.muse.startmarker.spike_mHa    = {'Spike_mHaT2_1'};
config{2}.muse.endmarker.spike_mHa      = {'Spike_mHaT2_1'};

config{2}.LFP.name                      = {'spike_mHa'};
config{2}.LFP.channel                   = {'_mHaT2_1','_mHaT2_2','_mHaT2_3','_mHaT2_4','_mHaT2_5'};
config{2}.LFP.hpfilter                  = 'no';
config{2}.LFP.hpfreq                    = 1;
config{2}.LFP.resamplefs                = 250;
config{2}.LFP.baseline                  = 'yes';
config{2}.LFP.baselinewindow.spike_mHa  = [-0.15, -0.05];

config{2}.circus.channel                        = {'mHaT2_1', 'mHaT2_3', 'mHaT2_4','mHaT2_6', 'mHaT2_7', 'mHaT2_8'};
config{2}.circus.reref                          = 'no';
config{2}.circus.refchan                        = '';
config{2}.circus.outputdir                      = 'SpykingCircus';
config{2}.circus.maxchan                        = 'phy';

