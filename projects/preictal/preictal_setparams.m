function [config] = preictal_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/vn_preictal';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\vn_preictal\';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');

%% config common for all patients
configcommon.name                  = {'Preictal'};
configcommon.circus.deadfilesuffix = '_SeizuresNotRemoved';
configcommon.muse.write            = true;
configcommon.statstime.timewin     = 10;
configcommon.statstime.slidestep   = 1;
configcommon.statstime.minbadtime  = 1;
configcommon.statstime.write       = true;
configcommon.spike.ISIbins         = 0:0.001:0.05;


%% patient 1
config{1}                     = configcommon;
config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.prefix              = 'pat_02256_0700_Crise1_m2mCi-';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'pat_02256_0700','eeg');                       %path to patient data
config{1}.directorylist{1}    = {'02256_2015-05-09_05-37','02256_2015-05-09_07-37'};                                               %list of folders to analyse
config{1}.circus.channel      = {'m2mCi_2','m2mCi_3','m2mCi_5','m2mCi_6','m2mCi_7','m2mCi_8'};       %name of the first electrode
config{1}.circus.reref        = 'yes';
config{1}.circus.refchan      = 'm2mCi_4';
config{1}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{1}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{1}.circus.hpfreq       = 0; % even when not using
config{1}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{1}.seizure_index       = 'last'; %Optional. index of the seizure to analyze, on the LAST dir. can be 'last' (default)
% %% patient 2
% config{2} = config{1};
% config{2}.prefix              = 'pat_02256_0700_Crise1_m2mCi-';   
% config{2}.circus.channel      = {'m2mCi_2','m2mCi_3','m2mCi_4','m2mCi_5','m2mCi_6','m2mCi_7','m2mCi_8'};       %name of the first electrode
% config{2}.circus.reref        = 'yes';
% config{2}.circus.refchan      = 'm2mCi_4';


%% patient 2
config{2}                     = configcommon;
config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.prefix              = 'pat_02256_0700_Crise2_m2mCi-';                                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'pat_02256_0700','eeg');                       %path to patient data
config{2}.directorylist{1}    = {'02256_2015-05-21_03-28','02256_2015-05-21_05-28'};                                               %list of folders to analyse
config{2}.circus.channel      = {'m2mCi_2','m2mCi_5','m2mCi_6','m2mCi_7','m2mCi_8'};       %name of the first electrode
config{2}.circus.reref        = 'yes';
config{2}.circus.refchan      = 'm2mCi_4';
config{2}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{2}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{2}.circus.hpfreq       = 0; % even when not using
config{2}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{2}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{2}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{2}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{2}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{2}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{2}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 

%% patient 3
config{3}                     = configcommon;
config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.prefix              = 'pat_02379_0828_Crise1_mHaBg-';                                                        %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'pat_02379_0828','eeg');                       %path to patient data
config{3}.directorylist{1}    = {'02379_2016-05-30_13-45','02379_2016-05-30_15-45'};                                               %list of folders to analyse
config{3}.circus.channel      = {'mHaBg_6'};       %name of the first electrode
config{3}.circus.reref        = 'no';
config{3}.circus.refchan      = '';
config{3}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{3}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{3}.circus.hpfreq       = 0; % even when not using
config{3}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{3}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{3}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{3}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{3}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{3}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{3}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 4
config{4}                     = configcommon;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.prefix              = 'pat_02379_0828_Crise1_mHa2d-';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'pat_02379_0828','eeg');                       %path to patient data
config{4}.directorylist{1}    = {'02379_2016-05-30_13-45','02379_2016-05-30_15-45'};                                               %list of folders to analyse
config{4}.circus.channel      = {'mHa2d_2'};       %name of the first electrode
config{4}.circus.reref        = 'no';
config{4}.circus.refchan      = '';
config{4}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{4}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{4}.circus.hpfreq       = 0; % even when not using
config{4}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{4}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{4}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{4}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{4}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{4}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{4}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 5
config{5}                     = configcommon;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.prefix              = 'pat_02379_0828_Crise2_mHaBg-';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'pat_02379_0828','eeg');                       %path to patient data
config{5}.directorylist{1}    = {'02379_2016-06-05_15-44','02379_2016-06-05_17-44'};                                               %list of folders to analyse
config{5}.circus.channel      = {'mHaBg_4'};       %name of the first electrode
config{5}.circus.reref        = 'no';
config{5}.circus.refchan      = '';
config{5}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{5}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{5}.circus.hpfreq       = 0; % even when not using
config{5}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{5}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{5}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{5}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{5}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{5}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{5}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 6
config{6}                     = configcommon;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.prefix              = 'pat_02379_0828_Crise2_mHa2d-';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'pat_02379_0828','eeg');                       %path to patient data
config{6}.directorylist{1}    = {'02379_2016-06-05_15-44','02379_2016-06-05_17-44'};                                               %list of folders to analyse
config{6}.circus.channel      = {'mHa2d_2'};       %name of the first electrode
config{6}.circus.reref        = 'no';
config{6}.circus.refchan      = '';
config{6}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{6}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{6}.circus.hpfreq       = 0; % even when not using
config{6}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{6}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{6}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{6}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{6}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{6}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{6}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 7
config{7}                     = configcommon;
config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.prefix              = 'pat_02599_1057_Crise1_mHaT2-';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'pat_02599_1057','eeg');                       %path to patient data
config{7}.directorylist{1}    = {'02599_2018-04-25_11-24','02599_2018-04-25_13-24'};                                               %list of folders to analyse
config{7}.circus.channel      = {'mHaT2_6'};       %name of the first electrode
config{7}.circus.reref        = 'no';
config{7}.circus.refchan      = '';
config{7}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{7}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{7}.circus.hpfreq       = 0; % even when not using
config{7}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{7}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{7}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{7}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{7}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{7}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{7}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 8
config{8}                     = configcommon;
config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.prefix              = 'pat_02599_1057_Crise2_mHaT2-';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'pat_02599_1057','eeg');                       %path to patient data
config{8}.directorylist{1}    = {'02599_2018-04-25_13-24','02599_2018-04-25_15-24'};                                               %list of folders to analyse
config{8}.circus.channel      = {'mHaT2_6'};       %name of the first electrode
config{8}.circus.reref        = 'no';
config{8}.circus.refchan      = '';
config{8}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{8}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{8}.circus.hpfreq       = 0; % even when not using
config{8}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{8}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{8}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{8}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{8}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{8}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{8}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 
%% patient 9
config{9}                     = configcommon;
config{9}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{9}.imagesavedir        = imagesavedir;
config{9}.prefix              = 'pat_02599_1057_Crise3_mHaT2-';                                                        %patient name. Must end by "-". namepatient-
config{9}.rawdir              = fullfile(rootpath_data,'pat_02599_1057','eeg');                       %path to patient data
config{9}.directorylist{1}    = {'02599_2018-04-28_11-23','02599_2018-04-28_13-23'};                                               %list of folders to analyse
config{9}.circus.channel      = {'mHaT2_6'};       %name of the first electrode
config{9}.circus.reref        = 'no';
config{9}.circus.refchan      = '';
config{9}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{9}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{9}.circus.hpfreq       = 0; % even when not using
config{9}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
config{9}.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
config{9}.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
config{9}.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
config{9}.bad.sample_list     = 'last'; %dernier marqueur crise END pris en compte (en cas de multiple crise sur un même fichier)
config{9}.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
config{9}.bad.removeseizures  = 'no'; %pour enlever les crises mettre 'yes' 