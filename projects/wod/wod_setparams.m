function [config] = wod_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')% peut etre changer avec les bonnes paths
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Antoine';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/16_ch/Extra_Neuralynx';
    rootpath_concatdata = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/16_ch/concatenated_LFP';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Antoine';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\16_ch\Extra_Neuralynx';
    rootpath_concatdata = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\16_ch\concatenated_LFP';
    os                  = 'windows';
else
    error('Platform not supported')
end 

concatdata_path= fullfile(rootpath_concatdata);
datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');
imagesavedir_data= {fullfile(imagesavedir,'TFR'), fullfile(imagesavedir,'band_depth'),fullfile(imagesavedir,'TFR','Log'),fullfile(imagesavedir,'LFHF_ratio')};
statsavedir=fullfile(rootpath_analysis,'stats');
freqstat_path= fullfile(statsavedir,'freq_data');
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all rats
configcommon.muse.templatemarker   = fullfile(datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

configcommon.name                  = {'WoD'};
configcommon.LFP.allchannel        = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'Puff'};
configcommon.LFP.name              = configcommon.name;

configcommon.muse.backupdir            = fullfile(datasavedir,'Backup_MuseMarker');
configcommon.muse.startmarker.WoD      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD             = [-900, 3400];
configcommon.epoch.pad.WoD             = 5;

configcommon.LFP.resamplefs                     = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                          = true; %save computed data to disk

configcommon.LFP.lpfilter_wod_detection         = 4;%Hz
configcommon.LFP.wod_toisearch                  = [-10 50]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [-1 50]; %s, were to find wor positive peak relative to the muse marker REAL VALUE -1 25
configcommon.LFP.hpfilter_wod_exclusion         = 1; %Hz

configcommon.timefreq.foi          = [1:2:100];%[1:2:100] is right value
configcommon.timefreq.foi_band     = {[1 9],[10 49],[55 80],[80 100]};%Hz
configcommon.timefreq.t_ftimwin.long    = 10;% in second, length of the time window 10 (right value)
configcommon.timefreq.t_ftimwin.short = 2; % in seconds
configcommon.timefreq.timestep.long     = 2.5;% in second, time between 2 sliding time windows. can be 'all' 2.5 right value
configcommon.timefreq.timestep.short     = 0.5;% in second, time between 2 sliding time windows. can be 'all'
configcommon.timefreq.movmeanwin   = [1,1,1,100,100,100];%in sample points, one value per analysis_name
configcommon.timefreq.tapsmofrq.long    = 2; %2
configcommon.timefreq.tapsmofrq.short    = 2; %2
configcommon.timefreq.toi.long           = [-900 3400];
configcommon.timefreq.toi.short           = [-900 600];
configcommon.timefreq.HF           = [80 100];
configcommon.timefreq.MF           = [55 80];
configcommon.timefreq.MLF          = [10 49];
configcommon.timefreq.LF           = [1 9];

configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Antoine.params');
configcommon.circus.writedeadfile  = 'no';
configcommon.circus.reref        = 'no';
configcommon.circus.refchan      = [];
configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.hpfreq       = 0; % even when not using
configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number


%% rat 4
config{4}                     = configcommon;
config{4}.concatdata_path     = concatdata_path;
config{4}.statsavedir         = statsavedir;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.imagesavedir_data   = imagesavedir_data;
config{4}.prefix              = 'Rat-19_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'2020_05_19_WOD');                       %path to patient data

config{4}.directorylist{1}    = {'2020-05-19_14-25', '2020-05-19_15-44', '2020-05-19_16-25'}; %liste de tous les fichiers, tous les protocoles
config{4}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP','Puff'};
config{4}.LFP.rename          = {'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{4}.LFP.chan_depth         = {565, 815, 1065, 1315, 1565, 1815};
config{4}.LFP.origin_WoD          = {'E14', 'E14'};
config{4}.LFP.origin_WoR          = {'E11', 'E11'};



config{4}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{4}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 5
config{5}                     = configcommon;
config{5}.concatdata_path     = concatdata_path;
config{5}.statsavedir         = statsavedir;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.imagesavedir_data   = imagesavedir_data;
config{5}.prefix              = 'Rat-25_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'2020_05_25_WOD');                       %path to patient data

config{5}.directorylist{1}    = {'2020-05-25_15-05'}; %liste de tous les fichiers, tous les protocoles
config{5}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Puff'};
config{5}.LFP.rename          = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{5}.LFP.origin_WoD      = {'E11', 'E11'};
config{5}.LFP.chan_depth      = {73, 323, 573, 823, 1073, 1323, 1573, 1823};


config{5}.circus.channel      = {'E04', 'E05','E06', 'E10', 'E11', 'E13'};       %name of the first electrode
config{5}.circus.rename      = {'E4', 'E5','E6', 'E10', 'E11', 'E13'};       %name of the first electrode


%% rat 6
config{6}                     = configcommon;
config{6}.concatdata_path     = concatdata_path;
config{6}.statsavedir         = statsavedir;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.imagesavedir_data   = imagesavedir_data;
config{6}.prefix              = 'Rat-27_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'2020_05_27_WOD');                       %path to patient data

config{6}.directorylist{1}    = {'2020-05-27_14-19', '2020-05-27_15-40', '2020-05-27_16-19'}; %liste de tous les fichiers, tous les protocoles
config{6}.LFP.channel         = { 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'Puff'};
config{6}.LFP.rename         = { 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E0'};
config{6}.LFP.chan_depth         = { 410,660,910, 1160, 1410, 1660, 1910, 2160};
config{6}.LFP.origin_WoD          = {'E8', 'E9'};
config{6}.LFP.origin_WoR          = {'E8', 'E12'};


config{6}.circus.channel      = {'E02','E03', 'E06', 'E08', 'E10', 'E11'};       %name of the first electrode
config{6}.circus.rename      = {'E2','E3', 'E6', 'E8', 'E10', 'E11'};       %name of the first electrode

%% rat 7
config{7}                     = configcommon;
config{7}.concatdata_path     = concatdata_path;
config{7}.statsavedir         = statsavedir;
config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.imagesavedir_data   = imagesavedir_data;
config{7}.prefix              = 'Rat-16_06_2020-';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'2020_06_16_WOD');                       %path to patient data

config{7}.directorylist{1}    = {'2020-06-16_14-35', '2020-06-16_15-57', '2020-06-16_16-35'}; %liste de tous les fichiers, tous les protocoles
config{7}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Puff'};
config{7}.LFP.rename          = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{7}.LFP.chan_depth        = {85, 335, 585, 835, 1085, 1335, 1585, 1835};
config{7}.LFP.origin_WoD          = {'E10', 'E10'};
config{7}.LFP.origin_WoR          = {'E13', 'E13'};


config{7}.circus.channel      = {'E03','E05', 'E11', 'E12', 'E14'};       %name of the first electrode
config{7}.circus.rename     = {'E3','E5', 'E11', 'E12', 'E14'};       %name of the first electrode

%% rat 8
config{8}                     = configcommon;
config{8}.concatdata_path     = concatdata_path;
config{8}.statsavedir         = statsavedir;
config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.imagesavedir_data   = imagesavedir_data;
config{8}.prefix              = 'Rat-22_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'2020_07_22_WOD2');                       %path to patient data

config{8}.directorylist{1}    = {'2020-07-22_13-28'}; %liste de tous les fichiers, tous les protocoles
config{8}.LFP.channel         = { 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{8}.LFP.rename          = {'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{8}.LFP.chan_depth         = {337, 587, 837, 1087, 1337, 1587, 1837};

config{8}.LFP.origin_WoD          = {'E10', 'E9'};
config{8}.LFP.origin_WoR          = {'E13', 'E13'};


config{8}.circus.channel      = {'E02', 'E05', 'E06', 'E11', 'E12'};       %name of the first electrode
config{8}.circus.rename    = {'E2', 'E5', 'E6', 'E11', 'E12'};       %name of the first electrode

%% rat 9
config{9}                     = configcommon;
config{9}.concatdata_path     = concatdata_path;
config{9}.statsavedir         = statsavedir;
config{9}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{9}.imagesavedir        = imagesavedir;
config{9}.imagesavedir_data   = imagesavedir_data;
config{9}.prefix              = 'Rat-28_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{9}.rawdir              = fullfile(rootpath_data,'2020_07_28_WOD');                       %path to patient data

config{9}.directorylist{1}    = {'2020-07-28_13-54', '2020-07-28_15-35', '2020-07-28_16-48', '2020-07-28_17-35'}; %liste de tous les fichiers, tous les protocoles
config{9}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP',  'Respi'};
config{9}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{9}.LFP.chan_depth         = {41, 291, 541, 791, 1041, 1291, 1541, 1791};

config{9}.LFP.origin_WoD          = {'E9', 'E9'};
config{9}.LFP.origin_WoR          = {'E12', 'E10'};


config{9}.circus.channel      = {'E07', 'E08', 'E09', 'E10', 'E11', 'E15'};       %name of the first electrode
config{9}.circus.rename      = {'E7', 'E8', 'E9', 'E10', 'E11', 'E15'};       %name of the first electrode

%% rat 10
config{10}                     = configcommon;
config{10}.concatdata_path     = concatdata_path;
config{10}.statsavedir         = statsavedir;
config{10}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{10}.imagesavedir        = imagesavedir;
config{10}.imagesavedir_data   = imagesavedir_data;
config{10}.prefix              = 'Rat-31_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{10}.rawdir              = fullfile(rootpath_data,'2020_07_31_WOD');                       %path to patient data

config{10}.directorylist{1}    = {'2020-07-31_13-19', '2020-07-31_14-50', '2020-07-31_15-19'}; %liste de tous les fichiers, tous les protocoles
config{10}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{10}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{10}.LFP.chan_depth        = {67, 317, 567, 817, 1067, 1317, 1567, 1817};

config{10}.LFP.origin_WoD          = {'E13', 'E11'};
config{10}.LFP.origin_WoR          = {'E13', 'E11'};


config{10}.circus.channel      = {'E11', 'E12', 'E14', 'E15'};       %name of the first electrode
config{10}.circus.rename      = {'E11', 'E12', 'E14', 'E15'};       %name of the first electrode

%% rat 11
config{11}                     = configcommon;
config{11}.concatdata_path     = concatdata_path;
config{11}.statsavedir         = statsavedir;
config{11}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{11}.imagesavedir        = imagesavedir;
config{11}.imagesavedir_data   = imagesavedir_data;
config{11}.prefix              = 'Rat-14_08_2020-';                                                        %patient name. Must end by "-". namepatient-
config{11}.rawdir              = fullfile(rootpath_data,'2020_08_14_WOD');                       %path to patient data

config{11}.directorylist{1}    = {'2020-08-14_13-26', '2020-08-14_14-57', '2020-08-14_15-26'}; %liste de tous les fichiers, tous les protocoles
config{11}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'Respi'};
config{11}.LFP.rename        = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E0'};
config{11}.LFP.chan_depth        = {126, 376, 626, 876, 1126, 1376, 1626};

config{11}.LFP.origin_WoD          = {'E10', 'E10'};
config{11}.LFP.origin_WoR          = {'E12', 'E12'};


config{11}.circus.channel      = {'E02', 'E03', 'E05', 'E09', 'E10', 'E11', 'E14'};       %name of the first electrode
config{11}.circus.rename      = {'E2', 'E3', 'E5', 'E9', 'E10', 'E11', 'E14'};       %name of the first electrode

%% rat 12
config{12}                     = configcommon;
config{12}.concatdata_path     = concatdata_path;
config{12}.statsavedir         = statsavedir;
config{12}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{12}.imagesavedir        = imagesavedir;
config{12}.imagesavedir_data   = imagesavedir_data;

config{12}.prefix              = 'Rat-10_09_2020-';                                                        %patient name. Must end by "-". namepatient-
config{12}.rawdir              = fullfile(rootpath_data,'2020_09_10_WOD');                       %path to patient data

config{12}.directorylist{1}    = {'2020-09-10_20-01'}; %liste de tous les fichiers, tous les protocoles
config{12}.LFP.channel         = {'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{12}.LFP.rename        = { 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{12}.LFP.chan_depth        = {434, 684, 934, 1184, 1434, 1684, 1934};


config{12}.circus.channel      = { 'E06', 'E10', 'E11','E12','E13' 'E14'};       %name of the first electrode
config{12}.circus.rename      = { 'E6', 'E10', 'E11','E12','E13', 'E14'};       %name of the first electrode

%% rat 13
config{13}                     = configcommon;
config{13}.concatdata_path     = concatdata_path;
config{13}.statsavedir         = statsavedir;
config{13}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{13}.imagesavedir        = imagesavedir;
config{13}.imagesavedir_data   = imagesavedir_data;

config{13}.prefix              = 'Rat-15_09_2020-';                                                        %patient name. Must end by "-". namepatient-
config{13}.rawdir              = fullfile(rootpath_data,'2020_09_15_WOD_2');                       %path to patient data

config{13}.directorylist{1}    = {'2020-09-15_16-49','2020-09-15_18-15'}; %liste de tous les fichiers, tous les protocoles
config{13}.LFP.channel         = {'E16LFP','E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{13}.LFP.rename        = { 'E16','E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{13}.LFP.chan_depth        = {25,275, 525, 775, 1025, 1275, 1525, 1775};

config{13}.circus.channel      = { 'E02','E04','E05','E07','E09', 'E10', 'E11','E12','E13' 'E14'};       %name of the first electrode
config{13}.circus.rename      = { 'E2','E4','E5','E7','E9', 'E10', 'E11','E12','E13', 'E14'};       %name of the first electrode

%% rat 14
config{14}                     = configcommon;
config{14}.concatdata_path     = concatdata_path;
config{14}.statsavedir         = statsavedir;
config{14}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{14}.imagesavedir        = imagesavedir;
config{14}.imagesavedir_data   = imagesavedir_data;

config{14}.prefix              = 'Rat-16_09_2020-';                                                        %patient name. Must end by "-". namepatient-
config{14}.rawdir              = fullfile(rootpath_data,'2020_09_16_WOD');                       %path to patient data

config{14}.directorylist{1}    = {'2020-09-16_17-44'}; %liste de tous les fichiers, tous les protocoles
config{14}.LFP.channel         = {'E16LFP','E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{14}.LFP.rename        = { 'E16','E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{14}.LFP.chan_depth         = {119,369, 619, 869,1119, 1369, 1619, 1869};


config{14}.circus.channel      = { 'E02','E03','E05','E07','E14' 'E15'};       %name of the first electrode
config{14}.circus.rename      = { 'E2','E3','E5','E7','E14', 'E15'};       %name of the first electrode

%% rat 15
config{15}                     = configcommon;
config{15}.concatdata_path     = concatdata_path;
config{15}.statsavedir         = statsavedir;
config{15}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{15}.imagesavedir        = imagesavedir;
config{15}.imagesavedir_data   = imagesavedir_data;
config{15}.epoch.toi.WoD             = [-900, 3200];
config{15}.timefreq.toi.long           = [-900 3200];

config{15}.prefix              = 'Rat-23_09_2020-';                                                        %patient name. Must end by "-". namepatient-
config{15}.rawdir              = fullfile(rootpath_data,'2020_09_23_WOD');                       %path to patient data

config{15}.directorylist{1}    = {'2020-09-23_17-56','2020-09-23_19-05','2020-09-23_22-06'}; %liste de tous les fichiers, tous les protocoles
config{15}.LFP.channel         = {'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'Respi'};
config{15}.LFP.rename        = {'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E0'};
config{15}.LFP.chan_depth        = {516, 766, 1016, 1266, 1516,1766 };


config{15}.circus.channel      = { 'E05','E08','E10','E11' 'E13'};       %name of the first electrode
config{15}.circus.rename      = { 'E5','E8','E10','E11', 'E13'};       %name of the first electrode

%% rat 16
config{16}                     = configcommon;
config{16}.concatdata_path     = concatdata_path;
config{16}.statsavedir         = statsavedir;
config{16}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{16}.imagesavedir        = imagesavedir;
config{16}.imagesavedir_data   = imagesavedir_data;

config{16}.prefix              = 'Rat-29_09_2020-';                                                        %patient name. Must end by "-". namepatient-
config{16}.rawdir              = fullfile(rootpath_data,'2020_09_29_WOD');                       %path to patient data

config{16}.directorylist{1}    = {'2020-09-29_15-38','2020-09-29_18-45'}; %liste de tous les fichiers, tous les protocoles
config{16}.LFP.channel         = {'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP','E09LFP', 'Respi'};
config{16}.LFP.rename        = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E0'};
config{16}.LFP.chan_depth        = {238, 488, 738, 988, 1238, 1488,1738};


config{16}.circus.channel      = { 'E02','E03','E04','E05','E06','E10','E11','E13' 'E14'};       %name of the first electrode
config{16}.circus.rename      = { 'E3','E4','E5','E6','E7','E11','E12','E14', 'E15'};       %name of the first electrode

%% rat 17
config{17}                     = configcommon;
config{17}.concatdata_path     = concatdata_path;
config{17}.statsavedir         = statsavedir;
config{17}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{17}.imagesavedir        = imagesavedir;
config{17}.imagesavedir_data   = imagesavedir_data;

config{17}.prefix              = 'Rat-07_10_2020-';                                                        %patient name. Must end by "-". namepatient-
config{17}.rawdir              = fullfile(rootpath_data,'2020_10_07_WOD');                       %path to patient data

config{17}.directorylist{1}    = {'2020-10-07_19-15'}; %liste de tous les fichiers, tous les protocoles
config{17}.LFP.channel         = { 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP','E09LFP', 'Respi'};
config{17}.LFP.rename        = { 'E15', 'E14', 'E13', 'E12','E11','E10', 'E0'};
config{17}.LFP.chan_depth        = { 432, 682, 932, 1182, 1432,1682};


config{17}.circus.channel      = { 'E01','E04','E05','E07','E08','E10','E11','E12' 'E13'};       %name of the first electrode
config{17}.circus.rename      = { 'E1','E4','E5','E7','E8','E10','E11','E12', 'E13'};       %name of the first electrode



