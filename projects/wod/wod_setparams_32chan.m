function [config] = wod_setparams_32chan

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Sofia';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/32_ch/Extra_Neuralynx';
    rootpath_concatdata = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/32_ch/concatenated_LFP';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Sofia';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\32_ch\Extra_Neuralynx';
    rootpath_concatdata = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\32_ch\concatenated_LFP';
    os                  = 'windows';
else
    error('Platform not supported')
end 

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');
imagesavedir_data= {fullfile(imagesavedir,'TFR'), fullfile(imagesavedir,'band_depth'),fullfile(imagesavedir,'band_cx')};
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all patients
configcommon.concatdata_path = rootpath_concatdata;
configcommon.muse.templatemarker   = fullfile(datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

configcommon.name                  = {'WoD'};
configcommon.LFP.allchannel        = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01','Events_0001'};
configcommon.LFP.name              = configcommon.name;

configcommon.muse.backupdir            = fullfile(datasavedir,'Backup_MuseMarker');
configcommon.muse.startmarker.WoD      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD             = [-900, 3400];
configcommon.epoch.pad.WoD             = 5;
configcommon.LFP.resamplefs                     = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                          = true; %save computed data to disk

configcommon.LFP.lpfilter_wod_detection         = 7;%Hz
configcommon.LFP.wod_toisearch                  = [-1 50]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [-1 25]; %s, were to find wor positive peak relative to the muse marker
configcommon.LFP.hpfilter_wod_exclusion         = 1; %Hz

configcommon.timefreq.foi          = 1:2:100;%[1:2:100] is ritght value
configcommon.timefreq.foi_band     = {[1 5],[10 20],[25 50],[70 90]};%Hz
configcommon.timefreq.t_ftimwin.long    = 10;% in second, length of the time window 10 (right value)
configcommon.timefreq.t_ftimwin.short = 4; % in seconds
configcommon.timefreq.timestep.long     = 2.5;% in second, time between 2 sliding time windows. can be 'all' 2.5 right value
configcommon.timefreq.timestep.short     = 1;% in second, time between 2 sliding time windows. can be 'all'
configcommon.timefreq.movmeanwin   = [1,1,1,100,100,100];%in sample points, one value per analysis_name
configcommon.timefreq.tapsmofrq.long    = 0; %2
configcommon.timefreq.tapsmofrq.short    = 0; %2
configcommon.timefreq.toi.long           = [-900 3400];
configcommon.timefreq.toi.short           = [-900 600];


configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Sofia.params');
configcommon.circus.writedeadfile  = 'no';
configcommon.circus.reref        = 'no';
configcommon.circus.refchan      = [];
configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.hpfreq       = 0; % even when not using
configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number



%% rat 1
config{1}                     = configcommon;
config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.imagesavedir_data   = imagesavedir_data;
config{1}.prefix              = 'Rat-2021_03_12';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'2021_03_12_WOD');                       %path to patient data

config{1}.directorylist{1}    = {'2021-03-12_14-01', '2021-03-12_15-36', '2021-03-12_17-05'}; %liste de tous les fichiers, tous les protocoles
%config{1}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{1}.LFP.channel         = {'E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{1}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02'};
config{1}.LFP.chan_depth      = {3138,	3038,	2938,	2838,	2738,	2638,	2538,	2438,	2338,	2238,	2138,	2038,	1938,	1838,	1738,	1638,	1538,	1438,	1338,	1238,	1138,	1038,	938,	838,	738,	638,	538,	438,	338,	238,	138};
config{1}.LFP.origin_WoD          = {'E14', 'E14'};
config{1}.LFP.origin_WoR          = {'E11', 'E11'};

config{1}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{1}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode


%% rat 2
config{2}                     = configcommon;
config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.imagesavedir_data   = imagesavedir_data;
config{2}.prefix              = 'Rat-2021_03_30';                                              %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'2021_03_30_WOD');                       %path to patient data

config{2}.directorylist{1}    = {'2021-03-30_14-10', '2021-03-30_15-22'}; %liste de tous les fichiers, tous les protocoles
%config{2}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{2}.LFP.channel         = { 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{2}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04'};
config{2}.LFP.chan_depth      = {2897.2,	2797.2,	2697.2,	2597.2,	2497.2,	2397.2,	2297.2,	2197.2,	2097.2,	1997.2,	1897.2,	1797.2,	1697.2,	1597.2,	1497.2,	1397.2,	1297.2,	1197.2,	1097.2,	997.2,	897.2,	797.2,	697.2,	597.2,	497.2,	397.2,	297.2,	197.2,	97.1999999999998};
config{2}.LFP.origin_WoD          = {'E14', 'E14'};
config{2}.LFP.origin_WoR          = {'E11', 'E11'};

config{2}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{2}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 3
config{3}                     = configcommon;
config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.imagesavedir_data   = imagesavedir_data;
config{3}.prefix              = 'Rat-2021_03_10';                                              %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'2021_03_10_WOD');                       %path to patient data

config{3}.directorylist{1}    = {'2021-03-10_15-53', '2021-03-10_17-25','2021-03-10_18-48'}; %liste de tous les fichiers, tous les protocoles
%config{3}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{3}.LFP.channel         = { 'E32','E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02'};
config{3}.LFP.rename          = {'E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{3}.LFP.chan_depth      = {3215,	3115,	3015,	2915,	2815,	2715,	2615,	2515,	2415,	2315,	2215,	2115,	2015,	1915,	1815,	1715,	1615,	1515,	1415,	1315,	1215,	1115,	1015,	915,	815,	715,	615,	515,	415,	315,	215};
config{3}.LFP.origin_WoD          = {'E14', 'E14'};
config{3}.LFP.origin_WoR          = {'E11', 'E11'};

config{3}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{3}.circus.rename       = {'E2', 'E4', 'E11'};

%% rat 4
config{4}                     = configcommon;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.imagesavedir_data   = imagesavedir_data;
config{4}.prefix              = 'Rat-2021_03_05';                                              %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'2021_03_05_WOD');                       %path to patient data

config{4}.directorylist{1}    = {'2021-03-05_15-20'}; %liste de tous les fichiers, tous les protocoles
%config{4}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{4}.LFP.channel         = { 'E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{4}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02'};
config{4}.LFP.chan_depth      = {3103,	3003,	2903,	2803,	2703,	2603,	2503,	2403,	2303,	2203,	2103,	2003,	1903,	1803,	1703,	1603,	1503,	1403,	1303,	1203,	1103,	1003,	903,	803,	703,	603,	503,	403,	303,	203,	103};
config{4}.LFP.origin_WoD          = {'E14', 'E14'};
config{4}.LFP.origin_WoR          = {'E11', 'E11'};

config{4}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{4}.circus.rename       = {'E2', 'E4', 'E11'};

%% rat 5
config{5}                     = configcommon;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.imagesavedir_data   = imagesavedir_data;
config{5}.prefix              = 'Rat-2021_02_18';                                              %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'2021_02_18_WOD');                       %path to patient data

config{5}.directorylist{1}    = {'2021-02-18_15-54','2021-02-18_17-54'}; %liste de tous les fichiers, tous les protocoles
%config{5}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{5}.LFP.channel         = { 'E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{5}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02'};
config{5}.LFP.chan_depth      = {3119.3,	3019.3,	2919.3,	2819.3,	2719.3,	2619.3,	2519.3,	2419.3,	2319.3,	2219.3,	2119.3,	2019.3,	1919.3,	1819.3,	1719.3,	1619.3,	1519.3,	1419.3,	1319.3,	1219.3,	1119.3,	1019.3,	919.3,	819.3,	719.3,	619.3,	519.3,	419.3,	319.3,	219.3,	119.3,};
config{5}.LFP.origin_WoD          = {'E14', 'E14'};
config{5}.LFP.origin_WoR          = {'E11', 'E11'};
config{5}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{5}.circus.rename       = {'E2', 'E4', 'E11'};

%% rat 6
config{6}                     = configcommon;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.imagesavedir_data   = imagesavedir_data;
config{6}.prefix              = 'Rat-2021_02_10';                                              %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'2021_02_10_WOD');                       %path to patient data
config{6}.directorylist{1}    = {'2021-02-10_15-22','2021-02-10_17-22','2021-02-10_18-34'}; %liste de tous les fichiers, tous les protocoles
%config{6}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{6}.LFP.channel         = { 'E32','E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{6}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
config{6}.LFP.chan_depth      = {3247.7,	3147.7,	3047.7,	2947.7,	2847.7,	2747.7,	2647.7,	2547.7,	2447.7,	2347.7,	2247.7,	2147.7,	2047.7,	1947.7,	1847.7,	1747.7,	1647.7,	1547.7,	1447.7,	1347.7,	1247.7,	1147.7,	1047.7,	947.7,	847.7,	747.7,	647.7,	547.7,	447.7,	347.7,	247.7,	147.7};
config{6}.LFP.origin_WoD          = {'E14', 'E14'};
config{6}.LFP.origin_WoR          = {'E11', 'E11'};

config{6}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{6}.circus.rename       = {'E2', 'E4', 'E11'};


