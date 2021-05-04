function [config] = wod_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')% peut etre changer avec les bonnes paths
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Pierre';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-wod/Ouabaine/Extra_Neuralynx';
    rootpath_concatdata = '/network/lustre/iss01/charpier/raw/rat-wod/Ouabaine/concatenated_LFP';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Pierre';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-wod\Ouabaine\Extra_Neuralynx';
    rootpath_concatdata = '\\lexport\iss01.charpier\raw\rat-wod\Ouabaine\concatenated_LFP';
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

%% config common for all patients
configcommon.muse.templatemarker   = fullfile(datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

configcommon.name                  = {'WoD'};
configcommon.LFP.allchannel        = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP'};
configcommon.LFP.name              = configcommon.name;

configcommon.muse.backupdir            = fullfile(datasavedir,'Backup_MuseMarker');
configcommon.muse.startmarker.WoD      = 'AD__START__';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD        = 'AD__END__';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD             = [-30, 20];
configcommon.epoch.pad.WoD             = 5;



configcommon.LFP.resamplefs                     = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                          = true; %save computed data to disk

configcommon.LFP.lpfilter_wod_detection         = 3;%Hz
configcommon.LFP.wod_toisearch                  = [0 40]; %s, were to find AD negative peak in AD marker range

configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Antoine.params');
configcommon.circus.writedeadfile  = 'no';
configcommon.circus.reref        = 'no';
configcommon.circus.refchan      = [];
configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.hpfreq       = 0; % even when not using
configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number


%% rat 1
config{1}                     = configcommon;
config{1}.concatdata_path     = concatdata_path;
config{1}.statsavedir         = statsavedir;
config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.imagesavedir_data   = imagesavedir_data;
config{1}.prefix              = 'Rat-17_03_2021-';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'2021_03_17_OUABA');                       %path to patient data

config{1}.directorylist{1}    = {'2021-03-17_18-05', '2021-03-17_18-26'}; %liste de tous les fichiers, tous les protocoles
config{1}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP','E11LFP'};
config{1}.LFP.rename          = {'E31', 'E30', 'E29', 'E28', 'E27', 'E26', 'E25', 'E24', 'E23', 'E22', 'E21', 'E20', 'E19', 'E18', 'E17',...
    'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10'};
config{1}.LFP.chan_depth      =  {161,261,361,461,561,661,761,861,961,1061,1161,1261,1361,1461,1561,1661,1761,1861,1961,2061,2161,2261} ;                                                           % true depth channels
config{1}.LFP.inject_depth    = 853;


config{1}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{1}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode


%% rat 2
config{2}                     = configcommon;
config{2}.concatdata_path     = concatdata_path;
config{2}.statsavedir         = statsavedir;
config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.imagesavedir_data   = imagesavedir_data;
config{2}.prefix              = 'Rat-19_03_2021-';                                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'2021_03_19_OUABA');                       %path to patient data

config{2}.directorylist{1}    = {'2021-03-19_17-52'}; %liste de tous les fichiers, tous les protocoles
config{2}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP','E11LFP'};
config{2}.LFP.rename          = {'E31','E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{2}.LFP.chan_depth      = {150,250,350,450,550,650,750,850,950,1050,1150,1250,1350,1450,1550,1650,1750,1850,1950,2050,2150,2250};                                                           % true depth channels
config{2}.LFP.inject_depth    = 850;


config{2}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{2}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 3
config{3}                     = configcommon;
config{3}.concatdata_path     = concatdata_path;
config{3}.statsavedir         = statsavedir;
config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.imagesavedir_data   = imagesavedir_data;
config{3}.prefix              = 'Rat-23_03_2021-';                                                        %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'2021_03_23_OUABA');                       %path to patient data

config{3}.directorylist{1}    = {'2021-03-23_13-56', '2021-03-23_14-11', '2021-03-23_16-44'}; %liste de tous les fichiers, tous les protocoles
config{3}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP'};
config{3}.LFP.rename          = {'E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{3}.LFP.chan_depth      = {625,725,825,925,1025,1125,1225,1325,1425,1525,1625,1725,1825,1925,2025,2125,2225,2325};                                                           % true depth channels
config{3}.LFP.inject_depth    = 850;


config{3}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{3}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 4
config{4}                     = configcommon;
config{4}.concatdata_path     = concatdata_path;
config{4}.statsavedir         = statsavedir;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.imagesavedir_data   = imagesavedir_data;
config{4}.prefix              = 'Rat-02_04_2021-';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'2021_04_02_OUABA');                       %path to patient data

config{4}.directorylist{1}    = {'2021-04-02_15-23'}; %liste de tous les fichiers, tous les protocoles
config{4}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP','E11LFP'};
config{4}.LFP.rename          = {'E31','E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{4}.LFP.chan_depth      = {150,250,350,450,550,650,750,850,950,1050,1150,1250,1350,1450,1550,1650,1750,1850,1950,2050,2150,2250};                                                           % true depth channels
config{4}.LFP.inject_depth    = 900;


config{4}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{4}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 5
config{5}                     = configcommon;
config{5}.concatdata_path     = concatdata_path;
config{5}.statsavedir         = statsavedir;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.imagesavedir_data   = imagesavedir_data;
config{5}.prefix              = 'Rat-06_04_2021-';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'2021_04_06_OUABA');                       %path to patient data

config{5}.directorylist{1}    = {'2021-04-06_13-36', '2021-04-06_16-28', '2021-04-06_17-39'}; %liste de tous les fichiers, tous les protocoles
config{5}.LFP.channel         = {'E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP','E11LFP','E10LFP','E9LFP','E8LFP','E7LFP'};
config{5}.LFP.rename          = {'E32','E31','E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{5}.LFP.chan_depth      = {101,201,301,401,501,601,701,801,901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901,2001,2101,2201,2301};                                                           % true depth channels
config{5}.LFP.inject_depth    = 900;


config{5}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{5}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 6
config{6}                     = configcommon;
config{6}.concatdata_path     = concatdata_path;
config{6}.statsavedir         = statsavedir;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.imagesavedir_data   = imagesavedir_data;
config{6}.prefix              = 'Rat-08_04_2021-';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'2021_04_08_OUABA');                       %path to patient data

config{6}.directorylist{1}    = {'2021-04-08_16-50'}; %liste de tous les fichiers, tous les protocoles
config{6}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP','E11LFP'};
config{6}.LFP.rename          = {'E31','E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{6}.LFP.chan_depth      = {166,266,366,466,566,666,766,866,966,1066,1166,1266,1366,1466,1566,1666,1766,1866,1966,2066,2166,2266};                                                           % true depth channels
config{6}.LFP.inject_depth    = 2000;


config{6}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{6}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 7
config{7}                     = configcommon;
config{7}.concatdata_path     = concatdata_path;
config{7}.statsavedir         = statsavedir;
config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.imagesavedir_data   = imagesavedir_data;
config{7}.prefix              = 'Rat-13_04_2021-';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'2021_04_13_OUABA');                       %path to patient data

config{7}.directorylist{1}    = {'2021-04-13_14-52','2021-04-13_16-30'}; %liste de tous les fichiers, tous les protocoles
config{7}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP'};
config{7}.LFP.rename          = {'E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12','E11','E10'};
config{7}.LFP.chan_depth      = {306,406,506,606,706,806,906,1006,1106,1206,1306,1406,1506,1606,1706,1806,1906,2006,2106,2206,2306};                                                           % true depth channels
config{7}.LFP.inject_depth    = 2000;


config{7}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{7}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode



%% rat 8
config{8}                     = configcommon;
config{8}.concatdata_path     = concatdata_path;
config{8}.statsavedir         = statsavedir;
config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.imagesavedir_data   = imagesavedir_data;
config{8}.prefix              = 'Rat-19_04_2021-';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'2021_04_19_OUABA');                       %path to patient data

config{8}.directorylist{1}    = {'2021-04-19_13-01'}; %liste de tous les fichiers, tous les protocoles
config{8}.LFP.channel         = {'E32LFP','E31LFP','E30LFP','E29LFP','E28LFP','E27LFP','E26LFP','E25LFP','E24LFP','E23LFP','E22LFP','E21LFP','E20LFP','E19LFP','E18LFP','E17LFP','E16LFP','E15LFP','E14LFP','E13LFP','E12LFP'};
config{8}.LFP.rename          = {'E32','E31','E30','E29','E28','E27','E26','E25','E24','E23','E22','E21','E20','E19','E18','E17','E16','E15','E14','E13','E12'};
config{8}.LFP.chan_depth      = {102,202,302,402,502,602,702,802,902,1002,1102,1202,1302,1402,1502,1602,1702,1802,1902,2002,2102};                                                           % true depth channels
config{8}.LFP.inject_depth    = 2000;


config{8}.circus.channel      = {'E04', 'E04', 'E11'};       %name of the first electrode
config{8}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode
