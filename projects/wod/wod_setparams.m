function [config] = wod_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-iso/DATA_Antoine/Extra_Neuralynx';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-iso\DATA_Antoine\Extra_Neuralynx';
    os                  = 'windows';
else
    error('Platform not supported')
end 

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all patients
configcommon.muse.templatemarker   = fullfile(datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

configcommon.name                  = {'WoD'};
configcommon.LFP.allchannel        = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP','Puff'};
configcommon.LFP.name              = configcommon.name;

configcommon.muse.startmarker.WoD      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD             = [-900, 3400];
configcommon.epoch.pad.WoD             = 5;
configcommon.LFP.resamplefs                     = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                          = true; %save computed data to disk

configcommon.LFP.lpfilter_wod_detection         = 10;%Hz
configcommon.LFP.wod_toisearch                  = [0 40]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [0 30]; %s, were to find wor positive peak relative to the muse marker

configcommon.timefreq.foi          = {[0.5:0.25:4.5],[4.5:0.25:10],[10:0.25:25],[25:0.25:50]};%Hz
configcommon.timefreq.t_ftimwin    = 5;% in second, length of the time window
configcommon.timefreq.timestep     = 2.5;% in second, time between 2 sliding time windows. can be 'all'
configcommon.recovery.movmeanwin   = 100;%in sample points

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
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.prefix              = 'Rat-19_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'2020_05_19_WOD');                       %path to patient data

config{4}.directorylist{1}    = {'2020-05-19_14-25', '2020-05-19_15-44', '2020-05-19_16-25'}; %liste de tous les fichiers, tous les protocoles
% config{4}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{4}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{4}.LFP.rename          = {'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3','E2', 'E1', 'E0'};

config{4}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{4}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 5
config{5}                     = configcommon;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.prefix              = 'Rat-25_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'2020_05_25_WOD');                       %path to patient data

config{5}.directorylist{1}    = {'2020-05-25_15-05'}; %liste de tous les fichiers, tous les protocoles
config{5}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Puff'};
config{5}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};


config{5}.circus.channel      = {'E04', 'E05','E06', 'E10', 'E11', 'E13'};       %name of the first electrode
config{5}.circus.rename      = {'E4', 'E5','E6', 'E10', 'E11', 'E13'};       %name of the first electrode


%% rat 6
config{6}                     = configcommon;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.prefix              = 'Rat-27_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'2020_05_27_WOD');                       %path to patient data

config{6}.directorylist{1}    = {'2020-05-27_14-19', '2020-05-27_15-40', '2020-05-27_16-19'}; %liste de tous les fichiers, tous les protocoles
config{6}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Puff'};
config{6}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{6}.circus.channel      = {'E02','E03', 'E06', 'E08', 'E10', 'E11'};       %name of the first electrode
config{6}.circus.rename      = {'E2','E3', 'E6', 'E8', 'E10', 'E11'};       %name of the first electrode

%% rat 7
config{7}                     = configcommon;
config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.prefix              = 'Rat-16_06_2020-';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'2020_06_16_WOD');                       %path to patient data

config{7}.directorylist{1}    = {'2020-06-16_14-35', '2020-06-16_15-57', '2020-06-16_16-35'}; %liste de tous les fichiers, tous les protocoles
config{7}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Puff'};
config{7}.LFP.rename          = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{7}.circus.channel      = {'E03','E05', 'E11', 'E12', 'E14'};       %name of the first electrode
config{7}.circus.rename     = {'E3','E5', 'E11', 'E12', 'E14'};       %name of the first electrode

%% rat 8
config{8}                     = configcommon;
config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.prefix              = 'Rat-22_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'2020_07_22_WOD2');                       %path to patient data

config{8}.directorylist{1}    = {'2020-07-22_13-28'}; %liste de tous les fichiers, tous les protocoles
config{8}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Respi'};
config{8}.LFP.rename          = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{8}.circus.channel      = {'E02', 'E05', 'E06', 'E11', 'E12'};       %name of the first electrode
config{8}.circus.rename    = {'E2', 'E5', 'E6', 'E11', 'E12'};       %name of the first electrode

%% rat 9
config{9}                     = configcommon;
config{9}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{9}.imagesavedir        = imagesavedir;
config{9}.prefix              = 'Rat-28_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{9}.rawdir              = fullfile(rootpath_data,'2020_07_28_WOD');                       %path to patient data

config{9}.directorylist{1}    = {'2020-07-28_13-54', '2020-07-28_15-35', '2020-07-28_16-48', '2020-07-28_17-35'}; %liste de tous les fichiers, tous les protocoles
config{9}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Respi'};
config{9}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{9}.circus.channel      = {'E07', 'E08', 'E09', 'E10', 'E11', 'E15'};       %name of the first electrode
config{9}.circus.rename      = {'E7', 'E8', 'E9', 'E10', 'E11', 'E15'};       %name of the first electrode

%% rat 10
config{10}                     = configcommon;
config{10}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{10}.imagesavedir        = imagesavedir;
config{10}.prefix              = 'Rat-31_07_2020-';                                                        %patient name. Must end by "-". namepatient-
config{10}.rawdir              = fullfile(rootpath_data,'2020_07_31_WOD');                       %path to patient data

config{10}.directorylist{1}    = {'2020-07-31_13-19', '2020-07-31_14-50', '2020-07-31_15-19'}; %liste de tous les fichiers, tous les protocoles
config{10}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Respi'};
config{10}.LFP.rename         = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{10}.circus.channel      = {'E11', 'E12', 'E14', 'E15'};       %name of the first electrode
config{10}.circus.rename      = {'E11', 'E12', 'E14', 'E15'};       %name of the first electrode

%% rat 11
config{11}                     = configcommon;
config{11}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{11}.imagesavedir        = imagesavedir;
config{11}.prefix              = 'Rat-14_08_2020-';                                                        %patient name. Must end by "-". namepatient-
config{11}.rawdir              = fullfile(rootpath_data,'2020_08_14_WOD');                       %path to patient data

config{11}.directorylist{1}    = {'2020-08-14_13-26', '2020-08-14_14-57', '2020-08-14_15-26'}; %liste de tous les fichiers, tous les protocoles
config{11}.LFP.channel         = {'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP', 'E02LFP', 'E01LFP', 'Respi'};
config{11}.LFP.rename        = {'E16', 'E15', 'E14', 'E13', 'E12', 'E11', 'E10', 'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1', 'E0'};

config{11}.circus.channel      = {'E02', 'E03', 'E05', 'E09', 'E10', 'E11', 'E14'};       %name of the first electrode
config{11}.circus.rename      = {'E2', 'E3', 'E5', 'E9', 'E10', 'E11', 'E14'};       %name of the first electrode



