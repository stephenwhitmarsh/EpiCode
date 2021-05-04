function [config] = DC_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Antoine';
    rootpath_data       = '/network/lustre/iss01/charpier/analyses/wod/Antoine/data/Intra';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Antoine';
    rootpath_data       = '\\lexport\iss01.charpier\analyses\wod\Antoine\data\Intra';
    os                  = 'windows';
else
    error('Platform not supported')
end 

datasavedir  =  rootpath_data;
imagesavedir =  fullfile(rootpath_analysis,'images','Intra');
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all patients

configcommon.name                  = {'WoD'};
configcommon.Intra.chanlist       = {'Vm', 'Im','EEG-S1-L'};
configcommon.Intra.name              = configcommon.name;

configcommon.Intra.resamplefs                    =     1000;%Downsampling factor i.e. keep 1/10 sample

configcommon.Intra.write                          = true; %save computed data to disk

configcommon.Intra.lpfilter                        = 3;%Hz
configcommon.Intra.wod_toisearch                   = [-20 20];%Hz


%% protocol 1
config{1}                     = configcommon;

config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.prefix              = '03_06_2019_n1';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'03_06_2019_WOD');                       %path to patient data

config{1}.directorylist{1}    = {'03_06_2019_n1'}; %liste de tous les fichiers, tous les protocoles
config{1}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{1}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{1}.Intra.dep          = 1561;

%% protocol 2
config{2}                     = configcommon;

config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.prefix              = '03_06_2019_n2';                                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'03_06_2019_WOD');                       %path to patient data

config{2}.directorylist{1}    = {'03_06_2019_n2'}; %liste de tous les fichiers, tous les protocoles
config{2}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{2}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{2}.Intra.dep          = 2150;

%% protocol 3
config{3}                     = configcommon;

config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.prefix              = '13_06_2019_n1';                                                        %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'13_06_2019_WOD');                       %path to patient data

config{3}.directorylist{1}    = {'13_06_2019_n1'}; %liste de tous les fichiers, tous les protocoles
config{3}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{3}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{3}.Intra.dep          = 1815;

%% protocol 4
config{4}                     = configcommon;

config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.prefix              = '13_06_2019_n2';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'13_06_2019_WOD');                       %path to patient data

config{4}.directorylist{1}    = {'13_06_2019_n2'}; %liste de tous les fichiers, tous les protocoles
config{4}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{4}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{4}.Intra.dep          = 2540;

%% protocol 5
config{5}                     = configcommon;

config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.prefix              = '13_06_2019_n3';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'13_06_2019_WOD');                       %path to patient data

config{5}.directorylist{1}    = {'13_06_2019_n3'}; %liste de tous les fichiers, tous les protocoles
config{5}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{5}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{5}.Intra.dep          = 1863;

%% protocol 6
config{6}                     = configcommon;

config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.prefix              = '24_05_2019_n1';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'24_05_2019_WOD');                       %path to patient data

config{6}.directorylist{1}    = {'24_05_2019_n1'}; %liste de tous les fichiers, tous les protocoles
config{6}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{6}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{6}.Intra.dep          = 1204;

%% protocol 7
config{7}                     = configcommon;

config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.prefix              = '24_05_2019_n2';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'24_05_2019_WOD');                       %path to patient data

config{7}.directorylist{1}    = {'24_05_2019_n2'}; %liste de tous les fichiers, tous les protocoles
config{7}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{7}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{7}.Intra.dep          = 1228;

%% protocol 8
config{8}                     = configcommon;

config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.prefix              = '24_05_2019_n3';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'24_05_2019_WOD');                       %path to patient data

config{8}.directorylist{1}    = {'24_05_2019_n3'}; %liste de tous les fichiers, tous les protocoles
config{8}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{8}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{8}.Intra.dep          = 1348;

%% protocol 9
config{9}                     = configcommon;

config{9}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{9}.imagesavedir        = imagesavedir;
config{9}.prefix              = '24_05_2019_n4';                                                        %patient name. Must end by "-". namepatient-
config{9}.rawdir              = fullfile(rootpath_data,'24_05_2019_WOD');                       %path to patient data

config{9}.directorylist{1}    = {'24_05_2019_n4'}; %liste de tous les fichiers, tous les protocoles
config{9}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{9}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{9}.Intra.dep          = 2154;

%% protocol 10
config{10}                     = configcommon;

config{10}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{10}.imagesavedir        = imagesavedir;
config{10}.prefix              = '16_02_2021_n1';                                                        %patient name. Must end by "-". namepatient-
config{10}.rawdir              = fullfile(rootpath_data,'16_02_2021_WOD');                       %path to patient data

config{10}.directorylist{1}    = {'16_02_2021_n1'}; %liste de tous les fichiers, tous les protocoles
config{10}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{10}.Intra.rename          = {'Intra_sup', 'current','ECoG'};
config{10}.Intra.dep          = 924;

%% protocol 11
config{11}                     = configcommon;

config{11}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{11}.imagesavedir        = imagesavedir;
config{11}.prefix              = '19_02_2021_n1';                                                        %patient name. Must end by "-". namepatient-
config{11}.rawdir              = fullfile(rootpath_data,'19_02_2021_WOD');                       %path to patient data

config{11}.directorylist{1}    = {'19_02_2021_n1'}; %liste de tous les fichiers, tous les protocoles
config{11}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{11}.Intra.rename          = {'Intra_sup', 'current','ECoG'};
config{11}.Intra.dep          = 558;

%% protocol 12
config{12}                     = configcommon;

config{12}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{12}.imagesavedir        = imagesavedir;
config{12}.prefix              = '16_03_2021_n1';                                                        %patient name. Must end by "-". namepatient-
config{12}.rawdir              = fullfile(rootpath_data,'16_03_2021_WOD');                       %path to patient data

config{12}.directorylist{1}    = {'16_03_2021_n1'}; %liste de tous les fichiers, tous les protocoles
config{12}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{12}.Intra.rename          = {'Intra_sup', 'current','ECoG'};
config{12}.Intra.dep          = 463;

%% protocol 13
config{13}                     = configcommon;

config{13}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{13}.imagesavedir        = imagesavedir;
config{13}.prefix              = '29_09_2020_n2';                                                        %patient name. Must end by "-". namepatient-
config{13}.rawdir              = fullfile(rootpath_data,'29_09_2020_WOD');                       %path to patient data

config{13}.directorylist{1}    = {'29_09_2020_n2'}; %liste de tous les fichiers, tous les protocoles
config{13}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{13}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{13}.Intra.dep          = 2046;

%% protocol 14
config{14}                     = configcommon;

config{14}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{14}.imagesavedir        = imagesavedir;
config{14}.prefix              = '10_07_2020_n1';                                                        %patient name. Must end by "-". namepatient-
config{14}.rawdir              = fullfile(rootpath_data,'10_07_2020_WOD');                       %path to patient data

config{14}.directorylist{1}    = {'10_07_2020_n1'}; %liste de tous les fichiers, tous les protocoles
config{14}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{14}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{14}.Intra.dep          = 1465;

%% protocol 15
config{15}                     = configcommon;

config{15}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{15}.imagesavedir        = imagesavedir;
config{15}.prefix              = '13_11_2020_n1';                                                        %patient name. Must end by "-". namepatient-
config{15}.rawdir              = fullfile(rootpath_data,'13_11_2020_WOD');                       %path to patient data

config{15}.directorylist{1}    = {'13_11_2020_n1'}; %liste de tous les fichiers, tous les protocoles
config{15}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{15}.Intra.rename          = {'Intra_sup', 'current','ECoG'};
config{15}.Intra.dep          = 784;

%% protocol 16
config{16}                     = configcommon;

config{16}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{16}.imagesavedir        = imagesavedir;
config{16}.prefix              = '13_11_2020_n2';                                                        %patient name. Must end by "-". namepatient-
config{16}.rawdir              = fullfile(rootpath_data,'13_11_2020_WOD');                       %path to patient data

config{16}.directorylist{1}    = {'13_11_2020_n2'}; %liste de tous les fichiers, tous les protocoles
config{16}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{16}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{16}.Intra.dep          = 1954;