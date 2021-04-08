function [config] = DC_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Antoine';
    rootpath_data       = '/network/lustre/iss01/charpier/analyses/wod/Antoine/data/DoubleDC';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Antoine';
    rootpath_data       = '\\lexport\iss01.charpier\analyses\wod\Antoine\data\DoubleDC';
    os                  = 'windows';
else
    error('Platform not supported')
end 

datasavedir  =  rootpath_data;
imagesavedir =  fullfile(rootpath_analysis,'images','DC');
test_filtsavedir= fullfile(imagesavedir,'DC','test_filter');
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all patients

configcommon.name                  = {'WoD'};
configcommon.DC.chanlist       = {'Vm', 'DC'};
configcommon.DC.name              = configcommon.name;

configcommon.DC.resamplefs                    =     1000;%Downsampling factor i.e. keep 1/10 sample

configcommon.DC.write                          = true; %save computed data to disk

configcommon.DC.lpfilter                        = 3;%Hz


%% protocol 1
config{1}                     = configcommon;

config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.prefix              = '08_12_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'08_12_2020_WOD');                       %path to patient data

config{1}.directorylist{1}    = {'08_12_2020_P1_300_1700'}; %liste de tous les fichiers, tous les protocoles
config{1}.DC.channel         = {'Vm', 'DC'};
config{1}.DC.rename          = {'DC_sup', 'DC_dep'};
config{1}.DC.sup_dep          = 300;
config{1}.DC.dep_dep          = 1700;

%% protocol 2
config{2}                     = configcommon;
config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.prefix              = '08_12_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'08_12_2020_WOD');                       %path to patient data

config{2}.directorylist{1}    = {'08_12_2020_P2_300_1200'}; %liste de tous les fichiers, tous les protocoles
config{2}.DC.channel         = {'Vm', 'DC'};
config{2}.DC.rename          = {'DC_sup', 'DC_dep'};
config{2}.DC.sup_dep          = 300;
config{2}.DC.dep_dep          = 1200;

%% protocol 3
config{3}                     = configcommon;
config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.prefix              = '08_12_2020_p3';                                                        %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'08_12_2020_WOD');                       %path to patient data

config{3}.directorylist{1}    = {'08_12_2020_P3_300_800'}; %liste de tous les fichiers, tous les protocoles
config{3}.DC.channel         = {'Vm', 'DC'};
config{3}.DC.rename          = {'DC_sup', 'DC_dep'};
config{3}.DC.sup_dep          = 300;
config{3}.DC.dep_dep          = 800;

%% protocol 4
config{4}                     = configcommon;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.prefix              = '10_12_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'10_12_2020_WOD');                       %path to patient data

config{4}.directorylist{1}    = {'10_12_2020_WOD_p1_300_1200'}; %liste de tous les fichiers, tous les protocoles
config{4}.DC.channel         = {'Vm', 'DC'};
config{4}.DC.rename          = {'DC_sup', 'DC_dep'};
config{4}.DC.sup_dep          = 300;
config{4}.DC.dep_dep          = 1200;

%% protocol 5
config{5}                     = configcommon;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.prefix              = '10_12_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'10_12_2020_WOD');                       %path to patient data

config{5}.directorylist{1}    = {'10_12_2020_WOD_p2_300_1200'}; %liste de tous les fichiers, tous les protocoles
config{5}.DC.channel         = {'Vm', 'DC'};
config{5}.DC.rename          = {'DC_sup', 'DC_dep'};
config{5}.DC.sup_dep          = 300;
config{5}.DC.dep_dep          = 1200;

%% protocol 6
config{6}                     = configcommon;
config{6}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{6}.imagesavedir        = imagesavedir;
config{6}.prefix              = '10_12_2020_p3';                                                        %patient name. Must end by "-". namepatient-
config{6}.rawdir              = fullfile(rootpath_data,'10_12_2020_WOD');                       %path to patient data

config{6}.directorylist{1}    = {'10_12_2020_WOD_p3_300_900'}; %liste de tous les fichiers, tous les protocoles
config{6}.DC.channel         = {'Vm', 'DC'};
config{6}.DC.rename          = {'DC_sup', 'DC_dep'};
config{6}.DC.sup_dep          = 300;
config{6}.DC.dep_dep          = 900;

%% protocol 7
config{7}                     = configcommon;
config{7}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{7}.imagesavedir        = imagesavedir;
config{7}.prefix              = '10_12_2020_p4';                                                        %patient name. Must end by "-". namepatient-
config{7}.rawdir              = fullfile(rootpath_data,'10_12_2020_WOD');                       %path to patient data

config{7}.directorylist{1}    = {'10_12_2020_WOD_p4_300_900'}; %liste de tous les fichiers, tous les protocoles
config{7}.DC.channel         = {'Vm', 'DC'};
config{7}.DC.rename          = {'DC_sup', 'DC_dep'};
config{7}.DC.sup_dep          = 300;
config{7}.DC.dep_dep          = 900;

%% protocol 8
config{8}                     = configcommon;
config{8}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{8}.imagesavedir        = imagesavedir;
config{8}.prefix              = '17_12_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{8}.rawdir              = fullfile(rootpath_data,'17_12_2020_WOD');                       %path to patient data

config{8}.directorylist{1}    = {'17_12_2020_WOD-p1_270_1000'}; %liste de tous les fichiers, tous les protocoles
config{8}.DC.channel         = {'Vm', 'DC'};
config{8}.DC.rename          = {'DC_sup', 'DC_dep'};
config{8}.DC.sup_dep          = 270;
config{8}.DC.dep_dep          = 1000;

%% protocol 9
config{9}                     = configcommon;
config{9}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{9}.imagesavedir        = imagesavedir;
config{9}.prefix              = '17_12_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{9}.rawdir              = fullfile(rootpath_data,'17_12_2020_WOD');                       %path to patient data

config{9}.directorylist{1}    = {'17_12_2020_WOD-p2_330_1000'}; %liste de tous les fichiers, tous les protocoles
config{9}.DC.channel         = {'Vm', 'DC'};
config{9}.DC.rename          = {'DC_sup', 'DC_dep'};
config{9}.DC.sup_dep          = 330;
config{9}.DC.dep_dep          = 1000;

%% protocol 10
config{10}                     = configcommon;
config{10}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{10}.imagesavedir        = imagesavedir;
config{10}.prefix              = '22_12_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{10}.rawdir              = fullfile(rootpath_data,'22_12_2020_WOD');                       %path to patient data

config{10}.directorylist{1}    = {'22_12_2020_WOD_p1_330_1000'}; %liste de tous les fichiers, tous les protocoles
config{10}.DC.channel         = {'Vm', 'DC'};
config{10}.DC.rename          = {'DC_sup', 'DC_dep'};
config{10}.DC.sup_dep          = 330;
config{10}.DC.dep_dep          = 1000;

%% protocol 11
config{11}                     = configcommon;
config{11}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{11}.imagesavedir        = imagesavedir;
config{11}.prefix              = '22_12_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{11}.rawdir              = fullfile(rootpath_data,'22_12_2020_WOD');                       %path to patient data

config{11}.directorylist{1}    = {'22_12_2020_WOD_p2_330_1000'}; %liste de tous les fichiers, tous les protocoles
config{11}.DC.channel         = {'Vm', 'DC'};
config{11}.DC.rename          = {'DC_sup', 'DC_dep'};
config{11}.DC.sup_dep          = 330;
config{11}.DC.dep_dep          = 1000;

%% protocol 12
config{12}                     = configcommon;
config{12}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{12}.imagesavedir        = imagesavedir;
config{12}.prefix              = '22_12_2020_p3';                                                        %patient name. Must end by "-". namepatient-
config{12}.rawdir              = fullfile(rootpath_data,'22_12_2020_WOD');                       %path to patient data

config{12}.directorylist{1}    = {'22_12_2020_WOD_p3_330_1000'}; %liste de tous les fichiers, tous les protocoles
config{12}.DC.channel         = {'Vm', 'DC'};
config{12}.DC.rename          = {'DC_sup', 'DC_dep'};
config{12}.DC.sup_dep          = 330;
config{12}.DC.dep_dep          = 1000;

%% protocol 13
config{13}                     = configcommon;
config{13}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{13}.imagesavedir        = imagesavedir;
config{13}.prefix              = '04_03_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{13}.rawdir              = fullfile(rootpath_data,'04_03_2020_WOD');                       %path to patient data

config{13}.directorylist{1}    = {'04_03_2020_WOD_p2_1304'}; %liste de tous les fichiers, tous les protocoles
config{13}.DC.channel         = {'Vm'};
config{13}.DC.rename          = { 'DC_dep'};
config{13}.DC.sup_dep          = NaN;
config{13}.DC.dep_dep          = 1304;

%% protocol 14
config{14}                     = configcommon;
config{14}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{14}.imagesavedir        = imagesavedir;
config{14}.prefix              = '12_03_2020_p2';                                                        %patient name. Must end by "-". namepatient-
config{14}.rawdir              = fullfile(rootpath_data,'12_03_2020_WOD');                       %path to patient data

config{14}.directorylist{1}    = {'12_03_2020_WOD_p2_1200'}; %liste de tous les fichiers, tous les protocoles
config{14}.DC.channel         = {'Vm'};
config{14}.DC.rename          = { 'DC_dep'};
config{14}.DC.sup_dep          = NaN;
config{14}.DC.dep_dep          = 1200;

%% protocol 15
config{15}                     = configcommon;
config{15}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{15}.imagesavedir        = imagesavedir;
config{15}.prefix              = '10_09_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{15}.rawdir              = fullfile(rootpath_data,'10_09_2020_WOD');                       %path to patient data

config{15}.directorylist{1}    = {'10_09_2020_WOD_p1_1688'}; %liste de tous les fichiers, tous les protocoles
config{15}.DC.channel         = {'Vm'};
config{15}.DC.rename          = { 'DC_dep'};
config{15}.DC.sup_dep          = NaN;
config{15}.DC.dep_dep          = 1688;

