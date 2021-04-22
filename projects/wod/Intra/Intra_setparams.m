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


%% protocol 1
config{1}                     = configcommon;

config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.prefix              = '08_12_2020_p1';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'08_12_2020_WOD');                       %path to patient data

config{1}.directorylist{1}    = {'08_12_2020_P1_300_1700'}; %liste de tous les fichiers, tous les protocoles
config{1}.Intra.channel         = {'Vm', 'Im','EEG-S1-L'};
config{1}.Intra.rename          = {'Intra_dep', 'current','ECoG'};
config{1}.Intra.dep          = 300;

