function [config] = wod_setparams_lfp

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix %CED library does not work on unix
%     rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod';
%     rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-iso/DATA_Antoine/Extra_Neuralynx';
%     os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-iso\DATA_Antoine\Extra_Neuralynx';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');

%% config common for all patients
configcommon.name                  = {'WoD'};
configcommon.datasavedir         = datasavedir;       %path where to save MuseStruct data
configcommon.imagesavedir        = imagesavedir;

%% rat 1
config{1}                     = configcommon;

config{1}.prefix              = 'Rat-04_03_2020-'; %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'2020_05_27_WOD');                       %path to patient data
config{1}.directorylist{1}    = {'Protocol1_27052020'}; %list of folders to analyse. For CED : WITHOUT THE EXTENSION
config{1}.directorylist{2}    = {'Protocol2_27052020'}; %list of folders to analyse. For CED : WITHOUT THE EXTENSION

config{1}.LFP.name            = {'0.5_4.5','4.5_10','10_25','25_50'};

config{1}.LFP.channel{1}            = {[],[], };


