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

%% config common for all patients
configcommon.name                  = {'WoD'};
configcommon.circus.writedeadfile  = 'no';

%% patient 1
config{1}                     = configcommon;
config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.prefix              = 'Rat-27_05_2020-';                                                        %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'2020_05_27_WOD');                       %path to patient data
config{1}.directorylist{1}    = {'2020-05-27_14-19','2020-05-27_15-40','2020-05-27_16-19'};                                               %list of folders to analyse
config{1}.circus.channel      = {'E10','E08','E06','E03','E02'};       %name of the first electrode
config{1}.circus.reref        = 'no';
config{1}.circus.refchan      = [];
config{1}.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
config{1}.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
config{1}.circus.hpfreq       = 0; % even when not using
config{1}.circus.postfix      = []; % after using circus-gui-matlab's SAVE number
