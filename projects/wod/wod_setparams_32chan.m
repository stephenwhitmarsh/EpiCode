function [config] = wod_setparams

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/32_ch/Extra_Neuralynx';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\32_ch\Extra_Neuralynx';
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

configcommon.timefreq.foi          = [1:2:100];%[1:2:100] is ritght value
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
config{1}.LFP.channel         = {'E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01','Events_0001'};
config{1}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02','Events_0001'};
config{1}.LFP.origin_WoD          = {'E14', 'E14'};
config{1}.LFP.origin_WoR          = {'E11', 'E11'};

config{1}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{1}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode





