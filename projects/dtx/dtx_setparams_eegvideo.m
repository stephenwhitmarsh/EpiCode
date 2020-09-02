
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

function [config] = dtx_setparams_eegvideo(config)

disp('setting parameters for EEG-video');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-EEG-VIDEO/';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-EEGRodents-Brainvision/';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-EEG-VIDEO';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-EEGRodents-Brainvision\';
    os                  = 'windows';
else
    error('Platform not supported')
end


%% Congig common for all rats

datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis);

configcommon.os                        = os;
% SlowWave_begin
% SlowWave_peak
configcommon.name                      = {'SlowWave_begin','SlowWave_EMGalign','SlowWavealign_EMG'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.continous                 = true;

configcommon.muse.startend             = {'SlowWave','SlowWave';'SlowWave_EMG','SlowWave_EMG'; 'SlowWave','SlowWave'};   % start and end Muse marker. For defining trials
configcommon.labels.macro              = {'M1G','M1D','PtA'};%do not put the emg channels here
%1, 2 and 3 associated with config.muse.startend
configcommon.epoch.toi{1}              = [-5, 25];
configcommon.epoch.toi{2}              = [-5, 25];
configcommon.epoch.toi{3}              = [-5, 25];
configcommon.epoch.toi{4}              = [-5, 25];
configcommon.epoch.pad{1}              = 10;
configcommon.epoch.pad{2}              = 10;
configcommon.epoch.pad{3}              = 10;
configcommon.epoch.pad{4}              = 10;

configcommon.align.name                = {'SlowWave_peak'};%{'SlowWave_R','SlowWave_R'};
configcommon.align.channel             = {'M1G','M1G'};      % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
configcommon.align.flip                = {'no','no'};
configcommon.align.abs                 = {'no','no'};
configcommon.align.method              = {'nearestmin','nearestmin'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp','lp'};%bp pour retirer les variations trop lentes de signal qui empêchent de comparer la baseline et de chercher des maxima locaux
configcommon.align.freq                = {5,5};          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no','no'};
configcommon.align.thresh.value        = [1, 1];
configcommon.align.thresh.method       = {'trial','trial'};%'medianbl','both';
configcommon.align.toiplot             = {[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.toiactive           = {[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline         = {[-10.5, -0.5],[-10.5, -0.5],[-10.5, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.reref               = 'no';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = 'all';
configcommon.align.notch               = 'no';
configcommon.align.demean              = 'yes';
%for detection of begin of event
configcommon.align.begin.name          = {'SlowWave_begin'};
configcommon.align.begin.maxtimeshift  = 0.4;
configcommon.align.begin.thresh        = 0.2; % percent of peak

configcommon.alignEMG.name             = {'SlowWave_EMGalign'}; %name of analysis
configcommon.alignEMG.channel          = {'EMG1'};
configcommon.alignEMG.toiplot          = {[-1,  2],[-1,  2]}; 
configcommon.alignEMG.toiactive        = {[-0.5, 0.5], [-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.alignEMG.toibaseline      = {[-1, -0.5], [-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.alignEMG.maxtimeshift     = 0.25; %if abs(timeshift) is superior, trial is considered as misdetected and remived 
configcommon.alignEMG.hpfreq           = 100;


configcommon.LFP.name                  = configcommon.name;
configcommon.LFP.emg                   = {'no','EMG1','EMG1'};%name of EMG channel associated with marker LFP.name. 'no' if no EMG associated to this seizure side 
configcommon.LFP.emgmarker             = {'no','no','SlowWave_EMG','SlowWave_EMG','SlowWave_EMG','SlowWave_EMG'};
configcommon.LFP.electrodetoplot       = {'M1G', 'M1G', 'M1G', 'M1G'};
configcommon.LFP.motorcortex           = {'M1G', 'M1G', 'M1G', 'M1G'};

configcommon.LFP.flip                  = true;
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 1;
configcommon.LFP.hpfiltord             = 4;
configcommon.LFP.hpfilttype            = 'but';
configcommon.LFP.resamplefs            = 1024;
configcommon.LFP.baseline              = 'no';
configcommon.LFP.baselinewindow        = {[-2, -1], [-2, -1], [-2, -1], [-2, -1]};
configcommon.LFP.slidestep             = 0.01;
configcommon.LFP.reref                 = 'no';
configcommon.LFP.rerefmethod           = 'avg';
configcommon.LFP.refchannel            = 'all';
configcommon.LFP.bsfilter              = 'no';
configcommon.LFP.bsfreq                = [1 30];
configcommon.LFP.lpfilter              = 'yes';
configcommon.LFP.lpfreq                = 30;
configcommon.LFP.lpfilttype            = 'fir';

configcommon.EMG.reref                 = 'yes';
configcommon.EMG.rerefmethod           = 'bipolar';
configcommon.EMG.refchannel            = 'EMG2';
configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 100;
configcommon.EMG.bsfilter              = 'yes';
configcommon.EMG.bsfreq                = [49 51];
configcommon.EMG.bsfiltord             = 3;
configcommon.EMG.envmethod             = 'rms';
configcommon.EMG.envparam              = 30;
configcommon.EMG.toi                   = [-5 5];

configcommon.LFP.TFR.doTFR                    = true;
configcommon.LFP.TFR.toi                      = [-15:0.01:35];
configcommon.LFP.TFR.baseline                 = [-10 -5];
configcommon.LFP.TFR.baselinetype             = 'relchange';

configcommon.stats.toibaseline                  = {[-10.5, -0.5],[-10.5, -0.5],[-10.5, -0.5]};%{[-4 -1],[-4 -1],[-4 -1],[-4 -1]};
configcommon.stats.alpha                        = 0.025;
configcommon.stats.toiploteeg                   = [-5 5];
configcommon.stats.toilatency                   = [-0.2 0.2];
configcommon.topoplot.timestep                  = 0.1;
configcommon.topoplot.toi                       = [-0.5 0.5];

%% Rodent 1
%NOTER QUE LES EMG SONT CHOISIS : REJETES SI PAS DE REPONSE EMG OU SI ARTEFACTS
config{1}                           = configcommon;
config{1}.prefix                    = 'DTX-EEGEMG-12-';
config{1}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-12');
config{1}.rawlabels.oldnames        = {'47','45','48','46','49'}; %for conversion deltamed to brainvision
config{1}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{1}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-12');       % where to print images

config{1}.injectiontime             = datetime('19-Feb-2020 10:30:00');

%config{irat}.directorylist
filelist = dir(config{1}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{1}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% Rodent 2 
config{2}                           = configcommon;
config{2}.prefix                    = 'DTX-EEGEMG-14-';
config{2}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-14');
config{2}.rawlabels.oldnames        = {'57','55','56','59'}; %for conversion deltamed to brainvision
config{2}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{2}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-14');       % where to print images

config{2}.injectiontime             = datetime('19-Feb-2020 15:20:00');

%config{irat}.directorylist
filelist = dir(config{2}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{2}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% Rodent 3
config{3}                           = configcommon;
config{3}.prefix                    = 'DTX-EEGEMG-16-';
config{3}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-16');
config{3}.rawlabels.oldnames        = {'52','50','51','54'}; %for conversion deltamed to brainvision
config{3}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{3}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-16');       % where to print images

config{3}.injectiontime             = datetime('03-Jun-2020 16:15:00');

%config{irat}.directorylist
filelist = dir(config{3}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{3}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% Rodent 4 - DC après 24h
% début à 17h00, 6h post injection, car changement de position à ce moment
% là. Pas enregistré en DC
config{4}                           = configcommon;
config{4}.prefix                    = 'DTX-EEGEMG-18-';
config{4}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-18');
config{4}.rawlabels.oldnames        = {'52','50','51','54'}; %for conversion deltamed to brainvision
config{4}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-18');       % where to print images

config{4}.injectiontime             = datetime('09-Jun-2020 11:30:00');

%config{irat}.directorylist
filelist = dir(config{4}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{4}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% FUSIONNER DC ET DELTAMED

%% Rodent 5 : Début Deltamed, fin sur ampli DC après 24h
%seulement la fin (DC) est convertie. voir la conversion de la première partie plus tard.
config{5}                           = configcommon;
config{5}.prefix                    = 'DTX-EEGEMG-19-Deltamed-';
config{5}.rawdir                    = 'W:\teamcharpier\Projet_P_Baudin\DTX\Rat\20200609-EEG-EMG-18-19-20-21';%fullfile(rootpath_data,'DTX-EEGEMG-19-Deltamed-');
config{5}.rawlabels.oldnames        = {'52','50','51','54'}; %for conversion deltamed to brainvision
config{5}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-19');       % where to print images

config{5}.injectiontime             = datetime('09-Jun-2020 14:45:00');

%config{irat}.directorylist
filelist = dir(config{5}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{5}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% Rodent 5 : Début Deltamed, fin sur ampli DC après 24h
%seulement la fin (DC) est convertie. voir la conversion de la première partie plus tard.
config{6}                           = configcommon;
config{6}.prefix                    = 'DTX-EEGEMG-19-DC-';
config{6}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-19-DC-');
config{6}.rawlabels.oldnames        = {'1-S1L','1-S1R','1-M1L','1-NTS'}; %for conversion deltamed to brainvision
config{6}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{6}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-19');       % where to print images

config{6}.injectiontime             = datetime('09-Jun-2020 14:45:00');

%config{irat}.directorylist
filelist = dir(config{6}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{6}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist

%% Rodent 7 
config{7}                           = configcommon;
config{7}.prefix                    = 'DTX-EEGEMG-20-DC-';
config{7}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-CHECK_ASUPPR_20');
config{7}.rawlabels.oldnames        = {'3-S1L','3-S1R','3-M1L','3-NTS'}; %for conversion deltamed to brainvision
config{7}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{7}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-20');       % where to print images

config{7}.injectiontime             = datetime('09-Jun-2020 16:15:00');

%config{irat}.directorylist
filelist = dir(config{7}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{7}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist


%% Rodent 7 
config{8}                           = configcommon;
config{8}.prefix                    = 'DTX-EEGEMG-20-Deltamed-';
config{8}.rawdir                    = 'W:\teamcharpier\Projet_P_Baudin\DTX\Rat\20200609-EEG-EMG-18-19-20-21';%fullfile(rootpath_data,'DTX-EEGEMG-CHECK_ASUPPR_20');
config{8}.rawlabels.oldnames        = {'57','55','56','59'}; %for conversion deltamed to brainvision
config{8}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
config{8}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-CHECK_ASUPPR_20');       % where to print images

config{8}.injectiontime             = datetime('09-Jun-2020 16:15:00');

%config{irat}.directorylist
filelist = dir(config{8}.rawdir);
i=0;
for ifile = 1:length(filelist)
    [~,~,file_extension] = fileparts(filelist(ifile).name);
    if strncmp(file_extension,'.eeg',4)
        i=i+1;
        config{8}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
    end
end
clear filelist
% A REMETTRE
% %% Rodent 6 AMPLI DC
% config{6}                           = configcommon;
% config{6}.prefix                    = 'DTX-EEGEMG-21-';
% config{6}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-21');
% config{6}.rawlabels.oldnames        = {'2-S1L','2-S1R','2-M1L','2-NTS'}; %for conversion deltamed to brainvision
% config{6}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
% config{6}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-24');       % where to print images
% 
% config{6}.injectiontime             = datetime('10-Jun-2020 15:00:00');
% 
% %config{irat}.directorylist
% filelist = dir(config{6}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{6}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist
% 
% 
% %% Rodent 7 AMPLI DC
% %commencé à marqué plus tard (pb d'enregistrement)
% config{7}                           = configcommon;
% config{7}.prefix                    = 'DTX-EEGEMG-22-';
% config{7}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-22');
% config{7}.rawlabels.oldnames        = {'1-S1L','1-S1R','1-M1L','1-NTS'}; %for conversion deltamed to brainvision
% config{7}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
% config{7}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-22');       % where to print images
% 
% config{7}.injectiontime             = datetime('29-Jun-2020 11:15:00');
% 
% %config{irat}.directorylist
% filelist = dir(config{7}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{7}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist
% 
% %% Rodent 8 AMPLI DC
% config{8}                           = configcommon;
% config{8}.prefix                    = 'DTX-EEGEMG-23-';
% config{8}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-23');
% config{8}.rawlabels.oldnames        = {'2-S1L','2-S1R','2-M1L','2-NTS'}; %for conversion deltamed to brainvision
% config{8}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
% config{8}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-23');       % where to print images
% 
% config{8}.injectiontime             = datetime('29-Jun-2020 14:55:00');
% 
% %config{irat}.directorylist
% filelist = dir(config{8}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{8}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist


% %% Rodent 7 TEST, A SUPPR AMPLI DC
% config{7}                           = configcommon;
% config{7}.prefix                    = 'DTX-EEGEMG-CHECK_ASUPPR_20-';
% config{7}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-CHECK_ASUPPR_20');
% config{7}.rawlabels.oldnames        = {'3-S1L','3-S1R','3-M1L','3-NTS'}; %for conversion deltamed to brainvision
% config{7}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
% config{7}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-CHECK_ASUPPR_20');       % where to print images
% 
% config{7}.injectiontime             = datetime('09-Jun-2020 16:15:00');
% 
% %config{irat}.directorylist
% filelist = dir(config{7}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{7}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist


% %% Rodent 4 AMPLI DC
% config{4}                           = configcommon;
% config{4}.prefix                    = 'DTX-EEGEMG-24-';
% config{4}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-24');
% config{4}.rawlabels.oldnames        = {'3-S1L','3-S1R','3-M1L','3-NTS'}; %for conversion deltamed to brainvision
% config{4}.rawlabels.newnames        = {'M1G','M1D', 'EMG1', 'EMG2'};
% config{4}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-24');       % where to print images
% 
% config{4}.injectiontime             = datetime('29-Jun-2020 16:45:00');
% 
% %config{irat}.directorylist
% filelist = dir(config{4}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{4}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist
% end


