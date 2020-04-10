
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

function [config] = dtx_setparams_eegvideo(config)

disp('setting parameters');

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
configcommon.name                      = {'SlowWave','SlowWave_EMG','InterIctal'};
configcommon.datasavedir               = datasavedir;         % where to write data
configcommon.continous                 = true;

configcommon.muse.startend             = {'SlowWave','SlowWave';'SlowWave_EMG','SlowWave_EMG'; 'SlowWave','SlowWave'};   % start and end Muse marker. For defining trials
configcommon.labels.macro              = {'M1G','M1D','PtA'};%do not put the emg channels here
%1, 2 and 3 associated with config.muse.startend
configcommon.epoch.toi{1}              = [-5, 25];
configcommon.epoch.toi{2}              = [-5, 25];
configcommon.epoch.toi{3}              = [-5, 25];
configcommon.epoch.pad{1}              = 10;
configcommon.epoch.pad{2}              = 10;
configcommon.epoch.pad{3}              = 10;

configcommon.align.name                = {'SlowWave'};%{'SlowWave_R','SlowWave_R'};
configcommon.align.channel             = {'M1G'};      % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
configcommon.align.flip                = {'no'};
configcommon.align.abs                 = {'no',};
configcommon.align.method              = {'nearestmin'};      % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp'};%bp pour retirer les variations trop lentes de signal qui empÍchent de comparer la baseline et de chercher des maxima locaux
configcommon.align.freq                = {2};          % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no'};
configcommon.align.thresh.value        = [1, 1];
configcommon.align.thresh.method       = {'trial'};%'medianbl','both';
configcommon.align.toiplot             = {[-1,  1]};     % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.toiactive           = {[-0.5, 0.5]};  % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline         = {[-1, -0.5],[-1, -0.5],[-1, -0.5]};   % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.reref               = 'no';
configcommon.align.rerefmethod         = 'avg';
configcommon.align.refchannel          = 'all';
configcommon.align.notch               = 'no';

configcommon.LFP.name                  = {'SlowWave','SlowWave_EMGalign','SlowWavealign_EMG'};%,'SlowWave_L'}; %alwa
configcommon.LFP.emg                   = {'no','EMG1','EMG1'};%name of EMG channel associated with marker LFP.name. 'no' if no EMG associated to this seizure side 
configcommon.LFP.emgmarker             = {'no','no','SlowWave_EMG','SlowWave_EMG','SlowWave_EMG','SlowWave_EMG'};
configcommon.LFP.electrodetoplot       = {'M1G', 'M1G', 'M1G'};
configcommon.LFP.motorcortex           = {'M1G', 'M1G', 'M1G'};

configcommon.LFP.flip                  = true;
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 1;
configcommon.LFP.hpfiltord             = 4;
configcommon.LFP.hpfilttype            = 'but';
configcommon.LFP.resamplefs            = 1024;
configcommon.LFP.baseline              = 'no';
configcommon.LFP.baselinewindow        = {[-2, -1], [-2, -1], [-2, -1]};
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
configcommon.EMG.refchannel            = {'EMG2'};
configcommon.EMG.hpfilter              = 'yes';
configcommon.EMG.hpfreq                = 10;
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



%% Rodent 1
%NOTER QUE LES EMG SONT CHOISIS : REJETES SI PAS DE REPONSE EMG OU SI ARTEFACTS
config{1}                           = configcommon;
config{1}.prefix                    = 'DTX-EEGEMG-12-';
config{1}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-12');
config{1}.rawlabels.oldnames        = {'47','45','48','46','49'}; %for conversion deltamed to brainvision
config{1}.labels.newnames           = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
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


% %% Rodent 2
% config{2}                           = configcommon;
% config{2}.prefix                    = 'DTX-EEGEMG-14-';
% config{2}.rawdir                    = fullfile(rootpath_data,'DTX-EEGEMG-14');
% config{2}.rawlabels                 = {'57','55','58','56','59'};
% config{2}.labels.macro              = {'M1G','M1D','PtA','EMG1','EMG2'};
% config{2}.imagesavedir              = fullfile(imagesavedir,'DTX-EEGEMG-12');       % where to print images
% config{2}.muse.startend             = {'SlowWave','SlowWave'; 'SlowWave_EMG','SlowWave_EMG'; 'Crise_Start','Crise_End'};   % start and end Muse marker. For defining trials
% config{2}.align.name                = {'SlowWave','SlowWave_EMG'};%{'SlowWave_R','SlowWave_R'};
% config{2}.align.channel             = {'M1G','M1G'};      % pattern to identify channel on which to based peak detection % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.LFP.name                  = {'SlowWave','SlowWave_EMG'};%,'SlowWave_L'}; %alwa
% config{2}.LFP.emg                   = {'no','EMG1'};%name of EMG channel associated with marker LFP.name. 'no' if no EMG associated to this seizure side 
% config{2}.EMG.refchannel            = {'EMG2'};
% 
% %config{irat}.directorylist
% filelist = dir(config{2}.rawdir);
% i=0;
% for ifile = 1:length(filelist)
%     [~,~,file_extension] = fileparts(filelist(ifile).name);
%     if strncmp(file_extension,'.eeg',4)
%         i=i+1;
%         config{2}.directorylist{1}{i}          =  filelist(ifile).name(1:end-4);
%     end
% end
% clear filelist

end


