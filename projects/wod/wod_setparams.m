
%% Setting parameters DTX project Paul Baudin

function [config] = wod_setparams(config)

disp('setting parameters');
muscale = 50; % multiunit scale

if ismac
    error('Platform not supported')
elseif isunix
%     rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/Analyses_Paul/';
%     rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
%     os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\WODtemp';
    rootpath_data       = '\\lexport\iss01.charpier\analyses\lgi1\DTX-INTRA\neurones-serum-phy';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data'); %removed of config{i} so we can more easily modify it
imagesavedir = rootpath_analysis;

%% recording test, to set the CED analysis strategy
                %VOIR SI BESOIN D'UTILISER lower(str) POUR CONSISTANCE NOMS 
%Peut être adapté la stratégie de récupération des markers si le noms ne
%sont pas consistants sur les différentes manips 
config{1}.prefix                    = 'WODtest-';

%For this project, Neuralync data are converted into Spike2 data for
%putting markers. 
%for getting Spike2 marker and convert it to be usable with Neuralynx data
config{1}.rawdir                 = fullfile(datasavedir,'data_converted');
config{1}.CEDrawdir              = fullfile(rootpath_data, 'test');
config{1}.directorylist{1}       =  {'40-01'}; %without the extension
config{1}.directorylist{2}       =  {'40-01'}; %
config{1}.labels.macro              = {'Vm','ECoG-M1G'};
config{1}.LFP.channel               = config{1}.labels.macro; %white space must be replaced by '_' in this field

%RE PASSER PAR NEURALYNX POUR ANALYSES ANTOINE

%config{1}.rawdir                    = 'dummy';
%config{1}.directorylist             = 'dummy';

config{1}.os                        = os;
config{1}.name                      = {'PA_moyen'};
config{1}.datasavedir               = datasavedir;

%trial between 2 files is ignored if eventindex == 1

% %1 align par mrk, 1 channel par mrk. 
% %boucle for pour remplir les paramètres communs
% config{1}.align.name                = {'all_WOD'};
% config{1}.align.channel             = 'dummy';%eval('MuseStruct{ipart}{idir}.{'E12LFP'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.flip                = {'no'};
% config{1}.align.abs                 = {'no'};
% config{1}.align.method              = {'min'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{1}.align.filter              = {'no'};
% config{1}.align.freq                = {5};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{1}.align.hilbert             = {'no'};
% config{1}.align.thresh.value        = [1, 1];
% config{1}.align.thresh.method       = {'trial','trial','trial'};%'medianbl','both';
% config{1}.align.toiplot{1}          = [-20,  20];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{1}        = [-10, 10];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{1}      = [-20, -10];
% % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% 


config{1}.LFP.electrodetopolot      = 'dummy'; %electrode by which the wod begins

config{1}.muse.startend             = {'PAmoyCt','PAmoyCt'};%; 'SlowWave','Crise_End'; 'Crise_End','SlowWave'};   % 'SlowWave','SlowWave'; for readLFP function : cut data ...s before SlowWave, and ...s after SlowWave
config{1}.muse.eventindex           = {[0 0]};%; [0 1]}; %index of the event related to the muse marker. ie : if is 1, take the next marker, else if it is 0, take the marker of the event

config{1}.LFP.name                  = {'MoyPA_ctrl'};
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 25000; %because sampling rate is 3200Hz
config{1}.LFP.baseline              = 'no';
config{1}.LFP.baselinewindow{1}     = [-2, -1];
config{1}.LFP.baselinewindow{2}     = [-2, -1];
config{1}.LFP.baselinewindow{3}     = [0, 1];
config{1}.LFP.slidestep             = 0.01;

config{1}.TFR.toi                      = [-15:0.01:35];
config{1}.TFR.baseline                 = 'no';%[-10 -5];
config{1}.TFR.baselinetype             = 'relchange';

% list of onset timing with respect to start-marker (s)
config{1}.epoch.toi{1}              = [-5, 25];
config{1}.epoch.toi{2}              = [-2, 1];
config{1}.epoch.toi{3}              = [1, -2];
config{1}.epoch.pad{1}              = 10;
config{1}.epoch.pad{2}              = 0.5;
config{1}.epoch.pad{3}              = 0.5;

config{1}.imagesavedir              = fullfile(imagesavedir, 'DTX5');       % where to print images

%config{1}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};

%config{1}.LFP.nr_chanCED            = [1];
%config{1}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};

end




