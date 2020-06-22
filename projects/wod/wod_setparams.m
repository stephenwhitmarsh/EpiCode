
%% Setting parameters DTX project Paul Baudin

function [config] = wod_setparams(config)

disp('setting WOD parameters');

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

%subject infos
config{1}.datasavedir               = datasavedir;
config{1}.imagesavedir               = imagesavedir;
config{1}.prefix                    = 'WODtest-';
config{1}.rawdir                    = fullfile(rootpath_data, 'test');
config{1}.directorylist{1}          =  {'40-01'}; %without the extension
% config{1}.directorylist{2}          =  {'40-01'}; 
%ajouter partname ?

% infos about the peak to align. Already aligned in Spike2, but those infos
% are usefull for plotting scripts
config{1}.align.toibaseline{1}      = [-0.2 -0.1];
config{1}.align.toiactive{1}       = [-0.005, 0.005]; 

%LFP analysis
% one field per LFP analysis
config{1}.LFP.name{1}               = 'PAMoy';

% preprocessing before defining trials
config{1}.LFP.hpfilter{1}           = 'no';
config{1}.LFP.hpfreq{1}             = 1;
config{1}.LFP.hpfiltord{1}          = [];%leave empty : default value
config{1}.LFP.hpfilttype{1}         = [];%leave empty : default valueco
config{1}.LFP.bsfilter{1}           = 'no';
config{1}.LFP.bsfreq{1}             = [1 30];
config{1}.LFP.lpfilter{1}           = 'no';
config{1}.LFP.lpfreq{1}             = 30;
config{1}.LFP.lpfilttype{1}         = 'fir';
config{1}.LFP.reref{1}              = 'no';
config{1}.LFP.rerefmethod{1}        = 'avg';
config{1}.LFP.refchannel{1}         = 'all';
config{1}.LFP.doresample{1}         = false;
config{1}.LFP.resamplefs{1}         = 'no';
config{1}.LFP.baseline{1}           = 'no';
config{1}.LFP.baselinewindow{1}     = [0, 5];%(before making trl)

%defining trials
config{1}.LFP.channel{1}           = {'Vm'};
config{1}.LFP.flip{1}               = [false]; %one logical per channel
config{1}.LFP.dorename{1}           = 'no'; %if trials have to be merged from several channels to one unique channels, xrite the name of the output channel (set 'no' to ignore)
config{1}.startend(1, 1:2)          = {'PAmoyCt','PAmoyCt'};%(if dorename, one per channel in the imarker line, ie : 3:4, 5:6 etc.)
config{1}.eventindex{1}             = [0 0]; %index of the event related to the marker. ie : if is 1, take the next marker, else if it is 0, take the marker of the event
config{1}.epoch.toi{1}              = [-0.2, 0.2];
config{1}.epoch.pad{1}              = 0;

%plotting
config{1}.LFP.electrodetoplot{1}        = 'Vm'; %some of my script plot only one channel
config{1}.TFR.doTFR{1}              = true;
config{1}.TFR.toi{1}                = [-0.2, 0.2];
config{1}.TFR.baseline{1}           = 'no';%[-10 -5];
config{1}.TFR.baselinetype{1}       = 'relchange';

% spike analysis
% Si spike analysis : utiliser CEDstruct pour writeSC et pour définir trl.
% Juste associer un fichier différent


end




