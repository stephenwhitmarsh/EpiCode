
%% Analyse manips WOD
%Paul Baudin

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\wod'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML
    CEDS64LoadLib('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\CEDMATLAB\CEDS64ML');

elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip;
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/wod'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
end
ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx

config = wod_setparams([]);
cfg=config{1};

for irat = 1:length(config)
    ipart = 1;
    
    %% Get LFP data
    
% 	For readCEDMarkers and CED2mat : Can only be done on windows because 
%   of CEDS64 library. For Linux, set force argument to false, it will 
%   load the results.


%si besoin : ajouter ipart aux noms d'images et de data enregistrées

%read all events
%directorylist CED
%datasavedir pour data coupées et event
    [CEDStruct]               = readCEDmarkers(config{irat}, true);
    dat_LFP                   = readCEDcontinuous(config{irat},CEDStruct,true,true);
    
   %remove artefacts voir mon script 
   %test plot lfp pour voir si pas d'erreur
   imarker = 1;
   dtx_plot_overdraw_allchannels(config{irat},dat_LFP, ipart, imarker,config{irat}.epoch.toi{imarker}, true);
   dtx_plot_morpho(config{irat},dat_LFP,ipart,imarker,config{irat}.LFP.electrodetoplot{imarker}, [-0.02 0.02],false, true);

    %% Plot
    
