function wod_project_antoine(irat)

% Ce script projet sert à calculer les données de chaque rat
% 'irat' est en input car il prendra la valeur du array slurm sur le cluster
% Rassembler les données de tous les rats pour faire des moyennes ou stats
% sur tous les rats : dans un autre script (si possible organisé comme
% celui-ci)

% a faire pour réorganiser les fonctions WOD : 
% - mettre en input : cfg, MuseStruct, LFP, et toutes les structures qui sont calculées 
% dans une autre fonction
% - retirer les addpath, les config = wod_setparams, les boucles irat
% - remplacer tous les config{irat} par cfg
% - si possible, sauvegarder les données en fin de script, et créer la 
% possibilité de les charger avec l'argument 'force'

%% set parameters
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparams;

%% analysis

%find concatenated LFP (see wod_concatenateLFP.m)
[~,dir_name]                       = fileparts(config{irat}.rawdir);
config{irat}.rawdir                = fullfile(config{irat}.concatdata_path);
config{irat}.directorylist{ipart}  = {dir_name};

%read Muse markers 
MuseStruct = readMuseMarkers(config{irat}, true);

%read LFP, append electrodes, and cut into trials according to Muse Markers
LFP = readLFP(config{irat}, MuseStruct, true);
LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
%end
%vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
    error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
end