
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

function [config] = dtx_intra_setparams(config)

disp('setting parameters for dtx intracellular data');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-INTRA/';
    rootpath_data       = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-INTRA';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-INTRA';
    rootpath_data       = '\\lexport\iss01.charpier\analyses\lgi1\DTX-INTRA';
else
    error('Platform not supported')
end

%% recording test, to set the CED analysis strategy
                %VOIR SI BESOIN D'UTILISER lower(str) POUR CONSISTANCE NOMS 
%Peut être adapté la stratégie de récupération des markers si le noms ne
%sont pas consistants sur les différentes manips 

%subject infos
config{1}.datasavedir               = fullfile(rootpath_analysis,'data');
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'Crises-spont-intra');
config{1}.prefix                    = 'Intra_allseizures-';
config{1}.rawdir                    = fullfile(rootpath_data, 'Crises-spont-intra');
config{1}.directorylist{1}          = {'MergeCrises'}; %without the extension

%LFP analysis
% one field per LFP analysis
config{1}.LFP.name               = {'ECoG', 'Intra'};

% preprocessing before defining trials
config{1}.LFP.lpfilter           = {'no','no'};
config{1}.LFP.lpfreq             = {40,40};
config{1}.LFP.lpfilttype         = {'fir', 'fir'};
config{1}.LFP.doresample         = {false,false};
config{1}.LFP.resamplefs         = {300, 300};
config{1}.LFP.demean             = {'yes','yes'};
config{1}.LFP.baselinewindow     = {[-1 -0.5],[-1 -0.5]};%(before making trl)

%not used :
config{1}.LFP.hpfilter           = {'no','no'};
config{1}.LFP.hpfreq             = {1,1};
config{1}.LFP.hpfiltord          = {[],[]};%leave empty : default value
config{1}.LFP.hpfilttype         = {[],[]};%leave empty : default valueco
config{1}.LFP.bsfilter           = {'no','no'};
config{1}.LFP.bsfreq             = {[1 30], []};
config{1}.LFP.reref              = {'no','no'};
config{1}.LFP.rerefmethod        = {'avg',[]};
config{1}.LFP.refchannel         = {'all',[]};

%defining trials
config{1}.LFP.channel{1}            = {'v_ECoG-M1G'};
config{1}.LFP.channel{2}            = {'v_Vm'};
config{1}.LFP.dorename           = {'no','no'}; %if trials have to be merged from several channels to one unique channels, write the name of the output channel (set 'no' to ignore)
config{1}.startend               = {'SlowWave','SlowWave';'SlowWave','SlowWave'};%(if dorename, one per channel in the imarker line, ie : 3:4, 5:6 etc.)
config{1}.epoch.toi              = {[-5, 5],[-5, 5]};
config{1}.epoch.pad              = {0,0};

% spike analysis
% Si spike analysis : utiliser CEDstruct pour writeSC et pour définir trl.
% Juste associer un fichier différent
