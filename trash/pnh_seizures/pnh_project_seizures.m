%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

%% Initialization

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
    %     addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    %     addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

% set up paths for FieldTrip
ft_defaults

% To fix bug for weird characters in reading neurlynx
feature('DefaultCharacterSet', 'CP1252') 

% go into debug mode when there is an error
dbstop if error

%% General analyses

% get parameters 
config = pnh_setparams_seizures;

% read channels of seizures onsets
seizures = readtable('\\lexport\iss01.charpier\analyses\vn_pnh\Seizure_description_rows_new.xlsx');


    
for ipatient = 2:3
    
    % select seizure from a single patient
    sel = seizures(seizures.patientnr == ipatient, :);

    % set up config for EpiCode functions
    config{ipatient}.directorylist{1}   = unique(sel.path)';
    config{ipatient}.LFP.channel        = [unique(sel.chan1)' unique(sel.control)'];
  
    % remove channels with diffuse / unclear onset from table
    ind = contains(config{ipatient}.LFP.channel, {'diffuse', '?', 'none'}) | contains(config{ipatient}.LFP.channel, {'diffuse', '?', 'none'});
    config{ipatient}.LFP.channel = config{ipatient}.LFP.channel(~ind);
    
    % set up config for EpiCode functions - directories to search or
    % seizure markers
    config{ipatient}.directorylist{1}   = unique(sel.path)';
  
    % add channels of seizures onsets to config  
    for ichan = 1 : size(config{ipatient}.LFP.channel, 2)
        config{ipatient}.LFP.channel{ichan} = ['_', config{ipatient}.LFP.channel{ichan}];
        b = config{ipatient}.LFP.channel{ichan}(1:end-1);
        n = num2str(str2double(config{ipatient}.LFP.channel{ichan}(end)) + 1);
        config{ipatient}.LFP.refchannel{ichan} = [b, n]; 
    end
    
    % read muse markers
    [MuseStruct{ipatient}] = readMuseMarkers(config{ipatient}, false);
    
    % read LFP data
    LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
    
    % add trialinfo
    LFP{ipatient}{1}.seizure.trialinfo.chan1    = string(sel.chan1);
    LFP{ipatient}{1}.seizure.trialinfo.control  = string(sel.control);
    LFP{ipatient}{1}.seizure.trialinfo.type     = string(sel.type);
    
    % remove seizure with diffuse / unclear onset from data
    ind = contains(LFP{ipatient}{1}.seizure.trialinfo.chan1, {'diffuse', '?', 'none'}) | contains(LFP{ipatient}{1}.seizure.trialinfo.control, {'diffuse', '?', 'none'});    
    cfg                         = [];
    cfg.trials                  = find(~ind);
    LFP{ipatient}{1}.seizure    = ft_selectdata(cfg, LFP{ipatient}{1}.seizure);
    
    % remove noisy seizures
    if ipatient == 1
        cfg = [];
        cfg.trials = 1:size(LFP{1}{1}.seizure.trial, 2);
        cfg.trials(18:19) = [];
        LFP{1}{1}.seizure = ft_selectdata(cfg, LFP{1}{1}.seizure);
    end
    if ipatient == 3
        cfg = [];
        cfg.trials = 1:size(LFP{3}{1}.seizure.trial, 2);
        cfg.trials(6) = [];
        LFP{3}{1}.seizure = ft_selectdata(cfg, LFP{3}{1}.seizure);
    end
    
    % calculate TFR for all seizures
    TFR{ipatient} = TFR_seizures(config{ipatient}, LFP{ipatient}, true);
    
    % plot overview of seizures, per channel of onset
    plotTimeCourses_seizures(config{ipatient}, LFP{ipatient}, TFR{ipatient});
    
end

for ipatient = 1:3

    
    % plot overview of seizures, per channel of onset
    plotTimeCourses_seizures(config{ipatient}, LFP{ipatient}, TFR{ipatient});
end

    % save data to excel
    writetable(LFP{ipatient}{1}.seizure.trialinfo, strcat(config{ipatient}.datasavedir, config{ipatient}.prefix, 'MUSEdata.xls'))



     

