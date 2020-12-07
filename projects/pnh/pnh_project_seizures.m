%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/altmany-export_fig-8b0ba13
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
    %     addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    %     addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses

config = pnh_setparams_seizures;

for ipatient = 1:3
    
    % read channels of seizures onsets
    seizures = readtable('\\lexport\iss01.charpier\analyses\vn_pnh\Seizure_description_rows.xlsx');
    sel      = seizures(seizures.patientnr == ipatient, :);
    config{ipatient}.directorylist{1} = unique(sel.path)';
    config{ipatient}.LFP.channel = unique(sel.chan1)';
    ind = strcmp(config{ipatient}.LFP.channel, 'diffuse');
    
    % add channels of seizures onsets to config 
    config{ipatient}.LFP.channel = config{ipatient}.LFP.channel(~ind);
    for ichan = 1 : size(config{ipatient}.LFP.channel, 2)
        config{ipatient}.LFP.channel{ichan} = ['_', config{ipatient}.LFP.channel{ichan}];
    end
    
    % read muse markers
    [MuseStruct{ipatient}] = readMuseMarkers(config{ipatient}, false);
    
    % read LFP data
    LFP{ipatient} = readLFP(config{ipatient}, MuseStruct{ipatient}, false);
    
    % add seizure onset channels to data
    if height(sel) < height(LFP{ipatient}{1}.seizure.trialinfo)
        LFP{ipatient}{1}.seizure.trialinfo.chan1(1:height(sel)) = string(sel.chan1);
        LFP{ipatient}{1}.seizure.trialinfo.chan2(1:height(sel)) = string(sel.chan2);
    else
        LFP{ipatient}{1}.seizure.trialinfo.chan1 = string(sel.chan1(1:height(LFP{ipatient}{1}.seizure.trialinfo)));
        LFP{ipatient}{1}.seizure.trialinfo.chan2 = string(sel.chan2(1:height(LFP{ipatient}{1}.seizure.trialinfo)));
    end
    
    % calculate TFR for all seizures
    TFR{ipatient} = TFR_seizures(config{ipatient}, LFP{ipatient}, false);
    
    % plot overview of seizures, per channel of onset
    plotTimeCourses_seizures(config{ipatient}, LFP{ipatient}, TFR{ipatient});
    
    % save data to excel
    writetable(LFP{ipatient}{1}.seizure.trialinfo, strcat(config{ipatient}.datasavedir, config{ipatient}.prefix, 'MUSEdata.xls'))

end





for ipatient = 1 : 3
end


    if ipatient == 1
        cfg = [];
        cfg.trials = 1 : size(LFP{1}.seizure.trial, 2);
        cfg.trials([23, 24, 27, 28]) = [];
        LFP{1}.seizure = ft_selectdata(cfg, LFP{1}.seizure);
    end
    
    
    % plot LFP timecourse examples for article
    
end
     

