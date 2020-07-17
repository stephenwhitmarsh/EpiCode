%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%


%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\  
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\      
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses

% load settings
config = hspike_setparams;


for ipatient = 2 : 4   
    
    % load settings
    config = hspike_setparams;
        
    [MuseStruct_orig]                       = readMuseMarkers(config{ipatient}, false);
    [MuseStruct_aligned]                    = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig, false);
    [clusterindx, LFP_cluster]              = clusterLFP(config{ipatient}, MuseStruct_aligned, true);
    
    
    templates                               = LFP_cluster{1}{1}.kmedoids{6};
    MuseStruct_template                     = MuseStruct_aligned;
    
    itemp = 1;
    config{ipatient}.template.name          = sprintf('template%d', itemp);
    [MuseStruct_template, indx, LFP_avg]    = detectTemplate(config{ipatient}, MuseStruct_template, templates{1}, true);
    
    config{ipatient}.align.name             = {'template1'};
    config{ipatient}.align.latency          = config{ipatient}.template.latency;
    config{ipatient}.align.write            = false;
    [MuseStruct_template_aligned]           = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_template, true);
    
    fname = fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'MuseStruct_alignedXcorr.mat']);
    save(fname, 'MuseStruct_template_aligned');
    
    [PSGtable]                              = PSG2table(config{ipatient}, MuseStruct_template_aligned, false);
    [t]                                     = plotHypnogram(config{ipatient}, MuseStruct_template_aligned);
    [marker, hypnogram]                     = hypnogramStats(config{ipatient}, MuseStruct_template_aligned, true);
    
end
