function hspike_cluster(ipatient)

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
config                                          = hspike_setparams;
[MuseStruct_orig{ipatient}]                     = readMuseMarkers(config{ipatient}, true);
[MuseStruct_aligned{ipatient}]                  = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, true);
[clusterindx{ipatient}, LFP_cluster{ipatient}]  = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, true);
[MuseStruct_template{ipatient}, ~,~, ~]         = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, true);

% add templates to config
for itemp = 1 : 6
    markername = sprintf("template%d", itemp);
    config{ipatient}.muse.startmarker.(markername)                                              = markername;
    config{ipatient}.muse.endmarker.(markername)                                                = markername;
    config{ipatient}.epoch.toi.(markername)                                                     = [-0.5  1];
    config{ipatient}.epoch.pad.(markername)                                                     = 0.5;
    config{ipatient}.LFP.baselinewindow.(markername)                                            = [-0.5  1];
    config{ipatient}.LFP.baselinewindow.(markername)                                            = [-0.5  1];
    config{ipatient}.LFP.name{itemp}                                                            = markername;
    config{ipatient}.hyp.markers{itemp}                                                         = markername;
end

[t{ipatient}]                                   = plotHypnogram(config{ipatient}, MuseStruct_template{ipatient});
[markers{ipatient}, hypnogram{ipatient}]        = hypnogramStats(config{ipatient}, MuseStruct_template{ipatient}, true);
[LFP{ipatient}]                                 = readLFP(config{ipatient}, MuseStruct_template{ipatient}, true);
[LFP_stage{ipatient}]                           = plotLFP_stages(config{ipatient}, LFP{ipatient}, true);
