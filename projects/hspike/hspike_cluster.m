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

[MuseStruct_orig{ipatient}]                                                                     = readMuseMarkers(config{ipatient}, false);
[MuseStruct_aligned{ipatient}]                                                                  = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
[~, LFP_cluster{ipatient}]                                                                      = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
[MuseStruct_template{ipatient}, ~,~, ~]                                                         = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, true);

switch ipatient
    case 1
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined1'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined2'};
    case 2
        markernames = {'combined1'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template4', 'combined1'; 'template6', 'combined1'};
    case 3
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined1'; 'template4', 'combined2'; 'template5', 'combined2'; 'template6', 'combined2'};
    case 4
        markernames = {'combined1', 'combined2', 'combined3'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined2'; 'template4', 'combined3'; 'template5', 'combined1';};
    case 5
        markernames = {'combined1', 'combined2', 'combined3'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined2'; 'template5', 'combined3'; 'template6', 'combined1';};
    case 6
        markernames = {'combined1'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template3', 'combined1'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined1'};
    case 7
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template2', 'combined1'; 'template3', 'combined1'; 'template4', 'combined1'; 'template5', 'combined2'; 'template6', 'combined2'};
end

MuseStruct_combined{ipatient} = editMuseMarkers(config{ipatient}, MuseStruct_template{ipatient});

% focus time period a bit more
config{ipatient}.epoch.toi.Hspike           = [-0.2  0.8];
config{ipatient}.epoch.pad.Hspike           = 0.5;
config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];

% from now on work on manual and combined templates
itemp = 2;
for markername = string(markernames)
    config{ipatient}.muse.startmarker.(markername)    = markername;
    config{ipatient}.muse.endmarker.(markername)      = markername;
    config{ipatient}.epoch.toi.(markername)           = [-0.2  0.8];
    config{ipatient}.epoch.pad.(markername)           = 0.5;
    config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
    config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
    config{ipatient}.LFP.name{itemp}                  = markername;
    config{ipatient}.hyp.markers{itemp}               = markername;
    itemp = itemp + 1;
end

config{ipatient}.LFP.write = true;

[t{ipatient}]                                                                                   = plotHypnogram(config{ipatient}, MuseStruct_combined{ipatient});
[marker{ipatient}, hypnogram{ipatient}]                                                         = hypnogramStats(config{ipatient}, MuseStruct_combined{ipatient}, true);
[LFP{ipatient}]                                                                                 = readLFP(config{ipatient}, MuseStruct_combined{ipatient}, true);
[~]                                                                                             = plotLFP_stages(config{ipatient}, LFP{ipatient}, marker{ipatient}, hypnogram{ipatient}, true);
