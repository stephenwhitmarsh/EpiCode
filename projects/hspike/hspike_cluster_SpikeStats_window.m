function hspike_cluster_SpikeStats_window(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

hspike_setpaths;

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

labels = ["BAD__START__", "BAD__END__", ...
    "PHASE_1__START__", "PHASE_1__END__", ...
    "PHASE_2__START__", "PHASE_2__END__", ...
    "PHASE_3__START__", "PHASE_3__END__", ...
    "REM__START__", "REM__END__", ...
    "AWAKE__START__", "AWAKE__END__", ...
    "PRESLEEP__START__", "PRESLEEP__END__", ...
    "POSTSLEEP__START__", "POSTSLEEP__END__"];

%% General analyses

config                                                      = hspike_setparams;
MuseStruct{ipatient}                                        = readMuseMarkers(config{ipatient}, false);
MuseStruct{ipatient}                                        = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
[~, LFP_cluster{ipatient}]                                  = clusterLFP(config{ipatient}, MuseStruct{ipatient}, false);
[MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}]   = detectTemplate(config{ipatient}, MuseStruct{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
[config{ipatient}, MuseStruct{ipatient}]                    = alignTemplates(config{ipatient}, MuseStruct{ipatient}, LFP_cluster_detected{ipatient});
MuseStruct{ipatient}                                        = updateMarkers(config{ipatient}, MuseStruct{ipatient}, labels);
MuseStruct{ipatient}                                        = padHypnogram(MuseStruct{ipatient});
[config{ipatient}, MuseStruct{ipatient}]                    = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

% load raw data
SpikeRaw{ipatient}          = readSpikeRaw_Phy(config{ipatient}, false);

% segment into trials/segments
config{ipatient}.spike.name = {'window'};
SpikeTrials{ipatient}       = readSpikeTrials(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, true);

% statistics
config{ipatient}.spike.name = {'window'};
SpikeStats{ipatient}        = spikeTrialStats(config{ipatient}, SpikeTrials{ipatient}, true);