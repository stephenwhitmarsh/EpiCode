function hspike_cluster_spikestats(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

hspike_setpaths;

%% General analyses

config                                                          = hspike_setparams;

% read latest markers after alignment
MuseStruct{ipatient}                                            = alignMuseMarkersXcorr(config{ipatient});

[clusterindx{ipatient}, ~, LFP_cluster{ipatient}]               = clusterLFP(config{ipatient}, MuseStruct{ipatient}, false);
[config{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}] = alignClusters(config{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6});

[MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}]       = detectTemplate(config{ipatient}, MuseStruct{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);

% combine templates into single IED detection markers
config{ipatient}.editmarkerfile.torename = {'template1', 'detection'; 'template2', 'detection'; 'template3', 'detection'; 'template4', 'detection'; 'template5', 'detection'; 'template6', 'detection'};
MuseStruct{ipatient} = editMuseMarkers(config{ipatient}, MuseStruct{ipatient});

% [config{ipatient}, MuseStruct{ipatient}]    = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

% template LFP
config{ipatient}.LFP.name   = {'detection'};
LFP{ipatient}               = readLFP(config{ipatient}, MuseStruct{ipatient}, false);


% load raw data
SpikeRaw{ipatient}                          = readSpikeRaw_Phy(config{ipatient}, false);

% segment into trials/segments
config{ipatient}.spike.name                 = {'detection'};
SpikeTrials{ipatient}                       = readSpikeTrials(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, true);
% SpikeTrials{ipatient}                       = readSpikeTrials(config{ipatient});

% statistics, not for templates
% config{ipatient}.spike.name                 = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6', 'window'};
% config{ipatient}.spike.name                 = {'window'};
% SpikeStats{ipatient}                        = spikeTrialStats(config{ipatient}, SpikeTrials{ipatient}, true);

% waveforms
% SpikeWaveforms{ipatient}                    = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, true);

% spike density, not for window
% config{ipatient}.spike.name           = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
SpikeDensity{ipatient}               = spikePSTH(config{ipatient}, SpikeTrials{ipatient}, true);
