function hspike_cluster_spikestats(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

hspike_setpaths;

%% General analyses
config                                      = hspike_setparams;
[MuseStruct{ipatient}, ~, ~]                = detectTemplate(config{ipatient});
[config{ipatient}, MuseStruct{ipatient}]    = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

% load raw data
SpikeRaw{ipatient}                          = readSpikeRaw_Phy(config{ipatient}, false);

% segment into trials/segments
config{ipatient}.spike.name                 = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6', 'window'};
SpikeTrials{ipatient}                       = readSpikeTrials(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, false);

% statistics
% config{ipatient}.spike.name                 = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6', 'window'};
% SpikeStats{ipatient}                        = spikeTrialStats(config{ipatient}, SpikeTrials{ipatient}, false);

% waveforms
% SpikeWaveforms{ipatient}                    = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, true);

% spike density, not for window
config{ipatient}.spike.psth.name     = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
SpikeDensity{ipatient}               = spikePSTH(config{ipatient}, SpikeTrials{ipatient}, true);
