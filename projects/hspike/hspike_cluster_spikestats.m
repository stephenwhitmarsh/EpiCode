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

% template LFP (used for correlation with PSTH)
config{ipatient}.LFP.name                   = {'combined'};
LFP{ipatient}                               = readLFP(config{ipatient}, MuseStruct{ipatient}, false);

% load raw data
SpikeRaw{ipatient}                          = readSpikeRaw_Phy(config{ipatient}, false);

% Combine different templates (based on clusters) into single template
config{ipatient}.editmarkerfile.torename = {'template1', 'combined'; ...
                                            'template2', 'combined'; ...
                                            'template3', 'combined'; ...
                                            'template4', 'combined'; ...
                                            'template5', 'combined'; ...
                                            'template6', 'combined'};                                        
MuseStruct{ipatient} = editMuseMarkers(config{ipatient}, MuseStruct{ipatient});

% segment into trials/segments
config{ipatient}.spike.name                 = {'combined'};

SpikeTrials{ipatient}                       = readSpikeTrials(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, true);
% SpikeTrials{ipatient}                       = readSpikeTrials(config{ipatient});

% statistics, not for templates
% config{ipatient}.spike.name                 = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6', 'window'};
% config{ipatient}.spike.name                 = {'window'};
% SpikeStats{ipatient}                        = spikeTrialStats(config{ipatient}, SpikeTrials{ipatient}, true);

% waveforms
% SpikeWaveforms{ipatient}                    = readSpikeWaveforms(config{ipatient}, SpikeRaw{ipatient}, true);

% spike density, not for window
config{ipatient}.spike.name                 = {'combined'};
SpikeDensity{ipatient}                      = spikePSTH(config{ipatient}, SpikeTrials{ipatient}, true);
