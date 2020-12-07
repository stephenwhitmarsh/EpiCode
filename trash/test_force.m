%test force
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development'));
addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
ft_defaults
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));

%load precomputed data for the test :
config = dtx_spikes_setparams;
SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(config{1}, [], [], false);
SpikeTrials_timelocked  = removeArtefactedTrials(config{1}, SpikeTrials_timelocked);

% do the test :
stats1 = spikeTrialDensity(config{1});
stats2 = spikeTrialDensity(config{1}, SpikeTrials_timelocked, false);
stats3 = spikeTrialDensity(config{1}, SpikeTrials_timelocked, true);

%for me, stats1, stats2 and stats3 are identical