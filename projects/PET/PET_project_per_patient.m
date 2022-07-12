function PET_project_per_patient(ipatient)

restoredefaultpath
if ispc
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\shared'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\external'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\templates'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\projects/PET'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\development'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\SPIKY'))
    addpath \\l2export\iss02.charpier\analyses\vn_pet\MatlabImportExport_v6.0.0
    addpath \\l2export\iss02.charpier\analyses\vn_pet\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/external'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/development'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/projects/PET'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/SPIKY'))
    addpath /network/lustre/iss02/charpier/analyses/vn_pet/fieldtrip
end
ft_defaults

disp('Too short trials, and those containing BAD marker are removed for the analysis')

% read configuration for all patients
config = PET_setparams;
ipart=1;

%% To Do

% Read Ripplelab, quantifier occurence, duree, amplitude, frequence des
% oscillations. Xcorr par raport aux IED ? A quelle distance des IED ont
% lieu chaque HFO ?
    
%% Analysis

% read muse markers
% MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
MuseStruct = readMuseMarkers(config{ipatient}, true); 

% merge the different IED markers
cfgtemp                  = [];
cfgtemp.newmarker        = 'IEDall';
cfgtemp.markers_to_merge = config{ipatient}.muse.startmarker.interIED;
MuseStruct               = merge_MuseMarkers(cfgtemp, MuseStruct);

% define inter-IED periods to epoch data
[config{ipatient}, MuseStruct] = PET_add_interIED_epochs(config{ipatient},...
    MuseStruct, 'IEDall', config{ipatient}.epoch.toi.interIED, true);

%cut too long epochs into shorter epochs
new_trial_length = config{ipatient}.min_trial_length.interIED + sum(abs(config{ipatient}.epoch.toi.interIED)) + 1;
max_trial_length = new_trial_length * 2;
MuseStruct = PET_cut_long_epochs(config{ipatient}, MuseStruct, 'interIED',max_trial_length, new_trial_length);

% align muse markers 
MuseStruct = alignMuseMarkersXcorr(config{ipatient}, MuseStruct, true);

% template LFP
LFP = readLFP(config{ipatient}, MuseStruct, true);
LFP = remove_too_short_trials(config{ipatient}, LFP);
LFP = remove_artefacted_trials(config{ipatient}, LFP);

TFR = TFRtrials(config{ipatient}, true, LFP);
FFT = FFTtrials(config{ipatient}, true, LFP);

%read spike data
spikeraw = readSpikeRaw_Phy(config{ipatient}, true);

%read spike waveforms
spikewaveforms = readSpikeWaveforms(config{ipatient}, spikeraw, true);

%spike waveform stats
waveformstats = spikeWaveformStats(config{ipatient}, spikewaveforms, true);

% epoch data in trials according to muse markers
spiketrials  = readSpikeTrials(config{ipatient}, MuseStruct, spikeraw, true);
spiketrials  = remove_too_short_trials(config{ipatient}, spiketrials);
spiketrials  = remove_artefacted_trials(config{ipatient}, spiketrials);

% calculate statistics per trial
spikestats = spikeTrialStats(config{ipatient}, spiketrials, true);

% calculate psth and significant changes per unit
statsPSTH = spikePSTH(config{ipatient}, spiketrials, true, LFP);

% lecture Ripple Lab

%% plot unit overview
MuseStruct_concat = concatenateMuseMarkers(config{ipatient}, MuseStruct, true);


PET_plot_overview_units(config{ipatient}, spiketrials, spikestats, statsPSTH,...
    LFP, TFR, FFT, spikewaveforms, waveformstats, MuseStruct_concat);

%% output a table with the data computed for each unit

% table during the "interIED" periods
summarized_neurons_table(config{ipatient}, 'interIED', true, spikestats, waveformstats, spiketrials);

% table during the "IED" periods
PET_neurons_table_timelocked(config{ipatient}, ipart, statsPSTH, spiketrials, 'IED', true);

%% plot psth stats
statsPSTH2 = statsPSTH;
statsPSTH2{1}.psth = rmfield(statsPSTH{1}.psth, 'interIED');
plot_PSTH_stats(config{ipatient}, statsPSTH2);

%% plot unit waveform
config{ipatient}.plotspike.plotraw    = true;
config{ipatient}.plotspike.suffix     = '_raw';
config{ipatient}.plotspike.isi_lim    = [0 0.025];
config{ipatient}.plotspike.invert     = true;
config{ipatient}.plotspike.img_format = "png";
plot_spike_waveforms(config{ipatient}, "interIED", waveformstats, spikestats, spikewaveforms);

config{ipatient}.plotspike.plotraw    = false;
config{ipatient}.plotspike.suffix     = '_avg';
config{ipatient}.plotspike.img_format = ["png", "pdf"];
plot_spike_waveforms(config{ipatient}, "interIED", waveformstats, spikestats, spikewaveforms);

%% plot example of spike trial
clear toi
toi.interIED = [2 7]; %secondes
toi.IED      = "all";
for markername = ["IED", "interIED"]
    itrial = round(size(spiketrials{ipart}.(markername).trialinfo, 1)/2);
    plot_spike_trial_example(config{ipatient}, spiketrials, waveformstats, ipart, markername, itrial, toi.(markername));
end

%% plot raw and avg LFP with and without alignment
PET_plot_LFP_aligned(config{ipatient}, ipart, "IED");

