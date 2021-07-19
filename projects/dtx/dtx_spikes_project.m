function dtx_spikes_project(irat)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    ft_defaults
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/development'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    ft_defaults
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/Nlx2Mat_release-v7_Dec2015/binaries
end

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config = dtx_spikes_setparams;

%% analysis of spyking circus output

%read markers and LFP
MuseStruct          = readMuseMarkers(config{irat}, false);
MuseStruct          = dtx_remove_wrong_seizure(config{irat},MuseStruct,false);
MuseStruct          = alignMuseMarkersPeaks(config{irat},MuseStruct, false); %align to the begin of the bigger channel
MuseStruct_concat   = concatenateMuseMarkers(config{irat}, MuseStruct, false);
seizure_infos       = dtx_stats_seizure_timings(config{irat}, MuseStruct_concat,1);
LFP = dtx_readLFP(config{irat}, MuseStruct, false);
[config{irat},LFP] = dtx_correctDTX2name(config{irat},LFP); %correct an error in channel name during acquisition, for rat DTX2

%read spike data
SpikeRaw  = readSpikeRaw_Phy(config{irat},false);

% segment into trials based on Muse markers
SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(config{irat}, MuseStruct, SpikeRaw, false);
SpikeTrials_timelocked  = removeArtefactedTrials(config{irat}, SpikeTrials_timelocked);

%align interictal to the end
for ipart = 1:size(SpikeTrials_timelocked)
    for markername = "Interictal"
        if ~isfield(SpikeTrials_timelocked{ipart}, markername)
            continue
        end
        for itrial = 1:size(SpikeTrials_timelocked{ipart}.(markername).trialinfo,1)
            endtime   = SpikeTrials_timelocked{ipart}.(markername).trialtime(itrial,2);
            SpikeTrials_timelocked{ipart}.(markername).trialinfo.offset(itrial)     = SpikeTrials_timelocked{ipart}.(markername).trialinfo.offset(itrial) - endtime * SpikeTrials_timelocked{ipart}.(markername).hdr.Fs;
            SpikeTrials_timelocked{ipart}.(markername).trialtime(itrial,:)          = SpikeTrials_timelocked{ipart}.(markername).trialtime(itrial,:)      - endtime;
            for i_unit = 1:size(SpikeTrials_timelocked{ipart}.(markername).time,2)
                trial_idx = SpikeTrials_timelocked{ipart}.(markername).trial{i_unit} == itrial;
                SpikeTrials_timelocked{ipart}.(markername).time{i_unit}(trial_idx)      = SpikeTrials_timelocked{ipart}.(markername).time{i_unit}(trial_idx)      - endtime;
                SpikeTrials_timelocked{ipart}.(markername).sample{i_unit}(trial_idx)    = SpikeTrials_timelocked{ipart}.(markername).sample{i_unit}(trial_idx)    - endtime * SpikeTrials_timelocked{ipart}.(markername).hdr.Fs;
                SpikeTrials_timelocked{ipart}.(markername).timestamp{i_unit}(trial_idx) = SpikeTrials_timelocked{ipart}.(markername).timestamp{i_unit}(trial_idx) - endtime * SpikeTrials_timelocked{ipart}.(markername).hdr.Fs * SpikeTrials_timelocked{ipart}.(markername).hdr.TimeStampPerSample;
            end
        end
        config{irat}.stats.bl.Interictal            = [min(min(SpikeTrials_timelocked{ipart}.(markername).trialtime)), -100];
        config{irat}.LFP.baselinewindow.Interictal  = [min(min(SpikeTrials_timelocked{ipart}.(markername).trialtime)), -100];
        config{irat}.epoch.toi.Interictal           = [min(min(SpikeTrials_timelocked{ipart}.(markername).trialtime)), 0];
        config{irat}.spike.toi.Interictal           = [min(min(SpikeTrials_timelocked{ipart}.(markername).trialtime)), 0];
    end
end

%compute TFR
TFR = dtx_TFRtrials(config{irat}, LFP, false);

% segment into equal periods
% before windowing, set seizures as artefacts so they are removed too
MuseStruct_withoutSeizures  = addMuseBAD(config{irat}, MuseStruct);
SpikeTrials_windowed        = readSpikeTrials_windowed(config{irat}, MuseStruct_withoutSeizures, SpikeRaw, false);
SpikeTrials_windowed        = removeArtefactedTrials(config{irat}, SpikeTrials_windowed);

%compute stats
SpikeDensity_timelocked = spikeTrialDensity(config{irat}, SpikeTrials_timelocked, true);
config{irat}.postfix = 'windowed';
SpikeStats_windowed     = spikeTrialStats(config{irat}, SpikeTrials_windowed, true);%, 'windowed');

%read spike waveforms
SpikeWaveforms              = readSpikeWaveforms(config{irat}, SpikeTrials_windowed, false);
SpikeWaveforms_stats        = spikeWaveformStats(config{irat}, SpikeWaveforms, false);

plotOverviewDTX(config{irat}, SpikeTrials_timelocked, SpikeTrials_windowed,...
    SpikeStats_windowed, SpikeDensity_timelocked, LFP, TFR, SpikeWaveforms,...
    SpikeWaveforms_stats, MuseStruct_concat, seizure_infos);

summarized_neurons_table(config{irat}, 1, true, SpikeStats_windowed, SpikeDensity_timelocked, SpikeWaveforms_stats);