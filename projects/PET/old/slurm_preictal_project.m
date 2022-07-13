function slurm_preictal_project(ielec)

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))
    %     addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_Dec_2019'))
    
    addpath \\lexport\iss01.charpier\analyses\vn_preictal\scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
    %     addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_Dec_2019'))
    addpath /network/lustre/iss01/charpier/analyses/vn_preictal/scripts/fieldtrip-20200607
end

ft_defaults

config = preictal_setparams;

%% prepare data for Spyking Circus
% write a new joblist for the cluster
% preictal_spikes_slurm_joblist

% read muse markers
MuseStruct = readMuseMarkers(config{ielec}, false);

% remove post ictal from the whole analysis,
% according to config (some 'patients' will have
% shorter postictal kept because of noise, see setparams)
MuseStruct = addMuseBAD(config{ielec}, MuseStruct);

% add (sliding) timewindow add window marker
[config{ielec}, MuseStruct] = addSlidingWindows(config{ielec}, MuseStruct);

% % template LFP
% config{ielec}.LFP.name   = {'window'};
% LFP               = readLFP(config{ielec}, MuseStruct, true);
% 
% % calculate FFT on sliding timewindow
% config{ielec}.FFT.name  = {'window'};
% FFT              = FFTtrials(config{ielec}, true);   % take a look i guess

%     % write artefacts to file
%     writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresNotRemoved');
%
%     %%
%     % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresNotRemoved
%     %%

%     % add arteafct marker from seizure start to seizure end
%     cfgtemp                       = [];
%     cfgtemp.bad.markerStart       = 'CriseStart';
%     cfgtemp.bad.markerEnd         = 'CriseEnd';
%     MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);
%
%     writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresRemoved');

%     %%
%     % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresRemoved
%     %%

%     % write parameters file for spyking circus
%     writeSpykingCircusParameters(config{ielec});
%
%     % write file list for spyking circus ---> LINUX
%     writeSpykingCircusFileList(config{ielec}, true);
%
%     %%
%     % Now do your spike sorting
%     %%
% 
% %% perform the analysis after spike sorting

%read spike data
SpikeRaw = readSpikeRaw(config{ielec}, true);

%read spike waveforms
SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, true);

% epoch data into windows
SpikeTrials                    = readSpikeTrials(config{ielec}, MuseStruct, SpikeRaw, true);

% calculate statistics per window
SpikeStats                    = spikeTrialStats(config{ielec}, SpikeTrials, false);

