function PET_project

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_PET\scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_PET\scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_PET\scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_PET\scripts\SPIKY_apr_2021'))
    addpath \\lexport\iss01.charpier\analyses\vn_PET\scripts\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_PET/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_PET/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_PET/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_PET/scripts/SPIKY_apr_2021'))
    addpath /network/lustre/iss01/charpier/analyses/vn_PET/scripts/fieldtrip
end
ft_defaults

% read configuration for all patients
config = PET_setparams;

% write a new joblist for the cluster
% PET_slurm_joblist

% select patient number
ipatient = 1;

% read muse markers
MuseStruct = readMuseMarkers(config{ipatient}, true);

% template LFP
LFP = readLFP(config{ipatient}, MuseStruct, true);

% write artefacts to file
writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct, true);

% write parameters file for spyking circus
writeSpykingCircusParameters(config{ipatient});

% write file list for spyking circus ---> LINUX (just your Virtual Machine)
writeSpykingCircusFileList(config{ielec}, true);

% Now do your spike sorting: spyking-circus SpykingCircus.params

% Now run Phy and evaluate units (good/bad/noise)

%read spike data
SpikeRaw = readSpikeRaw(config{ielec}, true);

%read spike waveforms
SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, true);

% epoch data into windows
SpikeTrials  = readSpikeTrials(config{ielec}, MuseStruct, SpikeRaw, true);

%% and so forth