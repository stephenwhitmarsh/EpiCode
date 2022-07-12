function PET_project

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\shared'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\external'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\templates'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\projects/PET'))
%     addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\SPIKY_apr_2021'))
    addpath \\lexport\iss02.charpier\analyses\vn_pet\MatlabImportExport_v6.0.0
    addpath \\lexport\iss02.charpier\analyses\vn_pet\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/external'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/projects/PET'))
%     addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/SPIKY_apr_2021'))
    addpath /network/lustre/iss02/charpier/analyses/vn_pet/fieldtrip
end
ft_defaults

% read configuration for all patients
config = PET_setparams;

% write a new joblist for the cluster
% PET_slurm_joblist

% select patient number
ipatient = 2;

% read muse markers
% MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
MuseStruct = readMuseMarkers(config{ipatient}, true); 

% template LFP
LFP = readLFP(config{ipatient}, MuseStruct, true);
 
% EXAMPLE: plot LFP
cfg = [];
LFPavg = ft_timelockanalysis(cfg, LFP{1}.spike_mAmT2);
figure;
plot(LFPavg.time, LFPavg.avg)

% write artefacts to file
writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct, true);

% write parameters file for spyking circus
writeSpykingCircusParameters(config{ipatient});

% write file list for spyking circus ---> LINUX (just your Virtual Machine)
writeSpykingCircusFileList(config{ipatient}, true);

% Now do your spike sorting: spyking-circus SpykingCircus.params

% Now run Phy and evaluate units (good/bad/noise)

%read spike data
SpikeRaw = readSpikeRaw_Phy(config{ipatient}, true);

%read spike waveforms
SpikeWaveforms = readSpikeWaveforms(config{ipatient}, SpikeRaw, true);

% epoch data into windows
SpikeTrials  = readSpikeTrials(config{ipatient}, MuseStruct, SpikeRaw, true);



%% and so forth