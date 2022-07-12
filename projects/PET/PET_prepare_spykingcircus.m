function PET_prepare_spykingcircus

restoredefaultpath
if ispc
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\shared'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\external'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\templates'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_pet\EpiCode\projects/PET'))
    addpath \\l2export\iss02.charpier\analyses\vn_pet\MatlabImportExport_v6.0.0
    addpath \\l2export\iss02.charpier\analyses\vn_pet\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/external'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_pet/EpiCode/projects/PET'))
    addpath /network/lustre/iss02/charpier/analyses/vn_pet/fieldtrip
end
ft_defaults

% read configuration for all patients
config = PET_setparams;

for patient = 1 %[1:4] => patients 1 2 3 et 4, [1 4] => 1 et 4, [1 4 7] => patients 1, 4 et 7
   
    % write a new joblist for the cluster
    PET_slurm_joblist
    
    % read muse markers
    % MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
    MuseStruct = readMuseMarkers(config{patient}, true);
    
    % write artefacts to file
    writeSpykingCircusDeadfiles(config{patient}, MuseStruct, true);
    
    % write parameters file for spyking circus
    writeSpykingCircusParameters(config{patient});
    
    % write file list for spyking circus
    writeSpykingCircusFileList(config{patient}, true);
    
end