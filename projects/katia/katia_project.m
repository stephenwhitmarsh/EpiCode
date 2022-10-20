function katia_project

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
config = katia_setparams;

for patient = 1
   

    
    % read muse markers
    % MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
    MuseStruct = readMuseMarkers(config{patient}, true);
    
    % add (sliding) timewindow
    [config{ipatient}, MuseStruct{ipatient}] = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});
    
    % write artefacts to file
    writeSpykingCircusDeadfiles(config{patient}, MuseStruct, true);
    
    % write parameters file for spyking circus
    writeSpykingCircusParameters(config{patient});
    
    % write file list for spyking circus
    writeSpykingCircusFileList(config{patient}, true);
    
    % write a new spyking-circus joblist for the cluster
    katia_slurm_joblist
end