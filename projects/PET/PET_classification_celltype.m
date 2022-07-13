% PET_classification_celltype

% Try to separate units between putative PN and putative IN according to
% their waveform shape.

% To use with the maximum quantity of data as possible (do not use 
% with only few patient, wait to have almost or all patients included) 

% To use AFTER PET_project_per_patient because it use data computed in this
% script.

% Ne pas utiliser les résultats si les plots ne sont pas convaincants (ne 
% montrent pas une distribution bimodale)

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

config = PET_setparams;
ipart=1;

celltype = classification_celltype(config, true, true);