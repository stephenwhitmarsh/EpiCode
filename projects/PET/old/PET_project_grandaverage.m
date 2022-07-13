function PET_project_grandaverage

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\shared'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\external'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\templates'))
    addpath (genpath('\\lexport\iss02.charpier\analyses\vn_pet\EpiCode\projects/PET'))
    addpath \\lexport\iss02.charpier\analyses\vn_pet\MatlabImportExport_v6.0.0
    addpath \\lexport\iss02.charpier\analyses\vn_pet\fieldtrip
    
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

% classification des units

% plots des parametres inter pointes selon le groupe (hypo/normo
% metabolique)

