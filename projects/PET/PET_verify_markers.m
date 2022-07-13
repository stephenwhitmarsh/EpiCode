function PET_verify_markers

% Create a table with the wrong marker timings in
% \\l2xport\iss02.charpier\analyse\vn_pet\data\wrong_markers

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

% read configuration for all patients
config = PET_setparams;

for patient = [2]% => patients 1 2 3 et 4, [1 4] => 1 et 4, [1 4 7] => patients 1, 4 et 7
    
    % read muse markers
    % MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
    MuseStruct = readMuseMarkers(config{patient}, true);
    
    %verify markers
    wrong_markers = verifymarkers(config{patient}, MuseStruct);
end