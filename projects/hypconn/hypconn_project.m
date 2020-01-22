%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis script for Hypconn project %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add correct paths

if isunix
    restoredefaultpath % clean path to avoid overloading functions
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hypconn/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    ft_defaults % add correct FieldTrip directories if not done yet.   
end

if ispc
    restoredefaultpath % clean path to avoid overloading functions
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hypconn
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    ft_defaults % add correct FieldTrip directories if not done yet.       
end

%% General analyses

% load settings for all patients
config = hypconn_setparams;

% loop over patients in settings file
for ipatient = 1:size(config,2)
      
    % read muse markers
    MuseStruct = readMuseMarkers(config{ipatient}, false);
    
    % read hypnogram as table, and write to file in data directory
    PSGtable   = PSG2table(config{ipatient}, MuseStruct, true);
    
end
