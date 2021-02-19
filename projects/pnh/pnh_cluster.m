h_function pnh_project_cluster(ipatient)

%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

%% Add path

restoredefaultpath

%% Stephen's paths

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/pnh/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\pnh
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
end
% 
% %% Valerio
% if isunix
%     addpath /network/lustre/iss01/charpier/analyses/valerio.frazzini/fieldtrip
%     addpath /network/lustre/iss01/charpier/analyses/valerio.frazzini/EpiCode/projects/pnh/
%     addpath /network/lustre/iss01/charpier/analyses/valerio.frazzini/EpiCode/shared/
%     addpath /network/lustre/iss01/charpier/analyses/valerio.frazzini/EpiCode/shared/utilities
%     addpath(genpath('/network/lustre/iss01/charpier/analyses/valerio.frazzini/scripts/releaseDec2015/'));
% end
% 
% if ispc
%     addpath \\lexport\iss01.charpier\analyses\valerio.frazzini\fieldtrip
%     addpath \\lexport\iss01.charpier\analyses\valerio.frazzini\EpiCode\projects\pnh
%     addpath \\lexport\iss01.charpier\analyses\valerio.frazzini\EpiCode\shared
%     addpath \\lexport\iss01.charpier\analyses\valerio.frazzini\EpiCode\shared\utilities
%     addpath \\lexport\iss01.charpier\analyses\valerio.frazzini\MatlabImportExport_v6.0.0
% end

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% load settings
config = pnh_setparams;

% read muse markers
% [MuseStruct{ipatient}] = readMuseMarkers(config{ipatient}, false);

% align markers
% [MuseStruct_aligned{ipatient}]          = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);

% % read LFP data (trial-by-trial)
% LFP{ipatient} = readLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
% 
% % make TFR
% TFR{ipatient} = TFRtrials(config{ipatient}, LFP{ipatient}, false);

% write data and parameters for spyking circus
% writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct_aligned{ipatient}, false);

config{ipatient}.circus.version                = 'fieldtrip';
config{ipatient}.circus.params.detection.peaks = 'positive';

writeSpykingCircusParameters(config{ipatient});
writeSpykingCircus(config{ipatient}, true);
