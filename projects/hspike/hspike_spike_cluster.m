function hspike_spike_cluster(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%


%% Add path
restoredefaultpath

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries/ 
end
if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% get settings
config = hspike_setparams;

% get the artefacts
[MuseStruct_orig{ipatient}] = readMuseMarkers(config{ipatient}, false);

[marker{ipatient}, hypnogram{ipatient}] = hypnogramStats(config{ipatient}, MuseStruct_orig{ipatient}, false);

% trim files to only those within a hypnogram
MuseStruct_trimmed = MuseStruct_orig;
for ipart = 1 : 3
    sel = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
    first = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
    last = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
    config{ipatient}.directorylist{ipart} = config{ipatient}.directorylist{ipart}(first:last);
    MuseStruct_trimmed{ipatient}{ipart} = MuseStruct_trimmed{ipatient}{ipart}(first:last);
    if size(config{ipatient}.directorylist{ipart}, 2) > 7
        config{ipatient}.directorylist{ipart} = config{ipatient}.directorylist{ipart}(end-6:end);
        MuseStruct_trimmed{ipatient}{ipart} = MuseStruct_trimmed{ipatient}{ipart}(end-6:end);
    end
end

% write data concatinated for SC, artefacts, and output sampleinfo per file
writeSpykingCircus(config{ipatient}, MuseStruct_trimmed{ipatient}, true, true);

