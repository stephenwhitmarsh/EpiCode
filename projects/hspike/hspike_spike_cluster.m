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
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries/ 
end
if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 
end

ft_defaults

% labels to update
labels = ["BAD__START__", "BAD__END__", "PHASE_1__START__", "PHASE_1__END__", "PHASE_2__START__", "PHASE_2__END__", "PHASE_3__START__", "PHASE_3__END__", "REM__START__", "REM__END__", "AWAKE__START__", "AWAKE__END__", "NO_SCORE__START__", "NO_SCORE__END__"];

for ipatient = 1 : 3
    
    config                      = hspike_setparams;
    MuseStruct{ipatient}        = readMuseMarkers(config{ipatient}, false);
    MuseStruct{ipatient}        = updateMarkers(config{ipatient}, MuseStruct{ipatient}, labels);
    [~, hypnogram{ipatient}, ~] = hypnogramMuseStats(config{ipatient});
    
    % trim files to only those within a hypnogram
    MuseStruct_trimmed{ipatient}  = MuseStruct{ipatient};
    config_trimmed{ipatient}      = config{ipatient};
    for ipart = 1 : 3
        sel     = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
        first   = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
        last    = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
        config_trimmed{ipatient}.directorylist{ipart}   = config{ipatient}.directorylist{ipart}(first:last);
        MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct{ipatient}{ipart}(first:last);
    end
    
    writeSpykingCircusDeadfiles(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, true);
    writeSpykingCircusFileList(config_trimmed{ipatient}, true);
    writeSpykingCircusParameters_new(config_trimmed{ipatient});
    
end


