function hspike_cluster(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/SPIKY_apr_2021'))
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\SPIKY_apr_2021'));
end
ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

labels = ["BAD__START__", "BAD__END__", "PHASE_1__START__", "PHASE_1__END__", "PHASE_2__START__", "PHASE_2__END__", "PHASE_3__START__", "PHASE_3__END__", "REM__START__", "REM__END__", "AWAKE__START__", "AWAKE__END__", "NO_SCORE__START__", "NO_SCORE__END__"];

%% General analyses, looping over patients
config = hspike_setparams;

MuseStruct{ipatient}                                        = readMuseMarkers(config{ipatient}, false);
MuseStruct{ipatient}                                        = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
[clusterindx{ipatient}, LFP_cluster{ipatient}]              = clusterLFP(config{ipatient}, MuseStruct{ipatient}, false);
[MuseStruct{ipatient}, ~,~, LFP_cluster_detected{ipatient}] = detectTemplate(config{ipatient}, MuseStruct{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
MuseStruct{ipatient}                                        = alignTemplates(config{ipatient}, MuseStruct{ipatient}, LFP_cluster_detected{ipatient});
MuseStruct{ipatient}                                        = updateMarkers(config{ipatient}, MuseStruct{ipatient}, labels);

% add (sliding) timewindow
[config{ipatient}, MuseStruct{ipatient}] = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

% template LFP
config{ipatient}.LFP.name   = {'window'};
LFP{ipatient}               = readLFP(config{ipatient}, MuseStruct{ipatient}, true);

% calculate FFT on sliding timewindow
config{ipatient}.FFT.name  = {'window'};
FFT{ipatient}              = FFTtrials(config{ipatient}, true);


