function hspike_cluster_template(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% preamble
feature('DefaultCharacterSet', 'CP1252')
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
hspike_setpaths;
labels = ["BAD__START__", "BAD__END__", ...
    "PHASE_1__START__", "PHASE_1__END__", ...
    "PHASE_2__START__", "PHASE_2__END__", ...
    "PHASE_3__START__", "PHASE_3__END__", ...
    "REM__START__", "REM__END__", ...
    "AWAKE__START__", "AWAKE__END__", ...
    "PRESLEEP__START__", "PRESLEEP__END__", ...
    "POSTSLEEP__START__", "POSTSLEEP__END__"];

%% set up config and musestruct, and read template
config                                      = hspike_setparams;
config{ipatient}                            = addparts(config{ipatient});
MuseStruct{ipatient}                        = readMuseMarkers(config{ipatient}, false);
MuseStruct{ipatient}                        = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
MuseStruct{ipatient}                        = padHypnogram(MuseStruct{ipatient});
[~, ~, LFP_cluster{ipatient}]               = clusterLFP(config{ipatient}, MuseStruct{ipatient}, false);
[MuseStruct{ipatient}, ~, ~]                = detectTemplate(config{ipatient}, MuseStruct{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
[config{ipatient}, MuseStruct{ipatient}]    = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

%% LFP averages of templates
config{ipatient}.LFP.postfix                = {'_all'};
config{ipatient}.LFP.name                   = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
LFPavg{ipatient}                            = readLFPavg(config{ipatient}, MuseStruct{ipatient}, true);

% %% LFP of sliding timewindow
% config{ipatient}.LFP.name                   = {'window'};
% config{ipatient}.LFP.postfix                = {'_all'};
% LFP{ipatient}                               = readLFP(config{ipatient}, MuseStruct{ipatient}, true);
% 
% %% FFT on sliding timewindow
% config{ipatient}.FFT.name                   = {'window'};
% config{ipatient}.FFT.postfix                = {'_all'};
% FFT{ipatient}                               = FFTtrials(config{ipatient}, true);
% 
%% read LFP of only first three parts
config{ipatient}.directorylist              = config{ipatient}.directorylist(1:3);
MuseStruct{ipatient}                        = MuseStruct{ipatient}(1:3);
config{ipatient}.LFP.postfix                = [];
config{ipatient}.LFP.name                   = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
LFP{ipatient}                               = readLFP(config{ipatient}, MuseStruct{ipatient}, true);

%% rereference
config{ipatient}.LFP.name                   = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
LFP{ipatient}                               = rerefLFP(config{ipatient}, MuseStruct{ipatient}(1:3), true);

