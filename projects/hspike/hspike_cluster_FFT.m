function hspike_cluster_FFT(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

hspike_setpaths;

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses

config                          = hspike_setparams;
config{ipatient}.FFT.name       = {'window'};
config{ipatient}.FFT.postfix    = {'_noWelch'};
FFT{ipatient}                   = FFTtrials(config{ipatient}, true);