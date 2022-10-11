function hspike_setpaths


restoredefaultpath
if isunix
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis    
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip    
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath /network/lustre/iss02/charpier/analyses/stephen.whitmarsh/scripts/BrainNetViewer_20191031
    addpath(genpath('/network/lustre/iss02/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));    
    addpath(genpath('/network/lustre/iss02/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss02/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'));  
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer    
end

if ispc
    addpath Z:\analyses\stephen.whitmarsh\git\fieldtrip
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\shared
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\external\subaxis    
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip    
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath Z:\analyses\stephen.whitmarsh\scripts\BrainNetViewer_20191031
    addpath(genpath('Z:\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('Z:\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          Z:\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
    addpath Z:\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer    
end

ft_defaults
