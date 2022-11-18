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
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis    
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip    
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\scripts\BrainNetViewer_20191031
    addpath(genpath('\\l2export\iss02.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\l2export\iss02.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\l2export\iss02.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
    addpath \\l2export\iss02.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer    
end

ft_defaults
