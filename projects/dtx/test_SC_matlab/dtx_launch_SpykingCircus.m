function dtx_launch_SpykingCircus(irat)

addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));

config = dtx_setparams_probe_spikes;

for ipart = 1:size(config{irat}.directorylist,2)
    
    subjdir     = config{irat}.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    filename    = [config{irat}.prefix,'p',num2str(ipart),'-multifile-',config{irat}.circus.channel{1},'.ncs'];
    dirname     = fullfile(config{irat}.datasavedir,subjdir,partdir);
    
    
    load_sc         = 'module load spyking-circus/0.9.9;';
    open_dir        = sprintf('cd %s;', dirname);
    launch_sc       = sprintf('spyking-circus %s -c 28;', filename);
    convert_results = sprintf('spyking-circus %s -m converting;', filename);
    
    eval(['!',load_sc,open_dir,launch_sc,convert_results]);
    
end
