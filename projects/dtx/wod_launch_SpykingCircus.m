function wod_launch_SpykingCircus(slurm_task_id)


addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607


ft_defaults
config = wod_setparams;

for irat = slurm_task_id
   for ipart = 1:size(config{irat}.directorylist,2)
    
    subjdir     = config{irat}.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    filename    = [config{irat}.prefix,'p',num2str(ipart),'-multifile-',config{irat}.circus.channel{1},'.ncs'];
    dirname     = fullfile(config{irat}.datasavedir,subjdir,partdir);
    
    eval(sprintf('!cd %s', dirname));
    eval(sprintf('!spyking-circus %s -c 28', filename));
    eval(sprintf('!spyking-circus %s -m converting', filename));
  
   end
    
    
end