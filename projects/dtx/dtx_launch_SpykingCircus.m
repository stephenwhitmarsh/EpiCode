function dtx_launch_SpykingCircus(slurm_task_id)


addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));

ft_defaults
config = wod_setparams;

for irat = slurm_task_id
   for ipart = 1:size(config{irat}.directorylist,2)
    
    subjdir     = config{irat}.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    filename    = [config{irat}.prefix,'p',num2str(ipart),'-multifile-',config{irat}.circus.channel{1},'.ncs'];
    dirname     = fullfile(config{irat}.datasavedir,subjdir,partdir);
    
    eval('!module load spyking-circus/0.9.9');
    eval(sprintf('!cd %s', dirname));
    eval(sprintf('!spyking-circus %s -c 28', filename));
    eval(sprintf('!spyking-circus %s -m converting', filename));
  
   end
    
    
end