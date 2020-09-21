function preictal_launch_SpykingCircus(slurm_task_id)

addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/preictal'))

config = preictal_setparams;

for ipatient = slurm_task_id
   for ipart = 1:size(config{ipatient}.directorylist,2)
    
    subjdir     = config{ipatient}.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    filename    = [config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs'];
    dirname     = fullfile(config{ipatient}.datasavedir,subjdir,partdir);
    
    eval('!module load spyking-circus/0.9.9');
    eval(sprintf('!cd %s', dirname));
    eval(sprintf('!spyking-circus %s -c 28', filename));
    eval(sprintf('!spyking-circus %s -m converting', filename));
  
   end
    
    
end