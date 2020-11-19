function hspike_cluster_SpykingCircus(ipatient)

%find config script :
addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode'));
config = hspike_setparams;

%launch spykling circus : 
for ipart = 1:size(config{ipatient}.directorylist,2)
    
    subjdir     = config{ipatient}.prefix(1:end-1);
    partdir     = ['p', num2str(ipart)];
    filename    = [config{ipatient}.prefix, 'p', num2str(ipart), '-multifile-', config{ipatient}.circus.channel{1}, '.ncs'];
    dirname     = fullfile(config{ipatient}.datasavedir, subjdir, partdir);    
    toexecute   = sprintf('! module load spyking-circus/0.9.9; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28', dirname, filename, filename);
    toexecute   = sprintf('! module load spyking-circus/0.9.9; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28', dirname, filename, filename);
   
    eval(toexecute);    
    
end
