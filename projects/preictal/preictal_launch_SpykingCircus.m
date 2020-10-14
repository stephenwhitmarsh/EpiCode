function preictal_launch_SpykingCircus(ipatient)

%find config script :
addpath(genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/preictal'));
config = preictal_setparams;

%launch spykling circus : 
for ipart = 1:size(config{ipatient}.directorylist,2)
    
    subjdir     = config{ipatient}.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    filename    = [config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs'];
    dirname     = fullfile(config{ipatient}.datasavedir,subjdir,partdir);
    
    load_sc         = 'module load spyking-circus/0.9.9;';
    open_dir        = sprintf('cd %s;', dirname);
    launch_sc       = sprintf('spyking-circus %s -c 28;', filename);
    convert_results = sprintf('spyking-circus %s -m converting;', filename);
    
    eval(['!',load_sc,open_dir,launch_sc,convert_results]);
    
end
