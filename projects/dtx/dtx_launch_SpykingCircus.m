function dtx_launch_SpykingCircus(irat, analysis)

% irat (mandatory) : number of the rat in the setparam script
% analysis (optional, default = 'all') : 
% - 'all' to perform the spike sorting + converting
% - 'extracting' to extract pre-computed templates and only do template
% matching
% - 'converting' to only convert spyking circus output to phy format 

arguments
    irat 
    analysis = ft_getopt([],'analysis','all');
end

addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));

config = dtx_setparams_probe_spikes;

switch analysis
    
    case 'all'
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

    case 'extracting'
        for ipart = 1:size(config{irat}.directorylist,2)
            
            subjdir     = config{irat}.prefix(1:end-1);
            partdir     = ['p',num2str(ipart)];
            filename    = [config{irat}.prefix,'p',num2str(ipart),'-multifile-',config{irat}.circus.channel{1},'.ncs'];
            dirname     = fullfile(config{irat}.datasavedir,subjdir,partdir);
            
            load_sc         = 'module load spyking-circus/0.9.9;';
            open_dir        = sprintf('cd %s;', dirname);
            launch_sc       = sprintf('spyking-circus %s -m whitening,extracting,fitting -c 28;', filename);
            convert_results = sprintf('spyking-circus %s -m converting;', filename);
            
            eval(['!',load_sc,open_dir,launch_sc,convert_results]);
            
        end
        
    case 'converting'
        for ipart = 1:size(config{irat}.directorylist,2)
            
            subjdir     = config{irat}.prefix(1:end-1);
            partdir     = ['p',num2str(ipart)];
            filename    = [config{irat}.prefix,'p',num2str(ipart),'-multifile-',config{irat}.circus.channel{1},'.ncs'];
            dirname     = fullfile(config{irat}.datasavedir,subjdir,partdir);
            
            load_sc         = 'module load spyking-circus/0.9.9;';
            open_dir        = sprintf('cd %s;', dirname);
            convert_results = sprintf('spyking-circus %s -m converting;', filename);
            
            eval(['!',load_sc,open_dir,convert_results]);
            
        end
        
    otherwise
        error('%s is not a possible option for analysis', analysis);
end