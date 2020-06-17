function writeSpykingCircusParameters(cfg)
% write .params and .prb file 
% .params file is based on the one in the template folder of EpiCode.
% Change this file to have new default values.
% write slurm files and create slurm dir to store the spyking circus
% outputs when launched on the cluster (work only for paths begining with
% '/network/lustre/iss01' if this script is launched on windows

[p, ~, ~]               = fileparts(mfilename('fullpath'));
[p, ~, ~]               = fileparts(p);
fname_params_default    = fullfile(p,'templates','SpykingCircus.params');

for ipart = 1 : size(cfg.directorylist,2)
    
    subjdir         = cfg.prefix(1:end-1);
    partdir         = ['p',num2str(ipart)];
    filename        = [cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1},'.params'];
    fname_params    = fullfile(cfg.datasavedir,subjdir,partdir,filename);
    nb_channels     = size(cfg.circus.channel,2);
    fname_prb       = ['Adtech_', num2str(nb_channels), 'chan.prb'];

    % read Spyking-Circus params file
    ini = IniConfig();
    ini.ReadFile(fname_params_default);
    
    % remove inline comments
    [sections, count_sections] = ini.GetSections();
    for sectioni = 1 : count_sections
        [keys, count_keys] = ini.GetKeys(sections{sectioni});
        for keysi = 1 : count_keys
            old = ini.GetValues(sections{sectioni}, keys{keysi});
            temp = split(old,{'#'});
            ini.SetValues(sections{sectioni}, keys{keysi}, temp{1});
        end
    end
    
    % adjust parameters
    h1 = ini.SetValues('data', {'file_format','stream_mode','mapping','suffix','overwrite','output_dir'}, {'neuralynx','None',fname_prb,'','False','SpykingCircus'});
    h2 = ini.SetValues('noedits', {'filter_done','artefacts_done','ground_done','median_done'}, {'False','False','False','False'});
    h3 = ini.SetValues('triggers', {'dead_file','dead_unit','ignore_times'}, {'SpykingCircus_artefacts_ms.dead','ms','True'});
    if any([h1; h2; h3] ~= 1), error('Something went wrong with adjusting parameters'); end
    
    status = ini.WriteFile(fname_params);
    if status == false
        error('Couldn''t write file');
    end
    ini.ToString()
    
    % write params file 
    writeProbeFile(nb_channels,fullfile(cfg.datasavedir,subjdir,partdir,fname_prb));
    
    %% write 2 slurm files : one for launching spyking circus, the other to only process template matching 
    
    %deal with issue that this script might be launched on pc wheread slurm script needs to be launched on unix
    if ispc
        temp = split(cfg.datasavedir,'\');
        datadir = '/network/lustre/iss01';
        for i = 1:size(temp,1)
            if ~isempty(temp{i}) && ~strcmp(temp{i}, 'lexport')
                temp{i} = strrep(temp{i}, 'iss01.', []);
                datadir = [datadir,'/',temp{i}];
            end
        end
    else
        datadir = cfg.datasavedir;
    end
    
    for islurm_file = [1 2] 
        if islurm_file == 1
            fname = fullfile(cfg.datasavedir,subjdir,partdir,'slurm_spyking-circus.sh');
        else
            fname = fullfile(cfg.datasavedir,subjdir,partdir,'slurm_spyking-circus_extracting.sh');
        end
        
        fid = fopen(fname,'w+');
        
        fprintf(fid,'#!/bin/bash\n');
        if islurm_file == 1
            fprintf(fid,'#SBATCH --job-name=spyking-circus\n');
        else
            fprintf(fid,'#SBATCH --job-name=sc_extracting\n');
        end
        fprintf(fid,'#SBATCH --partition=normal,bigmem\n');
        fprintf(fid,'#SBATCH --time=24:00:00\n');
        fprintf(fid,'#SBATCH --mem=120G\n');
        fprintf(fid,'#SBATCH --cpus-per-task=28\n');
        fprintf(fid,'#SBATCH --chdir=.\n');
        fprintf(fid,'#SBATCH --output=%s/slurm_output/output-%%j-%%x.txt\n',datadir);
        fprintf(fid,'#SBATCH --error=%s/slurm_output/error-%%j-%%x.txt\n\n',datadir);
        
        fprintf(fid,'module load spyking-circus/0.9.1\n');
        if islurm_file == 1
            fprintf(fid,'spyking-circus %s.ncs -c 28\n', filename(1:end-7));
        else
            fprintf(fid,'spyking-circus %s.ncs -m whitening,extracting,fitting -c 28\n', filename(1:end-7));
        end
        fprintf(fid,'spyking-circus %s.ncs -m converting -c 28\n', filename(1:end-7));
        fprintf(fid,'sleep 5\n');
        
        fclose(fid);
    end
    
    %create dir to store slurm text outputs and errors
    mkdir(fullfile(cfg.datasavedir,subjdir,partdir,'slurm_output'));
end

