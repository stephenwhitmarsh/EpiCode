function writeSpykingCircusParameters(cfg)

fname_params_default        = 'SpykingCircusDefaultSettings.params'; % in path, i.e. in /shared scripts

for ipart = 1 : size(cfg.directorylist,2)
    
    subjdir         = cfg.prefix(1:end-1);
    partdir         = ['p',num2str(ipart)];
    filename        = [cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1},'.params'];
    fname_params    = fullfile(cfg.datasavedir,subjdir,partdir,filename);
    fname_prb       = [cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1},'.prb'];
    
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
    h1 = ini.SetValues('data', {'file_format','stream_mode','mapping','suffix','overwrite','output_dir','ncs_pattern'}, {'neuralynx','None',fname_prb,'','False','SpykingCircus',''});
    h2 = ini.SetValues('noedits', {'filter_done','artefacts_done','ground_done','median_done'}, {'False','False','False','False'});
    h3 = ini.SetValues('triggers', {'dead_file','dead_unit','ignore_times'}, {'SpykingCircus_artefacts_samples.dead','timestep','True'});
    if any([h1; h2; h3] ~= 1), error('Something went wrong with adjusting parameters'); end
    
    status = ini.WriteFile(fname_params);
    if status == false
        error('Couldn''t write file');
    end
    ini.ToString()
    
    % write params file 
    nb_channels         = size(cfg.circus.channel,2);
    writeProbeFile(nb_channels,fname_prb);
    
end

