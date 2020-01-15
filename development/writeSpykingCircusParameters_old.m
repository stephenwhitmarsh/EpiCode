function writeSpykingCircusParameters(cfg)

for ipart = 1 : size(cfg,2)
    
    fname_prb               = fullfile(cfg.datasavedir,[cfg.prefix(1:end-1)],['p',num2str(ipart)],[cfg.prefix,'SpykingCircus.prb']);
    fname_params_default    = 'SpykingCircusDefaultSettings.params';
    [~, name, ~] = fileparts(cfg.fnames_ncs{ipart}{1});
    fname_params_patient    = fullfile(cfg.datasavedir,[name,'.params']);
    fname_params_patient    = fullfile(cfg.datasavedir,[cfg.prefix(1:end-1)],['p',num2str(ipart)],[name,'.params']);
 
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
    
    h1 = ini.SetValues('data', {'file_format','stream_mode','mapping','suffix','overwrite','output_dir','ncs_pattern'}, {'neuralynx','None',[cfg.prefix,'SpykingCircus.prb'],'','False','SpykingCircus',''});
    h2 = ini.SetValues('noedits', {'filter_done','artefacts_done','ground_done','median_done'}, {'False','False','False','False'});
    h3 = ini.SetValues('triggers', {'dead_file','dead_unit','ignore_times'}, {[cfg.prefix,'SpykingCircus_artefacts_samples.dead'],'timestep','True'});
    if any([h1; h2; h3] ~= 1), error('Something went wrong with adjusting parameters'); end
    
    ini.WriteFile(fname_params_patient);
    ini.ToString()
    
    % write params file
    
    nb_channels         = size(cfg.circus.channel,2);
    fid                 = fopen(fname_prb,'w+');
    
    fprintf(fid,'# Probe file for Spyking Circus \n');
    fprintf(fid,'# Created automatically by writeSpykingCircus.m \n \n');
    
    fprintf(fid,'total_nb_channels = %d;\n',nb_channels);
    fprintf(fid,'radius            = 10;\n\n');
    fprintf(fid,'channel_groups = {\n');
    fprintf(fid,'\t1: {\n');
    fprintf(fid,"\t\t'channels':range(0,%d),\n",nb_channels);
    fprintf(fid,"\t\t'geometry': {\n");
    for i = 1 : nb_channels
        fprintf(fid,'\t\t\t\t\t\t%d: [0, %d],\n',i-1,i*200);
    end
    fprintf(fid,'\t\t},\n');
    fprintf(fid,"\t\t'graph' : []\n");
    fprintf(fid,"\t}\n");
    fprintf(fid,"}\n");
    fclose(fid);
    
end