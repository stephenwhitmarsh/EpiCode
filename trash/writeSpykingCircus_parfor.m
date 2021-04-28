function cfg = writeSpykingCircus_parfor(cfg, MuseStruct, poolsize, force)

poolobj = parpool(poolsize);

fname_output = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_trialinfo_parts.mat']);

if exist(fname_output,'file') && force == false
    fprintf('\nLoading trialinfo: %s \n',fname_output);
    temp = load(fname_output,'cfg');
    cfg.sampleinfo          = temp.cfg.sampleinfo;
    cfg.deadfile_ms         = temp.cfg.deadfile_ms;
    cfg.deadfile_samples    = temp.cfg.deadfile_samples;
    cfg.fnames_ncs          = temp.cfg.fnames_ncs;
    
else
    
    % loop through different parts
    for ipart = 1 : size(MuseStruct,2)
        
        fprintf('\n*** Starting on part %d ***\n',ipart)
        
        % process channels separately 
        for ichan = 1 : size(cfg.circus.channel,2)
            
            % unslice variables for parfor
            channel     = cfg.circus.channel;
            reref       = cfg.circus.reref;
            rerefchan   = cfg.circus.refchan;
            hpfilter    = cfg.circus.hpfilter;
            
            % loop over all directories (time), concatinating channel
            parfor idir_parfor = 1:length(MuseStruct{ipart})
                
                temp                      = dir(fullfile(MuseStruct{ipart}{idir_parfor}.directory,['*',channel{ichan},'.ncs']));
                cfgtemp                   = [];
                cfgtemp.dataset           = fullfile(MuseStruct{ipart}{idir_parfor}.directory, temp.name);
                fprintf('LOADING %s\n',cfgtemp.dataset);
                dirdat{idir_parfor}              = ft_preprocessing(cfgtemp);
                
                if strcmp(reref,'yes')
                    temp                      = dir(fullfile(MuseStruct{ipart}{idir_parfor}.directory,['*',refchan,'.ncs']));
                    cfgtemp                   = [];
                    cfgtemp.dataset           = fullfile(MuseStruct{ipart}{idir_parfor}.directory, temp.name);
                    refdat                    = ft_preprocessing(cfgtemp);
                    dirdat{idir_parfor}.trial{1}     = dirdat{idir_parfor}.trial{1} - refdat.trial{1};
                    % clear refdat
                end
                
                % creation of artefacts caused by large offsets
                % should not happen with continuous data
                if strcmp(hpfilter,'yes')
                    cfgtemp                   = [];
                    cfgtemp.hpfilter          = hpfilter;
                    cfgtemp.hpfreq            = hpfreq;
                    dirdat{idir_parfor}              = ft_preprocessing(cfgtemp,dirdat{idir_parfor});
                end
                
                % truncate label to make them equal over files
                dirdat{idir_parfor}.label{1}     = dirdat{idir_parfor}.label{1}(end-6:end); % can be replaced by circus.channel
                
                % save sampleinfo to reconstruct data again after reading SC
                sampleinfo(idir_parfor,:)    = dirdat{idir_parfor}.sampleinfo;
                
            end % parfor
            cfg.sampleinfo = sampleinfo;
            
            % concatinate data over files
            chandat = dirdat{1};
            for idir = 2 : length(MuseStruct{ipart})
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
                chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
            end
            clear dirdat
            
            % write data in .ncs format
            temp                        = dir(fullfile(MuseStruct{ipart}{idir}.directory,['*',cfg.circus.channel{ichan},'.ncs']));
            hdr_micro                   = ft_read_header(fullfile(MuseStruct{ipart}{idir}.directory, temp.name));
            hdr                         = [];
            hdr.Fs                      = hdr_micro.Fs;
            hdr.nSamples                = size(chandat.trial{1},2);
            hdr.nSamplePre              = 0;
            hdr.nChans                  = 1;
            hdr.FirstTimeStamp          = 0;
            hdr.TimeStampPerSample      = hdr_micro.TimeStampPerSample;
            hdr.label                   = chandat.label;
            fname                       = fullfile(cfg.datasavedir,[cfg.prefix,'p',num2str(ipart),'-multifile-',chandat.label{1},'.ncs']);
            ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
            
            cfg.fnames_ncs{ipart}{ichan} = fname;
            clear chandat
        end % ichan
        
        %% write deadtime, i.e. artefact file for Spyking-Circus
        
        deadfile_ms         = [];
        deadfile_samples    = [];
        last_ms             = 0;
        last_samples        = 0;
        dirlist{ipart}      = [];
        
        for idir = 1 : size(MuseStruct{ipart},2)
            if isfield(MuseStruct{ipart}{idir},'markers')
                if isfield(MuseStruct{ipart}{idir}.markers,'BAD__START__')
                    if isfield(MuseStruct{ipart}{idir}.markers.BAD__START__,'synctime')
                        % check if there is an equal amount of start and end markers
                        if size(MuseStruct{ipart}{idir}.markers.BAD__START__.events,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.events,2) == 0
                            
                            fprintf('Great, recovered same number of start and end markers \n')
                            deadfile_ms         = [deadfile_ms;         MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                            deadfile_samples    = [deadfile_samples;    MuseStruct{ipart}{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{ipart}{idir}.markers.BAD__END__.offset'+last_samples];
                            hdr                 = ft_read_header(fullfile(MuseStruct{ipart}{idir}.directory,MuseStruct{ipart}{idir}.filenames{1}));
                            last_samples        = last_samples + hdr.nSamples;
                            last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                            
                        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.events,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.events,2) > 0
                            
                            fprintf('ERROR! more start than end found in %s \n',MuseStruct{ipart}{idir}.directory);
                            
                        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.events,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.events,2) < 0
                            
                            fprintf('ERROR! more end than start found in %s - CORRECTING \n',MuseStruct{ipart}{idir}.directory)
                            for i = 1 : 10
                                for itrial = 1 : length(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime)
                                    start(itrial) = MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(itrial);
                                    stop(itrial)  = MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(itrial);
                                end
                                
                                x = find(start > stop,1,'first');
                                if ~isempty(x)
                                    MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(x) = [];
                                    MuseStruct{ipart}{idir}.markers.BAD__END__.offset(x) = [];
                                    
                                end
                            end
                            
                            deadfile_ms         = [deadfile_ms;         MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                            deadfile_samples    = [deadfile_samples;    MuseStruct{ipart}{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{ipart}{idir}.markers.BAD__END__.offset'+last_samples];
                            
                            hdr                 = ft_read_header(fullfile(MuseStruct{ipart}{idir}.directory,MuseStruct{ipart}{idir}.filenames(1).name));
                            last_samples        = last_samples + hdr.nSamples;
                            last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                        end
                    else
                        fprintf('Found no artefacts for file: %s\n',MuseStruct{ipart}{idir}.directory);
                    end
                else
                    fprintf('Found no artefacts for file: %s\n',MuseStruct{ipart}{idir}.directory);
                end
            end
            dirlist{ipart} = [dirlist{ipart}; MuseStruct{ipart}{idir}.directory];
            fprintf('%d\n',idir);
        end % ipart
        
        % return info
        cfg.deadfile_ms{ipart}         = deadfile_ms;
        cfg.deadfile_samples{ipart}    = deadfile_samples;
        
        % write artefacts to txt file for spyking circus
        filename = fullfile(cfg.datasavedir,[cfg.prefix,'p',num2str(ipart),'-SpykingCircus_artefacts_ms.dead']);
        fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
        dlmwrite(filename,deadfile_ms,'delimiter','	','precision','%.4f');
        
        filename = fullfile(cfg.datasavedir,[cfg.prefix,'p',num2str(ipart),'-SpykingCircus_artefacts_samples.dead']);
        fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
        dlmwrite(filename,deadfile_samples,'delimiter','	','precision','%.4f');
        
        filename = fullfile(cfg.datasavedir,[cfg.prefix,'p',num2str(ipart),'-SpykingCircus_dirlist.txt']);
        fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
        dlmwrite(filename,dirlist{ipart});
        
        fid = fopen(filename,'w+');
        for r=1:size(dirlist,1)
            fprintf(fid,'%s\n',dirlist{ipart}(r,:));
        end
        fclose(fid);
    end
    
    % save trialinfo stuff
    save(fname_output,'cfg');
end

delete(poolobj);

