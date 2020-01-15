function cfg = writeSpykingCircus(cfg, MuseStruct, force, varargin)

if isempty(varargin)
    write = true;
else
    write = varargin{1};
end

fname_output = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_trialinfo.mat']);

if exist(fname_output,'file') && force == false
    fprintf('Loading trialinfo: %s \n',fname_output);
    temp = load(fname_output,'cfg');
    cfg.sampleinfo          = temp.cfg.sampleinfo;
    cfg.deadfile_ms         = temp.cfg.deadfile_ms;
    cfg.deadfile_samples    = temp.cfg.deadfile_samples;
    cfg.fnames_ncs          = temp.cfg.fnames_ncs;
    
else
    
    % count MICRO channels in first directory x cfg.label
    %     micro_filenrs = [];
    %     for ichan = 1 : size(MuseStruct{1}.filenames,1)
    %         for ilabel = 1 : size(cfg.circus.channel,2)
    %             if ~isempty(strfind(MuseStruct{1}.filenames{ichan},cfg.circus.channel{ilabel}))
    %                 micro_filenrs = [micro_filenrs, ichan];
    %             end
    %         end
    %     end
    %
    % process channels separately
    for ichan = 1 : size(cfg.circus.channel,2)
        
        % loop over all directories (time), reading one channel
        for idir = 1:length(MuseStruct)
            
            % FIXME due to the way the loop is set up, rereferencing takes
            % twice as long as needed
            
            if write
                
                % find full filename according to channel name
                temp                      = dir(fullfile(MuseStruct{idir}.directory,['*',cfg.circus.channel{ichan},'.ncs']));
                cfgtemp                   = [];
                cfgtemp.dataset           = fullfile(MuseStruct{idir}.directory, temp.name);
                fprintf('Loading data from %s\n',cfgtemp.dataset);
                
                % load data
                clear fname
                fname{1}                  = cfgtemp.dataset;
                dirdat{idir}              = ft_read_neuralynx_interp(fname);
                %                 dirdat{idir}              = ft_preprocessing(cfgtemp);
                
                % rereference if needed
                if strcmp(cfg.circus.reref,'yes')
                    temp                      = dir(fullfile(MuseStruct{idir}.directory,['*',cfg.circus.refchan,'.ncs']));
                    cfgtemp                   = [];
                    cfgtemp.dataset           = fullfile(MuseStruct{idir}.directory, temp.name);
                    clear fname
                    fname{1}                  = cfgtemp.dataset;
                    fprintf('Loading data from %s\n',cfgtemp.dataset);
                    
                    dirdat{idir}              = ft_read_neuralynx_interp(fname);
                    %                     refdat                    = ft_preprocessing(cfgtemp);
                    dirdat{idir}.trial{1}     = dirdat{idir}.trial{1} - refdat.trial{1};
                    clear refdat
                end
                
                % prevent creation of artefacts caused by large offsets
                % should not happen with continuous data
                if strcmp(cfg.circus.hpfilter,'yes')
                    cfgtemp                   = [];
                    cfgtemp.hpfilter          = cfg.circus.hpfilter;
                    cfgtemp.hpfreq            = cfg.circus.hpfreq;
                    dirdat{idir}              = ft_preprocessing(cfgtemp,dirdat{idir});
                end
                
                % truncate label to make them equal over 2-hour files
                dirdat{idir}.label{1}     = dirdat{idir}.label{1}(end-6:end);
                
            end % if write
            
            % create sampleinfo (nr of samples in file)
            temp                    = dir(fullfile(MuseStruct{idir}.directory,['*',cfg.circus.channel{ichan},'.ncs']));
            fprintf('Loading header from %s\n',fullfile(MuseStruct{idir}.directory, temp.name));

            hdr                     = ft_read_header(fullfile(MuseStruct{idir}.directory, temp.name));
            
            % save sampleinfo to reconstruct data again after reading SC
            cfg.sampleinfo(idir,:)  = [hdr.nSamplesPre hdr.nSamples];
            
        end
        
        % concatinate data over files
        if write
            chandat = dirdat{1};
            for idir = 2 : length(MuseStruct)
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
                chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
            end
        end
        
        % create filename for concatinated data - based on first directory
        temp                        = dir(fullfile(MuseStruct{1}.directory,['*',cfg.circus.channel{ichan},'.ncs']));
        fprintf('Loading header from %s\n',fullfile(MuseStruct{1}.directory, temp.name));
        hdrtemp                     = ft_read_header(fullfile(MuseStruct{1}.directory, temp.name));
        fname                       = fullfile(cfg.datasavedir,[cfg.prefix,'multifile-',cfg.labels.micro{ichan},'.ncs']);
        
        % save filenames to cfg (output of function)
        cfg.fnames_ncs{ichan}       = fname;
        
        % write data in .ncs format
        if write
            hdr                         = [];
            hdr.Fs                      = hdrtemp.Fs;
            hdr.nSamples                = size(chandat.trial{1},2);
            hdr.nSamplePre              = 0;
            hdr.nChans                  = 1;
            hdr.FirstTimeStamp          = 0;
            hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
            hdr.label                   = chandat.label;
            ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
        end
        clear fname
        clear chandat
        clear dirdat
        
    end % ichan
    
    %     if write
    %% write deadtime
    deadfile_ms         = [];
    deadfile_samples    = [];
    last_ms             = 0;
    last_samples        = 0;
    dirlist             = [];
    
    for idir = 1 : size(MuseStruct,2)
        if isfield(MuseStruct{idir},'markers')
            if isfield(MuseStruct{idir}.markers,'BAD__START__')
                if isfield(MuseStruct{idir}.markers.BAD__START__,'events')
                    % check if there is an equal amount of start and end markers
                    if size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) == 0
                        fprintf('Great, recovered same number of start and end markers \n')
                        
                        deadfile_ms         = [deadfile_ms;         MuseStruct{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                        deadfile_samples    = [deadfile_samples;    MuseStruct{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{idir}.markers.BAD__END__.offset'+last_samples];
                        hdr                 = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames{1}));
                        last_samples        = last_samples + hdr.nSamples;
                        last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                    elseif size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) > 0
                        fprintf('ERROR! more start than end found in %s \n',MuseStruct{idir}.directory);
                    elseif size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) < 0
                        fprintf('ERROR! more end than start found in %s - CORRECTING \n',MuseStruct{idir}.directory)
                        for i = 1 : 10
                            for itrial = 1 : length(MuseStruct{idir}.markers.BAD__START__.synctime)
                                start(itrial) = MuseStruct{idir}.markers.BAD__START__.synctime(itrial);
                                stop(itrial)  = MuseStruct{idir}.markers.BAD__END__.synctime(itrial);
                            end
                            
                            x = find(start > stop,1,'first');
                            if ~isempty(x)
                                MuseStruct{idir}.markers.BAD__END__.synctime(x) = [];
                                MuseStruct{idir}.markers.BAD__END__.offset(x) = [];
                                
                            end
                        end
                        
                        deadfile_ms         = [deadfile_ms;         MuseStruct{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                        deadfile_samples    = [deadfile_samples;    MuseStruct{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{idir}.markers.BAD__END__.offset'+last_samples];
                        
                        hdr                 = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(1).name));
                        last_samples        = last_samples + hdr.nSamples;
                        last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                    end
                else
                    fprintf('Found no artefacts for file: %s\n',MuseStruct{idir}.directory);
                end
            else
                fprintf('Found no artefacts for file: %s\n',MuseStruct{idir}.directory);
            end
        end
        dirlist = [dirlist; MuseStruct{idir}.directory];
        fprintf('%d\n',idir);
    end
    
    % return info
    cfg.deadfile_ms         = deadfile_ms;
    cfg.deadfile_samples    = deadfile_samples;
    
    % save data
    save(fname_output,'cfg');
    
    % write to analysis directory
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_artefacts_ms.dead']);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,deadfile_ms,'delimiter','	','precision','%.4f');
    
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_artefacts_samples.dead']);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,deadfile_samples,'delimiter','	','precision','%.0f');
    
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_dirlist.txt']);
    fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,dirlist);
    
    fid = fopen(filename,'w+');
    for r=1:size(dirlist,1)
        fprintf(fid,'%s\n',dirlist(r,:));
    end
    fclose(fid);
    %     end
end
