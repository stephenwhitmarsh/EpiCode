function cfg = writeSpykingCircus(cfg, MuseStruct, force)

fname_output = fullfile(cfg.datasavedir,[cfg.prefix,'alldata_trialinfo.mat']);

% load trialinfo to reconstruct timeline of spike data
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

            temp                      = dir(fullfile(MuseStruct{idir}.directory,['*',cfg.circus.channel{ichan},'.ncs']));
            cfgtemp                   = [];
            cfgtemp.dataset           = fullfile(MuseStruct{idir}.directory, temp.name);
            dirdat{idir}              = ft_preprocessing(cfgtemp);
            
            if strcmp(cfg.circus.reref,'yes')
                temp                      = dir(fullfile(MuseStruct{idir}.directory,['*',cfg.circus.refchan,'.ncs']));
                cfgtemp                   = [];
                cfgtemp.dataset           = fullfile(MuseStruct{idir}.directory, temp.name);
                refdat                    = ft_preprocessing(cfgtemp);
                dirdat{idir}.trial{1}     = dirdat{idir}.trial{1} - refdat.trial{1};
                clear refdat
            end
            
            % creation of false spikes caused by large offsets
            % should not happen with continuous data
            % cfgtemp                   = [];
            % cfgtemp.hpfilter          = 'yes';
            % cfgtemp.hpfreq            = 100;
            % dirdat{idir}              = ft_preprocessing(cfgtemp,dirdat{idir});
            
            % truncate label to make them equal over files
            dirdat{idir}.label{1}     = dirdat{idir}.label{1}(end-6:end);
            
            % save sampleinfo to reconstruct data again after reading SC
            cfg.sampleinfo(idir,:)    = dirdat{idir}.sampleinfo;
            
        end
        
        % concatinate data over files
        chandat = dirdat{1};
        for idir = 2 : length(MuseStruct)
            fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
            chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
            chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
        end
        clear dirdat
        
        % write data in .ncs format
        hdr_micro                   = ft_read_header(fullfile(MuseStruct{1}.directory,MuseStruct{1}.filenames{ichan}));
        hdr                         = [];
        hdr.Fs                      = hdr_micro.Fs;
        hdr.nSamples                = size(chandat.trial{1},2);
        hdr.nSamplePre              = 0;
        hdr.nChans                  = 1;
        hdr.FirstTimeStamp          = 0;
        hdr.TimeStampPerSample      = hdr_micro.TimeStampPerSample;
        hdr.label                   = chandat.label;
        fname                       = fullfile(cfg.datasavedir,[cfg.prefix,'all_data_',chandat.label{1},'.ncs']);
        ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
        
        cfg.fnames_ncs{ichan} = fname;  
        clear chandat
    end % ichan
    
    %% write deadtime
    deadfile_ms         = [];
    deadfile_samples    = [];
    last_ms             = 0;
    last_samples        = 0;
    dirlist             = [];
    
    for idir = 1 : size(MuseStruct,2)
        if isfield(MuseStruct{idir},'markers')
            if isfield(MuseStruct{idir}.markers,'BAD__START__')
                if isfield(MuseStruct{idir}.markers.BAD__START__,'synctime')
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
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'all_data_artefacts_ms.dead']);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,deadfile_ms,'delimiter','	','precision','%.4f');
    
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'all_data_artefacts_samples.dead']);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,deadfile_samples,'delimiter','	','precision','%.4f');
    
    filename = fullfile(cfg.datasavedir,[cfg.prefix,'all_data_dirlist.txt']);
    fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
    dlmwrite(filename,dirlist);
    
    fid = fopen(filename,'w');
    for r=1:size(dirlist,1)
        fprintf(fid,'%s\n',dirlist(r,:));
    end
    fclose(fid);
    
end
