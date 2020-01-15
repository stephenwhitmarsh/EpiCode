function [dat_micro] = readMicroFs_alldata(cfg,MuseStruct_micro)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'data_raw_microFs_.mat']);
if exist(fname,'file') && cfg.force == false
    fprintf('Loading precomputed data: %s \n',fname);
    load(fname,'dat_micro');    
else
    hasmarker = zeros(length(MuseStruct_micro),1);
    for idir = 1:length(MuseStruct_micro)
        if isfield(MuseStruct_micro{idir},'markers')
            if isfield(MuseStruct_micro{idir}.markers,(cfg.startend{1}))
                if ~isempty(MuseStruct_micro{idir}.markers.(cfg.startend{1}).events)
                    
                    % select MICRO files
                    micro_filenrs = [];
                    for ifile = 1 : size(MuseStruct_micro{idir}.filenames,1)
                        for ilabel = 1 : size(cfg.channel,2)
                            if ~isempty(strfind(MuseStruct_micro{idir}.filenames(ifile).name,cfg.channel{ilabel}))
                                micro_filenrs = [micro_filenrs, ifile];
                            end
                        end
                    end                 
               
                    % load trials for selected MICRO channels
                    filecounter = 1;
                    for ifile = micro_filenrs
                        
                        cfgtemp                         = [];
                        cfgtemp.dataset                 = fullfile(MuseStruct_micro{idir}.directory, MuseStruct_micro{idir}.filenames(ifile).name);
                        filedat_micro{ifile}            = ft_redefinetrial(cfgtemp,temp);
                        clear temp
                              
                        % truncate label to make them equal over files
                        filedat_micro{ifile}.label{1}   = filedat_micro{ifile}.label{1}(end-6:end);
                        
                        filecounter = filecounter + 1;
                    end
                    
                    % concatinate channels
                    dirdat_micro{idir}                  = ft_appenddata([],filedat_micro{micro_filenrs});
                    
                    % transfer trialinfo (from first file)
                    dirdat_micro{idir}.trialinfo        = filedat_micro{1}.trialinfo;
                    
                    clear filedat
                    
                    % flag for averaging
                    hasmarker(idir) = 1;
                end %if
            end %if
        end %if
    end %idir
    
    % concatinate data over trials
    dat_micro = ft_appenddata([],dirdat_micro{find(hasmarker)});
    dat_micro.fsample = hdr_micro.Fs;
    dat_micro.hdr = hdr_micro;
    clear dirdat*
    
    save(fname,'dat_micro','-v7.3');
end
