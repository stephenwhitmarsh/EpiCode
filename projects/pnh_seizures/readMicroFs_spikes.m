function [dat_micro] = readMicroFs_spikes(cfg,MuseStruct_micro)

fname = fullfile(cfg.datasavedir,['data_aligned_microFs_',cfg.label,'.mat']);
if exist(fname,'file') && cfg.force == false
    fprintf('Loading precomputed data: %s \n',fname);
    load(fname,'dat_micro');    
else
    trialcount  = 1;       % to have a continuously trialcount
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
                    
                    % create Fieldtrip trl
                    hdr_micro   = ft_read_header(fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames(1).name));                    
                    Startsample     = [];
                    Endsample       = [];
                    Trialnr = [];
                    for ievent = 1 : size(MuseStruct_micro{idir}.markers.(cfg.startend{1}).events,2)
                        Startsample         = [Startsample; MuseStruct_micro{idir}.markers.(cfg.startend{1}).offset(ievent) - cfg.prestim  * hdr_micro.Fs]; % this can be replaced by the offset field
                        Endsample           = [Endsample;   MuseStruct_micro{idir}.markers.(cfg.startend{2}).offset(ievent) + cfg.poststim * hdr_micro.Fs];
                        Offset              = -ones(size(Endsample)) * cfg.prestim * hdr_micro.Fs;
                        Trialnr             = [Trialnr; trialcount];
                        trialcount          = trialcount + 1;
                    end
                    
                    trl                     = [Startsample, Endsample, Offset];
                    trl(:,4)                = ones(size(trl,1),1) * idir;
                    trl(:,5)                = 1:size(trl,1);            % try to find trials that are missing aftewards
                    trl(:,6)                = Startsample;                      % sampleinfo to cut up concatinated data later
                    trl(:,7)                = Endsample;                        % sampleinfo to cut up concatinated data later
                    trl(:,8)                = Endsample-Startsample+1;          % sampleinfo to cut up concatinated data later
                    trl(:,9)                = Offset;
                    trl(:,10)               = Trialnr;
                    trl                     = trl(Startsample > 0 & Endsample < hdr_micro.nSamples,:); % so not to read before BOF or after EOFs
                          
                    % load trials for selected MICRO channels
                    filecounter = 1;
                    for ifile = micro_filenrs
                        
                        cfgtemp                         = [];
                        cfgtemp.hpfilter                = cfg.hpfilter;
                        cfgtemp.hpfreq                  = cfg.hpfreq;
                        cfgtemp.dataset                 = fullfile(MuseStruct_micro{idir}.directory, MuseStruct_micro{idir}.filenames(ifile).name);
                        temp                            = ft_preprocessing(cfgtemp);
                        
                        cfgtemp.trl                     = trl;
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
