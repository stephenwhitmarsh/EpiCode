function [dat_micro, dat_macro] = readEEG(cfg,MuseStruct_micro,MuseStruct_macro,force)

%% Make overview of data as segmented by markers aligned to (first) peak

fname = fullfile(cfg.datasavedir,[cfg.prefix,'data_aligned.mat']);

if exist(fname,'file') && force == false
    fprintf('**************************************************\n');
    fprintf('*** Loading precomputed EEG data: %s ***\n',fname);
    fprintf('**************************************************\n\n');
    
    load(fname,'dat_micro');
    load(fname,'dat_macro');
    
else
    fprintf('**************************************************\n');
    fprintf('(re-) computing EEG data: %s ***\n',fname);
    fprintf('**************************************************\n\n');
   
    for imarker = 1 : size(cfg.name,2)
        hasmarker = zeros(length(MuseStruct_micro),1);
        for idir = 1:length(MuseStruct_micro)
            if isfield(MuseStruct_micro{idir},'markers')
                if isfield(MuseStruct_micro{idir}.markers,(cfg.muse.startend{imarker,1}))
                    if ~isempty(MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events)
                        
                        % select MICRO files
                        micro_filenrs = [];
                        for ifile = 1 : size(MuseStruct_micro{idir}.filenames,2)
                            for ilabel = 1 : size(cfg.labels.micro,2)
                                if ~isempty(strfind(MuseStruct_micro{idir}.filenames{ifile},cfg.labels.micro{ilabel}))
                                    micro_filenrs = [micro_filenrs, ifile];
                                end
                            end
                        end
                        
                        % select MACRO files
                        macro_filenrs = [];
                        for ifile = 1 : size(MuseStruct_macro{idir}.filenames,2)
                            for ilabel = 1 : size(cfg.labels.macro,2)
                                if ~isempty(strfind(MuseStruct_macro{idir}.filenames{ifile},cfg.labels.macro{ilabel}))
                                    macro_filenrs = [macro_filenrs, ifile];
                                end
                            end
                        end
                        
                        % load trials for selected MICRO channels
                        for ifile = micro_filenrs
                            
                            cfgtemp = [];
                            cfgtemp.dataset                 = fullfile(MuseStruct_micro{idir}.directory, MuseStruct_micro{idir}.filenames{ifile});
                            cfgtemp.hpfilter                = cfg.LFP.hpfilter; 
                            cfgtemp.hpfreq                  = cfg.LFP.hpfreq;
                            temp                            = ft_preprocessing(cfgtemp);
                            
                            % create Fieldtrip trl
                            hdr_micro = ft_read_header(fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames{ifile}));
                            Startsample     = [];
                            Endsample       = [];
                            for ievent = 1 : size(MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events,2)
                                Startsample         = [Startsample; MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(ievent) + cfg.epoch.toi{imarker}(1) * hdr_micro.Fs - cfg.epoch.pad(imarker) * hdr_micro.Fs];
                                Endsample           = [Endsample;   MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,2}).offset(ievent) + cfg.epoch.toi{imarker}(2) * hdr_micro.Fs + cfg.epoch.pad(imarker) * hdr_micro.Fs];
                                Offset              = ones(size(Endsample)) * (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad(imarker)) * hdr_micro.Fs;                                  
                            end
                            
                            cfgtemp                         = [];
                            cfgtemp.trl                     = [Startsample, Endsample, Offset];
                            cfgtemp.trl(:,4)                = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                            cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr_micro.nSamples,:); % so not to read before BOF or after EOFs
                            cfgtemp.dataset                 = fullfile(MuseStruct_micro{idir}.directory, MuseStruct_micro{idir}.filenames{ifile});
                            filedat_micro{ifile}            = ft_redefinetrial(cfgtemp,temp);
                            clear temp
                            
                            % downsample data and correct baseline
                            cfgtemp                         = [];
                            cfgtemp.resamplefs              = cfg.LFP.resamplefs;
                            if strcmp(cfg.LFP.baseline,'no')
                                cfgtemp.demean              = 'no';
                            else
                                cfgtemp.demean              = 'yes';
                                cfgtemp.baselinewindow      = cfg.LFP.baseline;
                            end
                            filedat_micro{ifile}   = ft_resampledata(cfgtemp,filedat_micro{ifile});
                            
                            % truncate label to make them equal over files
                            filedat_micro{ifile}.label{1} = filedat_micro{ifile}.label{1}(end-6:end);
                            
                         end
                        
                        % load trials for selected MACRO channels
                        for ifile = macro_filenrs
                            
                            cfgtemp = [];
                            cfgtemp.dataset                 = fullfile(MuseStruct_macro{idir}.directory, MuseStruct_macro{idir}.filenames{ifile});
                            temp                            = ft_preprocessing(cfgtemp);
                            if strcmp(cfg.LFP.hpfilter,'yes')
                                fprintf('Filtering macrodata %d of %d, this can take a while...',idir,length(MuseStruct_macro));
                                temp.trial{1}               = bandpassFilter(temp.trial{1},temp.fsample,cfg.LFP.hpfreq,cfg.LFP.resamplefs);
                            end
                            
                            % create Fieldtrip trl
                            hdr_macro       = ft_read_header(fullfile(MuseStruct_macro{idir}.directory,MuseStruct_macro{idir}.filenames{ifile}));
                            Startsample     = [];
                            Endsample       = [];
                            for ievent = 1 : size(MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).events,2)
                                
                                Startsample         = [Startsample; MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(ievent) + cfg.epoch.toi{imarker}(1) * hdr_macro.Fs - cfg.epoch.pad(imarker) * hdr_macro.Fs];
                                Endsample           = [Endsample;   MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,2}).offset(ievent) + cfg.epoch.toi{imarker}(2) * hdr_macro.Fs + cfg.epoch.pad(imarker) * hdr_macro.Fs];
                                Offset              = ones(size(Endsample)) * (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad(imarker)) * hdr_macro.Fs;
                                
                            end
                            
                            % load and segment data in trials
                            cfgtemp                         = [];
                            cfgtemp.trl                     = [Startsample, Endsample, Offset];
                            cfgtemp.trl(:,4)                = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                            cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr_macro.nSamples,:); % so not to read before BOF or after EOFs
                            cfgtemp.dataset                 = fullfile(MuseStruct_macro{idir}.directory, MuseStruct_macro{idir}.filenames{ifile});
                            cfgtemp.hpfilter                = cfg.LFP.hpfilter;
                            cfgtemp.hpfreq                  = cfg.LFP.hpfreq;
                            filedat_macro{ifile}            = ft_redefinetrial(cfgtemp,temp);
                            clear temp
                            
                            % downsample data and correct baseline
                            cfgtemp                         = [];
                            cfgtemp.resamplefs              = cfg.LFP.resamplefs;
                            if strcmp(cfg.LFP.baseline,'no')
                                cfgtemp.demean              = 'no';
                            else
                                cfgtemp.demean              = 'yes';
                                cfgtemp.baselinewindow      = cfg.LFP.baseline;
                            end
                            filedat_macro{ifile}   = ft_resampledata(cfgtemp,filedat_macro{ifile});
                            
                            % truncate label to make them equal over files
                            filedat_macro{ifile}.label{1} = filedat_macro{ifile}.label{1}(end-6:end);
                                          
                        end
                        % concatinate channels, separately for MICRO/MACRO
                        cfgtemp                             = [];
                        cfgtemp.keepsampleinfo              = 'no';
                        dirdat_micro{idir}                  = ft_appenddata(cfgtemp,filedat_micro{micro_filenrs});
                        dirdat_macro{idir}                  = ft_appenddata(cfgtemp,filedat_macro{macro_filenrs});
                        dirdat_micro{idir}.trialinfo(:,2)   = 1:size(dirdat_micro{idir}.trial,2);
                        dirdat_macro{idir}.trialinfo(:,2)   = 1:size(dirdat_macro{idir}.trial,2);
                        clear filedat*
                        
                        % flag for averaging
                        hasmarker(idir) = 1;
                    end
                end
            end
        end
        
        % concatinate data over trials
        dat_micro{imarker} = ft_appenddata([],dirdat_micro{find(hasmarker)});
        dat_macro{imarker} = ft_appenddata([],dirdat_macro{find(hasmarker)});
        
        % add samplerate
        dat_micro{imarker}.fsample = cfg.LFP.resamplefs;
        dat_macro{imarker}.fsample = cfg.LFP.resamplefs;
        
        clear dirdat*
    end
    
    save(fname,'dat_micro','dat_macro');
end
