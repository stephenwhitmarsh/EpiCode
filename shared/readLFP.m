function [LFP] = readLFP(cfg,MuseStruct,force,savedat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dat_micro, dat_macro] = readLFP(cfg,MuseStruct_micro,MuseStruct_macro,force,savedat)
%
% Reads data from macro and micro electrodes, epoched according to markers extracted from Muse,
% and downsampled to same samplerate.
%
% Necessary fields:





% Note:
% Names of markers that contain a space (' ') or minus ('-') will be
% replaced by an underscore ('_').
%
% Dependencies: writeMuseMarkers.m, dir2.m, recent FieldTrip version
%
% (c) Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'data_aligned.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed LFP data ***\n');
    fprintf('************************************\n\n');
    
    load(fname_out,'LFP');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing LFP data ***\n');
    fprintf('********************************\n\n');
    
    %get format to adapt script for each format
    %specificities :
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
    
    
    % select those markers to load
    %     markerlist = [];
    %     for i = 1 : size(cfg.LFP.name,2)
    %         if ismember(cfg.LFP.name{i},cfg.name)
    %             markerlist = [markerlist, i];
    %         end
    %     end
    
    for ipart = 1:length(MuseStruct)
        
        %         for imarker = markerlist
        for imarker = 1 : size(cfg.LFP.name,2)
            
            
            fprintf('For marker %s\n',cfg.LFP.name{imarker});
            
            hasmarker = false(length(MuseStruct{ipart}),1);
            
            for idir = 1:length(MuseStruct{ipart}) %according to cfg.directorylist
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,(cfg.muse.startend{imarker,1}))
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                            if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                            
                            % create trial segmentation common to resampled
                            % data. Neuralynx : same markers for all files
                            % of one dir.
                            Startsample             = [];
                            Endsample               = [];
                            Stage                   = [];
                            Offset                  = [];
                            trialnr                 = [];
                            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)
                                
                                ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent) * cfg.LFP.resamplefs);
                                idx =  find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime * cfg.LFP.resamplefs) >= ss,1,'first');
                                es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(idx) * cfg.LFP.resamplefs);
                                
                                if ~isempty(es) % && (es - ss) * hdr_micro.Fs < 4 %% find a way to check for Paul's data
                                    
                                    Startsample(ievent) = ss + cfg.epoch.toi{imarker}(1) * cfg.LFP.resamplefs - cfg.epoch.pad{imarker} * cfg.LFP.resamplefs;
                                    Endsample(ievent)   = es + cfg.epoch.toi{imarker}(2) * cfg.LFP.resamplefs + cfg.epoch.pad{imarker} * cfg.LFP.resamplefs;
                                    Offset(ievent)      = (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad{imarker}) * cfg.LFP.resamplefs;
                                    trialnr(ievent)     = ievent;
                                    Stage(ievent)       = -1;
                                    
                                    % find overlap with hypnogram markers
                                    for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                                        if isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                                            for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                                                y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i) * cfg.LFP.resamplefs;
                                                y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i) * cfg.LFP.resamplefs;
                                                if (y1 < ss) && (ss < y2)
                                                    fprintf('Found "%s" overlapping with "%s" : adding to trialinfo: ',cfg.name{imarker},cell2mat(hyplabel));
                                                    switch cell2mat(hyplabel)
                                                        case 'PHASE_1'
                                                            fprintf('%d\n',1);
                                                            Stage(ievent) = 1;
                                                        case 'PHASE_2'
                                                            fprintf('%d\n',2);
                                                            Stage(ievent) = 2;
                                                        case 'PHASE_3'
                                                            fprintf('%d\n',3);
                                                            Stage(ievent) = 3;
                                                        case 'REM'
                                                            fprintf('%d\n',4);
                                                            Stage(ievent) = 4;
                                                        case 'AWAKE'
                                                            fprintf('%d\n',0);
                                                            Stage(ievent) = 0;
                                                        case 'NO_SCORE'
                                                            fprintf('%d\n',0);
                                                            Stage(ievent) = 0;
                                                        otherwise
                                                            error('Unexpected label name in Hypnogram\n');
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end % ~isempty(es)
                            end
                            
                            
                            
                            % loop over files
                            if isNeuralynx
                                nfile = size(cfg.LFP.channel,2); %one file per channel
                            elseif isMicromed
                                nfile = 1; %only one file with all electrodes
                                fname                             = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
                            elseif isBrainvision
                                nfile = 1; %only one file with all electrodes
                                fname                             = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
                            end
                            
                            for ifile = 1 : nfile
                                
                                %load data
                                if isNeuralynx
                                    temp                    = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.LFP.channel{ifile},'*.ncs']));
                                    fname{1}                = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                                    dat                     = ft_read_neuralynx_interp(fname);
                                    %hdr                     = ft_read_header(fname);
                                elseif isMicromed || isBrainvision
                                    cfgtemp.dataset                   = fname;
                                    cfgtemp.channel                   = cfg.labels.macro';
                                    dat                               = ft_preprocessing(cfgtemp);
                                end
                                
                                % Different preprocessing steps according
                                % to cfg
                                
                                if isfield(cfg.LFP, 'reref')
                                    if strcmp(cfg.LFP.reref, 'yes')
                                        cfgtemp                       = [];
                                        cfgtemp.reref                 = cfg.LFP.reref;
                                        cfgtemp.rerefmethod           = cfg.LFP.rerefmethod;
                                        cfgtemp.refchannel            = cfg.LFP.refchannel;
                                        dat                           = ft_preprocessing(cfgtemp,dat);
                                    end
                                end
                                
                                if isfield(cfg.LFP, 'bsfilter')
                                    if strcmp(cfg.LFP.bsfilter, 'yes')
                                        cfgtemp                       = [];
                                        cfgtemp.bsfilter              = cfg.LFP.bsfilter;
                                        cfgtemp.bsfreq                = cfg.LFP.bsfreq;
                                        dat                           = ft_preprocessing(cfgtemp,dat);
                                    end
                                end
                                
                                if isfield(cfg.LFP, 'hpfilter')
                                    if strcmp(cfg.LFP.hpfilter, 'yes')
                                        cfgtemp                       = [];
                                        cfgtemp.hpfilter              = cfg.LFP.hpfilter;
                                        cfgtemp.hpfreq                = cfg.LFP.hpfreq;
                                        cfgtemp.hpfiltord             = cfg.LFP.hpfiltord;
                                        cfgtemp.hpfilttype            = cfg.LFP.hpfilttype;
                                        dat                           = ft_preprocessing(cfgtemp,dat);
                                    end
                                end
                                
                                if isfield(cfg.LFP, 'lpfilter')
                                    if strcmp(cfg.LFP.lpfilter, 'yes')
                                        cfgtemp                       = [];
                                        cfgtemp.lpfilter              = cfg.LFP.lpfilter;
                                        cfgtemp.lpfreq                = cfg.LFP.lpfreq;
                                        cfgtemp.lpfilttype            = cfg.LFP.lpfilttype;
                                        dat                           = ft_preprocessing(cfgtemp,dat);
                                    end
                                end
                                
                                
                                
                                % Load EMG data if any, and append with EEG
                                % data
                                if isfield(cfg.LFP, 'emg')
                                    if ~isempty(cfg.LFP.emg)
                                        if ~strcmp(cfg.LFP.emg{imarker},'no')
                                            
                                            
                                            loadrefemg = false;
                                            if isfield(cfg.EMG, 'reref')
                                                if strcmp(cfg.EMG.reref, 'yes')
                                                    loadrefemg = true;
                                                end
                                            end
                                            
                                            if loadrefemg
                                                cfgtemp                   = [];
                                                cfgtemp.channel           = {cfg.LFP.emg{imarker}, cfg.EMG.refchannel{1}};%load the emg associated with eeg marker, and the ref
                                                cfgtemp.dataset           = fname;
                                                cfgtemp.reref             = cfg.EMG.reref;
                                                cfgtemp.rerefmethod       = cfg.EMG.rerefmethod;
                                                cfgtemp.refchannel        = cfg.EMG.refchannel;
                                                data_EMG                  = ft_preprocessing(cfgtemp);
                                            else
                                                cfgtemp                       = [];
                                                cfgtemp.channel               = cfg.LFP.emg{imarker};%load only the emg associated with eeg marker
                                                cfgtemp.dataset               = fname;
                                                data_EMG                      = ft_preprocessing(cfgtemp);
                                            end
                                            
                                            if isfield(cfg.EMG, 'bsfilter')
                                                if strcmp(cfg.EMG.bsfilter, 'yes')
                                                    cfgtemp                   = [];
                                                    cfgtemp.bsfilter          = cfg.EMG.bsfilter;
                                                    cfgtemp.bsfreq            = cfg.EMG.bsfreq;
                                                    cfgtemp.bsfiltord         = cfg.EMG.bsfiltord;
                                                    data_EMG                  = ft_preprocessing(cfgtemp, data_EMG);
                                                end
                                            end
                                            
                                            if isfield(cfg.EMG, 'hpfilter')
                                                if strcmp(cfg.EMG.hpfilter, 'yes')
                                                    cfgtemp                   = [];
                                                    cfgtemp.hpfilter          = cfg.EMG.hpfilter;
                                                    cfgtemp.hpfreq            = cfg.EMG.hpfreq;
                                                    data_EMG                  = ft_preprocessing(cfgtemp, data_EMG);
                                                end
                                            end
                                            
                                            %Appendata EEG EMG
                                            cfgtemp                       = [];
                                            cfgtemp.keepsampleinfo        = 'yes';
                                            dat                           = ft_appenddata(cfgtemp, dat, data_EMG);
                                            
                                        end
                                    end
                                end
                                
                                
                                
                                % downsample data and correct baseline
                                cfgtemp                         = [];
                                cfgtemp.resamplefs              = cfg.LFP.resamplefs;
                                if strcmp(cfg.LFP.baseline,'no')
                                    cfgtemp.demean              = 'no';
                                else
                                    cfgtemp.demean              = 'yes';
                                    cfgtemp.baselinewindow      = cfg.LFP.baselinewindow{imarker};
                                end
                                dat                             = ft_resampledata(cfgtemp,dat);
                                
                                % create Fieldtrip trl
                                cfgtemp                         = [];
                                cfgtemp.trl                     = round([Startsample; Endsample; Offset]');
                                cfgtemp.trl(:,4)                = trialnr;
                                cfgtemp.trl(:,6)                = idir;
                                cfgtemp.trl(:,7)                = Stage;
                                cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < length(dat.trial{1}),:); % so not to read before BOF or after EOFs
                                filedat{ifile}                  = ft_redefinetrial(cfgtemp,dat);
                                clear dat
                                
                                if isNeuralynx
                                    % same label over files
                                    filedat{ifile}.label{1}         = cfg.LFP.channel{ifile};
                                end
                                
                                % flag for averaging
                                hasmarker(idir)                 = true;
                            end
                            
                            % concatinate channels
                            cfgtemp                             = [];
                            cfgtemp.keepsampleinfo              = 'no';
                            dirdat{idir}                        = ft_appenddata(cfgtemp,filedat{:});
                            clear filedat*
                            
                            end
                        end
                    end
                end
            end % idir
            
            if exist('dirdat','var') %in case there is no marker in the data
                
                % concatinate data of different datasets (over trials)
                LFP{ipart}{imarker} = ft_appenddata([],dirdat{find(hasmarker)});
                clear dirdat*
                
                % add samplerate
                LFP{ipart}{imarker}.fsample = cfg.LFP.resamplefs;
                
            else
                LFP{ipart}{imarker} = [];
                fprintf('%s part %d : No data with marker ''%s''\n',cfg.prefix(1:end-1), ipart, cfg.LFP.name{imarker});
            end
            
        end % imarker
        
    end % ipart
end

if savedat
    save(fname_out,'LFP','-v7.3');
end
