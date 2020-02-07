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
    markerlist = [];
    for i = 1 : size(cfg.name,2)
        if ismember(cfg.name{i},cfg.name)
            markerlist = [markerlist, i];
        end
    end
    
    for ipart = 1:length(MuseStruct)
        
        for imarker = markerlist
            
            hasmarker = false(length(MuseStruct{ipart}),1);
            
            for idir = 1:length(MuseStruct{ipart}) %according to cfg.directorylist
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,(cfg.muse.startend{imarker,1}))
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                            
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

                                    Startsample(ievent) = ss + cfg.epoch.toi{imarker}(1) * cfg.LFP.resamplefs - cfg.epoch.pad(imarker) * cfg.LFP.resamplefs;
                                    Endsample(ievent)   = es + cfg.epoch.toi{imarker}(2) * cfg.LFP.resamplefs + cfg.epoch.pad(imarker) * cfg.LFP.resamplefs;
                                    Offset(ievent)      = (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad(imarker)) * cfg.LFP.resamplefs;
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
                            elseif isMicromed || isBrainvision
                                nfile = 1; %only one file with all electrodes
                            end
                            
                            for ifile = 1 : nfile
                                
                                %load data
                                if isNeuralynx
                                    temp                    = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.LFP.channel{ifile},'*.ncs']));
                                    fname{1}                = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                                    dat                     = ft_read_neuralynx_interp(fname);
                                elseif isMicromed
                                    dat                     = preprocessing_eeg_emg(cfg, ipart, idir, true);
                                end
                                
                                % filter with FT
                                cfgtemp                 = [];
                                cfgtemp.hpfilter        = cfg.LFP.hpfilter;
                                cfgtemp.hpfreq          = cfg.LFP.hpfreq;
                                dat                     = ft_preprocessing(cfgtemp,dat);
                                
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
            end % idir
            
            % concatinate data of different datasets (over trials)
            LFP{ipart}{imarker} = ft_appenddata([],dirdat{find(hasmarker)});
            clear dirdat*
            
            % add samplerate
            LFP{ipart}{imarker}.fsample = cfg.LFP.resamplefs;
            
        end % imarker
        
    end % ipart
end

if savedat
    save(fname_out,'LFP','-v7.3');
end
