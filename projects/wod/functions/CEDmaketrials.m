function [CEDtrials] = CEDmaketrials(cfg,CEDStruct,force,savedat)


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
% ET PLUS
%
% Dependencies: writeMuseMarkers.m, dir2.m, recent FieldTrip version
%
% (c) Paul Baudin
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'data_trials.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed LFP data ***\n');
    fprintf('************************************\n\n');
    
    load(fname_out,'LFP');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing LFP data ***\n');
    fprintf('********************************\n\n');
    
    for ipart = 1:length(CEDStruct)
        
        %         for imarker = markerlist
        for imarker = 1 : size(cfg.LFP.name,2)
            
            
            fprintf('For marker %s\n',cfg.LFP.name{imarker});
            
            hasmarker = false(length(CEDStruct{ipart}),1);
            
            for idir = 1:length(CEDStruct{ipart}) %according to cfg.directorylist
                if isfield(CEDStruct{ipart}{idir},'markers')
                    if isfield(CEDStruct{ipart}{idir}.markers,(cfg.startend{imarker,1}))
                        if isfield(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,1}),'synctime')
                            if ~isempty(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,1}).synctime)
                            
                            
                            
                            % loop over files
                                nfile = size(cfg.LFP.channel,2); %one file per channel
 
                            
                            for ichannel = 1 : size(cfg.LFP.channel{imarker},2)
                                
                                
                                %load data
                                dat = CEDreadcontinous(cfg,cfg.LFP.channel{imarker}{ichannel}, ipart,idir);
                                  
                                % Different preprocessing steps according to cfg
                                
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
                                                              
                                
                                % downsample data 
                                % be carefull that if data are not
                                % resampled, all channels of all dirs MUST
                                % have the same sampling frequency
                                if cfg.LFP.doresample{imarker}
                                    cfgtemp                         = [];
                                    cfgtemp.resamplefs              = cfg.LFP.resamplefs;
                                    dat                             = ft_resampledata(cfgtemp,dat);
                                    trialsFs                     	= cfg.LFP.resamplefs;
                                else
                                    trialsFs                        = dat.fsample;
                                end
                                
                                % correct baseline
                                if strcmp(cfg.LFP.baseline,'yes')
                                    cfgtemp.demean              = 'yes';
                                    cfgtemp.baselinewindow      = cfg.LFP.baselinewindow{imarker};
                                    dat                         = ft_preprocessing(cfgtemp,dat);
                                end
                       
                                % create Fieldtrip trl
                                
                                
                                Startsample             = [];
                                Endsample               = [];
                                Stage                   = [];
                                Offset                  = [];
                                trialnr                 = [];
                                for ievent = 1 : size(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,1}).synctime,2)
                                    
                                    ss  = round(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,1}).synctime(ievent) * trialsFs);
                                    idx =  find(round(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,2}).synctime * trialsFs) >= ss,1,'first');
                                    es  = round(CEDStruct{ipart}{idir}.markers.(cfg.startend{imarker,2}).synctime(idx) * trialsFs);
                                    
                                    if ~isempty(es) % && (es - ss) * hdr_micro.Fs < 4 %% find a way to check for Paul's data
                                        
                                        Startsample(ievent) = ss + cfg.epoch.toi{imarker}(1) * trialsFs - cfg.epoch.pad{imarker} * trialsFs;
                                        Endsample(ievent)   = es + cfg.epoch.toi{imarker}(2) * trialsFs + cfg.epoch.pad{imarker} * trialsFs;
                                        Offset(ievent)      = (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad{imarker}) * trialsFs;
                                        trialnr(ievent)     = ievent;
                                        
                                    end
                                end
                                
                                cfgtemp                         = [];
                                cfgtemp.trl                     = round([Startsample; Endsample; Offset]');
                                cfgtemp.trl(:,4)                = trialnr;
                                cfgtemp.trl(:,5)                = idir;
                                cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < length(dat.trial{1}),:); % so not to read before BOF or after EOFs
                                filedat{ichannel}               = ft_redefinetrial(cfgtemp,dat);
                                clear dat
                                
                
                                % flag for averaging
                                hasmarker(idir)                 = true;
                            end %ichannel
                            
                            % concatinate channels
                            cfgtemp                             = [];
                            cfgtemp.keepsampleinfo              = 'no';
                            dirdat{idir}                        = ft_appenddata(cfgtemp,filedat{:});
                            clear filedat
                            
                            end
                        end
                    end
                end
            end % idir
            
            if exist('dirdat','var') %in case there is no marker in the data
                
                % concatinate data of different datasets (over trials)
                CEDtrials{ipart}{imarker} = ft_appenddata([],dirdat{hasmarker});
                clear dirdat*
                
                % add samplerate
                CEDtrials{ipart}{imarker}.fsample = trialsFs;
                
            else
                CEDtrials{ipart}{imarker} = [];
                fprintf('%s part %d : No data with marker ''%s''\n',cfg.prefix(1:end-1), ipart, cfg.LFP.name{imarker});
            end
            
        end % imarker
        
    end % ipart
end

if savedat
    save(fname_out,'CEDtrials','-v7.3');
end
