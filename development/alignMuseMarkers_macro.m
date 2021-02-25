function [MuseStruct_micro_aligned,MuseStruct_macro_aligned] = alignMuseMarkers_macro(cfg, MuseStruct_micro, MuseStruct_macro, force)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_micro_aligned.mat']);
if exist(fname,'file') && force == false
    fprintf('*******************************\n');
    fprintf('** Loading results alignment **\n');
    fprintf('*******************************\n\n');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_micro_aligned.mat']),'MuseStruct_micro_aligned');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_macro_aligned.mat']),'MuseStruct_macro_aligned');
else
    
    if force == true
        fprintf('*********************************\n');
        fprintf('** Forced redoing of alignment **\n');
        fprintf('*********************************\n\n');
    else
        fprintf('**************\n');
        fprintf('** Aligning **\n');
        fprintf('**************\n\n');
    end
    
    % initialize
    MuseStruct_micro_aligned = MuseStruct_micro;
    MuseStruct_macro_aligned = MuseStruct_macro;
    
    for imarker = 1 : size(cfg.name,2)
        
        % find data directories that have the required event
        markerindx = [];
        for idir = 1 : size(MuseStruct_micro,2)
            if isfield(MuseStruct_micro{idir},'markers')
                if isfield(MuseStruct_micro{idir}.markers,cfg.muse.startend{imarker,1})
                    if ~isempty(MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events)
                        fprintf('Found event timings for "%s" in %s\n',cfg.muse.startend{imarker,1},MuseStruct_micro{idir}.directory);
                        markerindx = [markerindx idir];
                    end
                end
            end
        end
        if isempty(markerindx)
            fprintf('Did not find any events for marker %s\n',cfg.name{imarker});
        end
        
        % find datafilename corresponding to channel
        channelindx = [];
        for ifile = 1 : size(MuseStruct_micro{idir}.filenames,2)
            if strfind(MuseStruct_micro{idir}.filenames{ifile},cfg.align.channel{imarker})
                channelindx = [channelindx ifile];
                fprintf('Found channel with pattern "%s" in %s\n',cfg.align.channel{imarker},MuseStruct_micro{idir}.filenames{channelindx});
                channeltype = 'micro';
            end
        end
        for ifile = 1 : size(MuseStruct_macro{idir}.filenames,2)
            if strfind(MuseStruct_macro{idir}.filenames{ifile},cfg.align.channel{imarker})
                channelindx = [channelindx ifile];
                fprintf('Found channel with pattern "%s" in %s\n',cfg.align.channel{imarker},MuseStruct_macro{idir}.filenames{channelindx});
                channeltype = 'macro';
            end
        end
        if isempty(channelindx)
            fprintf('Did not find any events for marker %s\n',cfg.align.channel{imarker});
        end
        if size(channelindx) > 1
            fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.align.channel{imarker});
        end
        
        % loop over all data directories that have the marker
        for idir = markerindx
            
            fprintf('Aligning data file %s\n',fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames{channelindx}));
            
            if strcmp(channeltype,'micro')
                cfgtemp             = [];
                cfgtemp.dataset     = fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames{channelindx});
                dat_orig            = ft_preprocessing(cfgtemp);
            else
                cfgtemp             = [];
                cfgtemp.dataset     = fullfile(MuseStruct_macro{idir}.directory,MuseStruct_macro{idir}.filenames{channelindx});
                dat_orig            = ft_preprocessing(cfgtemp);
            end
            
            dat_sel = dat_orig;
            
            % flip
            if strcmp(cfg.align.flip{imarker},'yes')
                dat_sel.trial{1} = -dat_sel.trial{1};
            end
            
            % abs
            if strcmp(cfg.align.abs{imarker},'yes')
                dat_sel.trial{1} = abs(dat_sel.trial{1});
            end
            
            % filtering data of selected channel for peak detection
            switch cfg.align.filter{imarker}
                case 'lp'
                    fprintf('Lowpass Filtering data for alignment. This can take a while \n');
                    cfgtemp             = [];
                    cfgtemp.lpfilter    = 'yes';
                    cfgtemp.lpfreq      = cfg.align.freq{imarker};
                    cfgtemp.lpfilttype  = 'fir';
                    dat_filt            = ft_preprocessing(cfgtemp,dat_sel);
                    %                 dat_filt.trial{1}   = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,1.0,cfg.lpfreq);
                case 'bp'
                    fprintf('Bandpass Filtering the data for alignment. This can take a while \n');
                    dat_filt            = dat_sel;
                    dat_filt.trial{1}   = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.align.freq{imarker}(1),cfg.align.freq{imarker}(2));
                    dat_filt.trial{1}   = envelope(dat_filt.trial{1});
                case 'none'
                    dat_filt            = dat_sel;

            end
            if strcmp(cfg.align.hilbert{imarker},'yes')
                fprintf('Converting to Hilbert envelope.\n');
                dat_filt.trial{1}   = envelope(dat_filt.trial{1});
            end
            
            % segment in trials around events - specific and only for the search
            cfgtemp             = [];
            cfgtemp.trl         = [MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).offset'+cfg.align.toiplot{imarker}(1)*dat_filt.fsample, ...
                                   MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).offset'+cfg.align.toiplot{imarker}(2)*dat_filt.fsample];
            cfgtemp.trl(:,3)    = ones(size(cfgtemp.trl,1),1)*+cfg.align.toiplot{imarker}(1)*dat_filt.fsample;
            cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
            cfgtemp.trl         = round(cfgtemp.trl);
            dat_filt_trl        = ft_redefinetrial(cfgtemp,dat_filt);
            dat_sel_trl         = ft_redefinetrial(cfgtemp,dat_sel);
            dat_orig_trl        = ft_redefinetrial(cfgtemp,dat_orig);
            
            % peak detection based on baseline and active period
            for itrial = 1 : size(dat_filt_trl.trial,2)
                
                % time indexes for baseline and active period
                t1_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.align.toibaseline{imarker}(1),1,'first');
                t2_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.align.toibaseline{imarker}(2),1,'last');
                t1_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.align.toiactive{imarker}(1),1,'first');
                t2_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.align.toiactive{imarker}(2),1,'last');
                
                % find peaks in active period
                [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));         
                
                % find peaks in baseline period
                [peaks_bl{itrial},~]               = findpeaks(dat_filt_trl.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
                if isempty(peaks_bl{itrial})
                    peaks_bl{itrial} = 0;
                end
                
                % select maximum peak in baseline
                if ~isempty(peaks_bl{itrial})
                    max_peaks_bl(itrial)           = max(peaks_bl{itrial});
                else
                    max_peaks_bl(itrial)           = nan;
                end
            end
            
            haspeak     = ones(1,size(dat_filt_trl.trial,2));
            
            for itrial = 1 : size(dat_filt_trl.trial,2)
                
                % threshold based on median of max peaks in all baseline periods
                peaks_sel_avg               = peaks_ac{itrial} > nanmean(max_peaks_bl) * cfg.align.thresh(imarker);
                peaks_ac_sel_avg{itrial}    = peaks_ac{itrial}(peaks_sel_avg);
                locs_ac_sel_avg{itrial}     = locs_ac{itrial}(peaks_sel_avg);
                
                % threshold based on median of max peaks in trail-by-trial baseline
                peaks_sel_trl               = peaks_ac{itrial} > max(peaks_bl{itrial}) * cfg.align.thresh(imarker);
                peaks_ac_sel_trl{itrial}    = peaks_ac{itrial}(peaks_sel_trl);
                locs_ac_sel_trl{itrial}     = locs_ac{itrial}(peaks_sel_trl);
                
                if isempty(locs_ac_sel_avg{itrial})
                    haspeak(itrial) = 0;
                    fprintf('Could not find peak in trial number %d. \n',itrial);
                end
            end
            
            dat_sel_aligned = dat_sel_trl;
            dat_orig_aligned = dat_orig_trl;
            hdr_macro       = ft_read_header(fullfile(MuseStruct_macro{idir}.directory,MuseStruct_macro{idir}.filenames{1})); % base on first file
            
            for itrial = 1 : size(dat_filt_trl.trial,2)
                if haspeak(itrial)
                    
                    % find relevant peak according to method 
                    if strcmp(cfg.align.method{imarker},'nearest')
                        [~, ip(itrial)] = min(abs(dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1))); % peak-time closest to zero
                    elseif strcmp(cfg.align.method{imarker},'max')
                        [~, ip(itrial)] = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1)); % max peak              
                    elseif strcmp(cfg.align.method{imarker},'first')
                        ip(itrial) = 1; % first peak after start search toi
                    end
                    
                    % align time axis to relevant peak
                    timeshift                       = dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                    dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                    dat_orig_aligned.time{itrial}   = dat_orig_trl.time{itrial} - timeshift;
    
                    sampleshift_micro               = round(timeshift*dat_filt_trl.fsample); % round to number of samples to shift
                    MuseStruct_micro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).offset(itrial)         = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(itrial) + sampleshift_micro; % note different direction than time
                    MuseStruct_micro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).samplehift(itrial)    = sampleshift_micro;
                    MuseStruct_micro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                    
                    sampleshift_macro               = round(timeshift*hdr_macro.Fs); % round to number of samples to shift
                    MuseStruct_macro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).offset(itrial)         = MuseStruct_macro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(itrial) + sampleshift_macro; % note different direction than time
                    MuseStruct_macro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).samplehift(itrial)    = sampleshift_macro;
                    MuseStruct_macro_aligned{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                    
                end
            end
            
            hasartefact     = zeros(1,size(dat_filt_trl.trial,2));
            
            if isfield(MuseStruct_micro{idir}.markers,'BAD__START__')
                % check if there is an equal amount of start and end markers
                if size(MuseStruct_micro{idir}.markers.BAD__START__.events,2)-size(MuseStruct_micro{idir}.markers.BAD__END__.events,2) == 0
                    fprintf('Great, recovered same number of starts and end markers \n')
                    for itrial = 1 : size(dat_filt_trl.trial,2)
                        fprintf('Checking trialnr. %d for artefacts... ',itrial);
                        trialinterval = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(itrial) : MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,2}).offset(itrial);
                        
                        for iart = 1 : size(MuseStruct_micro{idir}.markers.BAD__START__.events,2)
                            artinterval = MuseStruct_micro{idir}.markers.BAD__START__.offset(iart) : MuseStruct_micro{idir}.markers.BAD__END__.offset(iart);
                            if intersect(trialinterval,artinterval)
                                fprintf('Found artefact %d \n',iart);
                                hasartefact(itrial) = 1;
                            end
                        end
                        if hasartefact(itrial) == 0
                            fprintf('Ok\n');
                        end
                    end
                    
                elseif size(MuseStruct_micro{idir}.markers.BAD__START__.events,2)-size(MuseStruct_micro{idir}.markers.BAD__END__.events,2) > 0
                    fprintf('ERROR! more starts than ends found in %s \n',MuseStruct_micro{idir}.directory);
                elseif size(MuseStruct_micro{idir}.markers.BAD__START__.events,2)-size(MuseStruct_micro{idir}.markers.BAD__END__.events,2) < 0
                    fprintf('ERROR! more ends than starts found in %s - CORRECTING \n',MuseStruct_micro{idir}.directory)
                    
                    for i = 1 : 10
                        for itrial = 1 : length(MuseStruct_micro{idir}.markers.BAD__START__.synctime)
                            start(itrial) = MuseStruct_micro{idir}.markers.BAD__START__.synctime(itrial);
                            stop(itrial)  = MuseStruct_micro{idir}.markers.BAD__END__.synctime(itrial);
                        end
                        
                        x = find(start > stop,1,'first');
                        if ~isempty(x)
                            MuseStruct_micro{idir}.markers.BAD__END__.synctime(x) = [];
                            MuseStruct_micro{idir}.markers.BAD__END__.offset(x) = [];
                            
                        end
                    end
                    
                    for itrial = 1 : size(dat_filt_trl.trial,2)
                        fprintf('Checking trialnr. %d for artefacts... ',itrial);
                        trialinterval = MuseStruct_micro{idir}.markers.(cfg.muse.startend{1}).offset(itrial) : MuseStruct_micro{idir}.markers.(cfg.muse.startend{2}).offset(itrial);
                        
                        for iart = 1 : size(MuseStruct_micro{idir}.markers.BAD__START__.events,2)
                            artinterval = MuseStruct_micro{idir}.markers.BAD__START__.offset(iart) : MuseStruct_micro{idir}.markers.BAD__END__.offset(iart);
                            if intersect(trialinterval,artinterval)
                                fprintf('Found artefact %d \n',iart);
                                hasartefact(itrial) = 1;
                            end
                        end
                        if hasartefact(itrial) == 0
                            fprintf('Ok\n');
                        end
                    end
                    
                end
            else
                fprintf('No artefacts found! \n');
            end
            
            fprintf('Removing %d artefact(s) \n',sum(hasartefact));
            fprintf('Removing %d event(s) for which I cannot find a peak \n',sum(~haspeak));
            
            MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset    = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).offset(~hasartefact & haspeak);
            MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).clock     = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).clock(~hasartefact & haspeak);
            MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).synctime  = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(~hasartefact & haspeak);
            MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events    = MuseStruct_micro{idir}.markers.(cfg.muse.startend{imarker,1}).events(~hasartefact & haspeak);
            
            fig             = figure;
            fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap      
            h               = 1000;
            
            subplot(1,3,1);
            hold;
            for itrial = 1 : size(dat_orig_aligned.trial,2)
                if haspeak(itrial)
                    color = 'k';
                    t       = 0;
                    line([t,t],[(itrial-0.5)*h,(itrial+0.5)*h],'color','r');
                    plot(dat_orig_aligned.time{itrial},dat_orig_aligned.trial{itrial} + itrial*h,'color',color);
                else
                    color = 'r';
                end
                if hasartefact(itrial)
                    color = 'c';
                else
                    color = 'k';
                end
                
            end
            title('Original data');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            
            subplot(1,3,2);
            hold;
            for itrial = 1 : size(dat_filt_trl.trial,2)
                if haspeak(itrial)
                    color = 'k';
                    t     = dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                    line([t,t],[(itrial-0.5)*h,(itrial+0.5)*h],'color','r');
                else
                    color = 'r';
                end
                if hasartefact(itrial)
                    color = 'c';
                else
                    color = 'k';
                end
                plot(dat_filt_trl.time{itrial},dat_filt_trl.trial{itrial} + itrial*h,'color',color);
            end
            title('Peak detection');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            ax = axis;
            fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
            fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
            
            subplot(1,3,3);
            hold;
            for itrial = 1 : size(dat_sel_trl.trial,2)
                if haspeak(itrial)
                    color = 'k';
                    t       = dat_sel_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                    line([t,t],[(itrial-0.5)*h,(itrial+0.5)*h],'color','r');
                else
                    color = 'r';
                end
                if hasartefact(itrial)
                    color = 'c';
                else
                    color = 'k';
                end
                plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
            end
            title('Aligned data');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            ax = axis;
            fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
            fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                      
            % print to file
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'alignment_',cfg.name{imarker},'_',num2str(idir),'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'alignment_',cfg.name{imarker},'_',num2str(idir),'.png']),'-r600');
            close all
        end
    end
    % save data
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_micro_aligned.mat']),'MuseStruct_micro_aligned');
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_macro_aligned.mat']),'MuseStruct_macro_aligned');
end

