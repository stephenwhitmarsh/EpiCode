function [MuseStruct] = alignMuseMarkers_EMG(cfg,MuseStruct, force)
ft_debug off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revoir texte, par rapport à alignMuseMarkers
% remove trials of hasartefacts
% no contraction during bl
% sd of raw. find on filt
% Paul Baudin, helped by Stephen Whitmarsh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg,'alignEMG')
    fprintf('cfg.alignEMG does not exist, nothing done\n');
    return
end

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_EMGaligned.mat']);

if exist(fname,'file') && force == false
    fprintf('***********************************\n');
    fprintf('** Loading results alignment EMG **\n');
    fprintf('***********************************\n\n');
    load(fname,'MuseStruct');
    return
elseif exist(fname,'file') && force == true
    fprintf('*************************************\n');
    fprintf('** Forced redoing of alignment EMG **\n');
    fprintf('*************************************\n\n');
else
    fprintf('******************\n');
    fprintf('** Aligning EMG **\n');
    fprintf('******************\n\n');
end

%get format to adapt script for each format
%specificities :
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% select those markers to align
markerlist = [];
for i = 1 : size(cfg.name,2)
    if ismember(cfg.name{i}, cfg.alignEMG.name)
        markerlist = [markerlist, i];
    end
end

% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)
    
    imarkeralign = 0;
    for imarker = markerlist
        imarkeralign = imarkeralign+1;
        
        if ~strcmp(cfg.alignEMG.channel{imarkeralign},'no')
            
            % find data directories that have the required event
            markerindx = [];
            for idir = 1 : size(MuseStruct{ipart},2)
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                            if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                                fprintf('Found %d event timings for %s \n',size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock,2),cfg.muse.startend{imarker,1});
                                markerindx = [markerindx idir];
                            end
                        end
                    end
                end
            end
            if isempty(markerindx)
                fprintf('Did not find any events for marker %s\n',cfg.name{imarker});
            end
            
            % loop over all data directories that have the marker
            for idir = markerindx
                
                %% find datafilename corresponding to requested electrode
                if isNeuralynx
                    filelist = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
                    channelnr = 0;
                    for ifile = 1 : size(filelist,1)
                        if strfind(filelist(ifile).name,cfg.alignEMG.channel{imarkeralign})
                            fprintf('Found channel with pattern "%s" in %s\n',cfg.alignEMG.channel{imarkeralign},filelist(ifile).name);
                            dataset = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(ifile).name);
                            channelnr = channelnr + 1;
                        end
                    end
                    if channelnr == 0
                        fprintf('Did not find of channel pattern %s!\n',cfg.alignEMG.channel{imarkeralign});
                    end
                    if channelnr > 1
                        fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignEMG.channel{imarkeralign});
                    end
                    
                elseif isMicromed
                    dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
                    
                elseif isBrainvision
                    dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
                    
                end
                
                %load header
                hdr = ft_read_header(dataset);
                
                %% load EMG data
                %not done yet : for Neuralynx data, if need of rereferencing,
                %find the path of the channel to reref
                
                loadrefemg = false;
                if isfield(cfg.EMG, 'reref')
                    if strcmp(cfg.EMG.reref, 'yes')
                        loadrefemg = true;
                    end
                end
                
                if loadrefemg
                    cfgtemp                   = [];
                    cfgtemp.channel           = {cfg.alignEMG.channel{imarkeralign}, cfg.EMG.refchannel};%load the emg to align and the ref
                    cfgtemp.dataset           = dataset;
                    cfgtemp.reref             = cfg.EMG.reref;
                    cfgtemp.rerefmethod       = cfg.EMG.rerefmethod;
                    cfgtemp.refchannel        = cfg.EMG.refchannel;
                    data_EMG                  = ft_preprocessing(cfgtemp);
                else
                    cfgtemp                       = [];
                    cfgtemp.channel               = cfg.alignEMG.channel{imarkeralign};%load only the emg associated with eeg marker
                    cfgtemp.dataset               = dataset;
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
                        cfgtemp.hpfreq            = cfg.alignEMG.hpfreq;
                        data_EMG                  = ft_preprocessing(cfgtemp, data_EMG);
                    end
                end
                
                %remove reference channel if any
                cfgtemp         = [];
                cfgtemp.channel = cfg.alignEMG.channel{imarkeralign};
                dat_sel         = ft_selectdata(cfgtemp,data_EMG);
                if size(dat_sel.trial{1},1) == 0
                    fprintf('Did not find of channel pattern %s!\n',cfg.alignEMG.channel{imarkeralign});
                end
                if size(dat_sel.trial{1},1) > 1
                    fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignEMG.channel{imarkeralign});
                end
                
                % segment in trials for alignment
                startsample         = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime*dat_sel.fsample + cfg.alignEMG.toiplot{imarkeralign}(1)*dat_sel.fsample + hdr.nSamplesPre)';
                endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime*dat_sel.fsample + cfg.alignEMG.toiplot{imarkeralign}(2)*dat_sel.fsample + hdr.nSamplesPre)';
                offset              = round(ones(size(startsample))*+cfg.alignEMG.toiplot{imarkeralign}(1)*dat_sel.fsample);
                cfgtemp             = [];
                cfgtemp.trl         = [startsample, endsample, offset];
                cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                dat_sel_trl         = ft_redefinetrial(cfgtemp,dat_sel);
                clear dat_sel
                
                %% filter EMG data for detection
                
                % rectify EMG
                cfgtemp = [];
                cfgtemp.rectify = 'yes';
                dat_filt_trl = ft_preprocessing(cfgtemp, dat_sel_trl);
                
                % compute integral
                for itrial = 1: size(dat_filt_trl.trial,2)
                    dat_filt_trl.trial{itrial} = cumtrapz(dat_filt_trl.trial{itrial});
                end
                
                %remove linear trend
                cfgtemp = [];
                cfgtemp.detrend = 'yes';
                dat_filt_trl = ft_preprocessing(cfgtemp,dat_filt_trl);
                
                %flip for peak detection
                for itrial = 1: size(dat_filt_trl.trial,2)
                    dat_filt_trl.trial{itrial} =-(dat_filt_trl.trial{itrial});
                end
                
                %% detection and alignment
                
                %initialize variables for alignment
                dat_sel_aligned     = dat_sel_trl;
                dat_filt_aligned    = dat_filt_trl; %for plot
                haspeak         = true(1,size(dat_filt_trl.trial,2));
                hasartefact     = false(1,size(dat_filt_trl.trial,2));
                keeptrial = 0;
                
                % time indexes for baseline and active period
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    t1_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.alignEMG.toibaseline{imarkeralign}(1),1,'first');
                    t2_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.alignEMG.toibaseline{imarkeralign}(2),1,'last');
                    t1_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.alignEMG.toiactive{imarkeralign}(1),1,'first');
                    t2_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.alignEMG.toiactive{imarkeralign}(2),1,'last');
                end
                
                %find max peak of dat_filt
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    
                    % find peaks in active period
                    [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}:t2_ac_indx{itrial}));
                    
                    % find peaks in baseline period
                    [peaks_bl{itrial},~]               = findpeaks(dat_filt_trl.trial{itrial}(t1_bl_indx{itrial}:t2_bl_indx{itrial}));
                    if isempty(peaks_bl{itrial})%if no peak, take avg bl value
                        peaks_bl{itrial} = nanmean(dat_filt_trl.trial{itrial}(t1_bl_indx{itrial}:t2_bl_indx{itrial}));
                    end
                    
                    % select maximum peak in baseline
                    max_peaks_bl{itrial}           = max(peaks_bl{itrial});
                end
                
                %select peak
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    %remove peaks lower than baseline
                    peaks_sel               = peaks_ac{itrial} > max(peaks_bl{itrial});
                    peaks_ac_sel{itrial}    = peaks_ac{itrial}(peaks_sel);
                    locs_ac_sel{itrial}     = locs_ac{itrial}(peaks_sel);
                    
                    %choose maximum of selected peak
                    if ~isempty(peaks_ac_sel{itrial})
                        [peaks_ac_sel{itrial}, loc_idx]    = max(peaks_ac_sel{itrial});
                        locs_ac_sel{itrial}                = locs_ac_sel{itrial}(loc_idx);
                        timeshift               = dat_filt_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    else
                        locs_ac_sel{itrial}     = [];
                        haspeak(itrial) = false;
                        fprintf('Could not find peak in trial number %d \n',itrial);
                        timeshift = 0;
                    end
                    
                    % align data
                    dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                    dat_filt_aligned.time{itrial}   = dat_filt_trl.time{itrial} - timeshift; %for plot
                    
                    % if no detection in a certain window, set as artefact
                    if abs(timeshift) > cfg.alignEMG.maxtimeshift
                        hasartefact(itrial) = true;
                    end
                    
                    %keep only trials with good detection
                    if ~hasartefact(itrial) && haspeak(itrial)
                        keeptrial = keeptrial +1;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial) + timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial)          = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial) + seconds(timeshift);
                    end
                end
                
                %remove supplemental trials in case some are ignored
                %because of hasartefact or ~haspeak
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(1:keeptrial);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(1:keeptrial);
                
                %% Plot alignement
                
                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    h_temp(itrial) = max(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}: t2_ac_indx{itrial})); %amplitude of peak vs baseline
                end
                h               = mean(h_temp)/2;
                
                
                subplot(2,2,1);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    color = 'k';
                    t     = dat_filt_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line_position =   dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    if hasartefact(itrial)
                        color = 'c';
                    end
                    if ~haspeak(itrial)
                        color = 'r';
                    end
                    plot(dat_filt_trl.time{itrial},dat_filt_trl.trial{itrial} + itrial*h,'color',color);
                    xlim([-0.5 0.5])
                end
                title('Detection : integral of trial + detrend + flip (find max peak)');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.alignEMG.toiplot{imarkeralign});
                ax = axis;
                fill([cfg.alignEMG.toiactive{imarkeralign}(1), cfg.alignEMG.toiactive{imarkeralign}(2), cfg.alignEMG.toiactive{imarkeralign}(2), cfg.alignEMG.toiactive{imarkeralign}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.alignEMG.toibaseline{imarkeralign}(1), cfg.alignEMG.toibaseline{imarkeralign}(2), cfg.alignEMG.toibaseline{imarkeralign}(2), cfg.alignEMG.toibaseline{imarkeralign}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 0 1],'edgecolor','none','facealpha',0.1);
                
                subplot(2,2,2);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    color = 'k';
                    t       = 0;
                    line_position = dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    if isempty(line_position); t=[]; end
                    if ~hasartefact(itrial) && haspeak(itrial)
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                        plot(dat_filt_aligned.time{itrial},dat_filt_aligned.trial{itrial}+ itrial*h,'color',color);
                    end
                end
                title('Alignment');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.alignEMG.toiplot{imarkeralign});
                ax = axis;
                fill([cfg.alignEMG.toiactive{imarkeralign}(1), cfg.alignEMG.toiactive{imarkeralign}(2), cfg.alignEMG.toiactive{imarkeralign}(2), cfg.alignEMG.toiactive{imarkeralign}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.alignEMG.toibaseline{imarkeralign}(1), cfg.alignEMG.toibaseline{imarkeralign}(2), cfg.alignEMG.toibaseline{imarkeralign}(2), cfg.alignEMG.toibaseline{imarkeralign}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 0 1],'edgecolor','none','facealpha',0.1);
                
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    h_temp(itrial) = max(dat_sel_trl.trial{itrial}(t1_ac_indx{itrial}: t2_ac_indx{itrial})); %amplitude of peak vs baseline
                end
                h               = mean(h_temp);
                
                subplot(2,2,3);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    color = 'k';
                    t       = dat_sel_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line_position = dat_sel_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    
                    if hasartefact(itrial)
                        color = 'c';
                    end
                    if ~haspeak(itrial)
                        color = 'r';
                    end
                    plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
                end
                title('Original data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.alignEMG.toiplot{imarkeralign});
                
                subplot(2,2,4);
                hold;
                for itrial = 1 : size(dat_sel_aligned.trial,2)
                    color = 'k';
                    t       = 0;
                    line_position = dat_sel_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    if isempty(line_position); t=[]; end
                    if ~hasartefact(itrial) && haspeak(itrial)
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                        plot(dat_sel_aligned.time{itrial},dat_sel_aligned.trial{itrial} + itrial*h,'color',color);
                    end
                    
                end
                title('Aligned data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.alignEMG.toiplot{imarkeralign});
                
                %% print to file
                
                % check if images directory exists, if not create
                if ~isfolder(cfg.imagesavedir)
                    ft_notice('creating directory %s', cfg.imagesavedir);
                    mkdir(cfg.imagesavedir);
                end
                
                % check if aligment subdirectory exists, if not create
                if ~isfolder(fullfile(cfg.imagesavedir,'alignment_emg'))
                    ft_notice('creating directory %s', fullfile(cfg.imagesavedir,'alignment_emg'));
                    mkdir(fullfile(cfg.imagesavedir,'alignment_emg'));
                end
                
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'alignment_emg',[cfg.prefix,'p',num2str(ipart),'_',cfg.alignEMG.name{imarkeralign},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.pdf']));
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'alignment_emg',[cfg.prefix,'p',num2str(ipart),'_',cfg.alignEMG.name{imarkeralign},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.png']),'-r600');
                
                close all
                
            end % idir
        end
    end % imarker
end % ipart

%% save data
save(fname,'MuseStruct');
end

