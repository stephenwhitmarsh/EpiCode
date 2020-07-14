function [MuseStruct] = alignMuseMarkers_EMGthresh(cfg,MuseStruct, force)
ft_debug off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revoir texte, par rapport à alignMuseMarkers
% remove trials of hasartefacts
% no contraction during bl
% sd of raw. find on filt
% Paul Baudin, helped by Stephen Whitmarsh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_EMGaligned.mat']);
%function for later detection 
zci                     = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %find time cross zero


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


% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)
    
    
    for imarker = 1 : size(cfg.alignEMG.name,2)%markerlist
        if ~strcmp(cfg.alignEMG.channel{imarker},'no')
            
            % find data directories that have the required event
            markerindx = [];
            for idir = 1 : size(MuseStruct{ipart},2)
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.alignEMG.EMGmarker{imarker})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}),'synctime')
                            if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).synctime)
                                fprintf('Found %d event timings for %s \n',size(MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).clock,2),cfg.alignEMG.EMGmarker{imarker});
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
                        if strfind(filelist(ifile).name,cfg.alignEMG.channel{imarker})
                            fprintf('Found channel with pattern "%s" in %s\n',cfg.alignEMG.channel{imarker},filelist(ifile).name);
                            dataset = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(ifile).name);
                            channelnr = channelnr + 1;
                        end
                    end
                    if channelnr == 0
                        fprintf('Did not find of channel pattern %s!\n',cfg.alignEMG.channel{imarker});
                    end
                    if channelnr > 1
                        fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignEMG.channel{imarker});
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
                    cfgtemp.channel           = {cfg.alignEMG.channel{imarker}, cfg.EMG.refchannel};%load the emg to align and the ref
                    cfgtemp.dataset           = dataset;
                    cfgtemp.reref             = cfg.EMG.reref;
                    cfgtemp.rerefmethod       = cfg.EMG.rerefmethod;
                    cfgtemp.refchannel        = cfg.EMG.refchannel;
                    data_EMG                  = ft_preprocessing(cfgtemp);
                else
                    cfgtemp                       = [];
                    cfgtemp.channel               = cfg.alignEMG.channel{imarker};%load only the emg associated with eeg marker
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
                        cfgtemp.hpfreq            = cfg.EMG.hpfreq;
                        data_EMG                  = ft_preprocessing(cfgtemp, data_EMG);
                    end
                end
                
                %remove reference channel if any
                cfgtemp         = [];
                cfgtemp.channel = cfg.alignEMG.channel{imarker};
                dat_sel         = ft_selectdata(cfgtemp,data_EMG);
                if size(dat_sel.trial{1},1) == 0
                    fprintf('Did not find of channel pattern %s!\n',cfg.alignEMG.channel{imarker});
                end
                if size(dat_sel.trial{1},1) > 1
                    fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignEMG.channel{imarker});
                end
                
                % segment in trials for alignment
                startsample         = round(MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).synctime*dat_sel.fsample + cfg.alignEMG.toiplot{imarker}(1)*dat_sel.fsample + hdr.nSamplesPre)';
                endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).synctime*dat_sel.fsample + cfg.alignEMG.toiplot{imarker}(2)*dat_sel.fsample + hdr.nSamplesPre)';
                offset              = round(ones(size(startsample))*+cfg.alignEMG.toiplot{imarker}(1)*dat_sel.fsample);
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
                
                %envelope of emg
                for itrial = 1 : size(dat_filt_trl.trial, 2)
                    [dat_filt_trl.trial{itrial}, ~] = envelope(dat_filt_trl.trial{itrial}(1,:),cfg.EMG.envparam,cfg.EMG.envmethod);
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
                    t1_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.alignEMG.toibaseline{imarker}(1),1,'first');
                    t2_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.alignEMG.toibaseline{imarker}(2),1,'last');
                    t1_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.alignEMG.toiactive{imarker}(1),1,'first');
                    t2_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.alignEMG.toiactive{imarker}(2),1,'last');
                end
                
                 for itrial = 1 : size(dat_filt_trl.trial,2)
                     mean_bl(itrial)         = nanmean(dat_filt_trl.trial{itrial}(t1_bl_indx{itrial}:t2_bl_indx{itrial}));
                     max_ac(itrial)                  = max(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}:t2_ac_indx{itrial}));
                     if max_ac(itrial)<mean_bl(itrial) * 2; haspeak(itrial) = false; end
                 end
                 med_all_bl = nanmedian(mean_bl);
                 
                 for itrial = 1 : size(dat_filt_trl.trial,2)
                    if haspeak(itrial)
                        %Paul find last cross threshold
                         
                        thresh(itrial)          = (max_ac(itrial) - med_all_bl) * 0.1;
%                         peak_loc(itrial)        = locs_ac_sel{itrial}(ip(itrial))+t1_ac_indx(itrial)-1;
                        %find cross thresh between begin of bl and detected peak
                        index_tresh{itrial}     = zci(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}:t2_ac_indx{itrial}) - med_all_bl - thresh(itrial)) + t1_ac_indx{itrial} - 1;
%                         plot(dat_filt_trl.time{itrial}, dat_filt_trl.trial{itrial});
%                         plot([dat_filt_trl.time{itrial}(index_tresh{itrial}), dat_filt_trl.time{itrial}(index_tresh{itrial})], [-100 100]);
%                         plot([dat_filt_trl.time{itrial}(peak_loc(itrial)), dat_filt_trl.time{itrial}(peak_loc(itrial))], [400 300]);
%                         plot(dat_filt_trl.time{itrial}(t1_ac_indx(itrial):peak_loc(itrial)), dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):peak_loc(itrial)) - mean_bl(itrial) - thresh(itrial));
                        if ~isempty(index_tresh{itrial} )
                            % align time axis to relevant (peak) index 
                            index_tresh{itrial}     = index_tresh{itrial}(end-1);
                            timeshift                       = dat_filt_trl.time{itrial}(index_tresh{itrial});
                            dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                            dat_filt_aligned.time{itrial}   = dat_filt_trl.time{itrial} - timeshift; %for plot

                            if abs(timeshift) > 0.3
                                hasartefact(itrial) = true;
                            end 

                            if ~hasartefact
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial) + timeshift;
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial)          = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial) + seconds(timeshift);
                            end
                        else
                            haspeak(itrial) = false;
                        end
                    end
                end
                
                %remove supplemental trials in case some are ignored
                %because of hasartefact or ~haspeak
                MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).synctime = MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).synctime(1:keeptrial);
                MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).clock = MuseStruct{ipart}{idir}.markers.(cfg.alignEMG.EMGmarker{imarker}).clock(1:keeptrial);
                
                %% Plot alignement
                
                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    h_temp(itrial) = max(dat_sel_trl.trial{itrial}(t1_ac_indx{itrial}: t2_ac_indx{itrial})); %amplitude of peak vs baseline
                end
                h               = mean(h_temp);
           

                subplot(2,2,1);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t     = dat_filt_trl.time{itrial}(index_tresh{itrial});
                        line_position = dat_filt_trl.trial{itrial}(index_tresh{itrial});
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
%                     else
%                         color = 'k';
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
                
                subplot(2,2,2);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = 0;
                        line_position = dat_filt_trl.trial{itrial}(index_tresh{itrial});
                        if isempty(line_position); t=[]; end
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r'); 
                        
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
%                     else
%                         color = 'k';
                    end
                   plot(dat_filt_aligned.time{itrial},dat_filt_aligned.trial{itrial}+ itrial*h,'color',color); 
                    
                end
                title('Alignment');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                
               
                subplot(2,2,3);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = dat_sel_trl.time{itrial}(index_tresh{itrial});
                        line_position = dat_sel_trl.trial{itrial}(index_tresh{itrial});
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
%                     else
%                         color = 'k';
                    end
                    plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
                end
                title('Original data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                ax = axis;
                %fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                %fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                
                subplot(2,2,4);
                hold;
                for itrial = 1 : size(dat_sel_aligned.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = 0;
                        line_position = dat_sel_trl.trial{itrial}(index_tresh{itrial});
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                        if isempty(line_position); t=[]; end
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
%                     else
%                         color = 'k';
                    end
                    plot(dat_sel_aligned.time{itrial},dat_sel_aligned.trial{itrial} + itrial*h,'color',color);
          
                end
                title('Aligned data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                
                %% print to file
                
                % check if images directory exists, if not create
                if ~isfolder(cfg.imagesavedir)
                    ft_notice('creating directory %s', cfg.imagesavedir);
                    mkdir(cfg.imagesavedir);
                end
                
                % check if aligment subdirectory exists, if not create
                if ~isfolder(fullfile(cfg.imagesavedir,'alignment_emg_tresh'))
                    ft_notice('creating directory %s', fullfile(cfg.imagesavedir,'alignment_emg_tresh'));
                    mkdir(fullfile(cfg.imagesavedir,'alignment_emg_tresh'));
                end
                
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'alignment_emg_tresh',[cfg.prefix,'p',num2str(ipart),'_',cfg.alignEMG.name{imarker},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.pdf']));
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'alignment_emg_tresh',[cfg.prefix,'p',num2str(ipart),'_',cfg.alignEMG.name{imarker},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.png']),'-r600');
                
                close all
                
            end % idir
        end
    end % imarker
end % ipart

%% save data
save(fname,'MuseStruct');
end

