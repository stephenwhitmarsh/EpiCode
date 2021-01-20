function plotdata_stim(cfg,dat_stim, micromed_markers)

% select channels for plotting
i50 = 1;
i1 = 1;
clear dat_50* dat_1*

markers_50hz = micromed_markers(micromed_markers.Frequency == 50,:);
markers_1hz = micromed_markers(micromed_markers.Frequency == 1,:);

dat_50hz            = ft_appenddata([],dat_stim{find(micromed_markers.Frequency == 50)});
% dat_1hz             = ft_appenddata([],dat_stim{find(micromed_markers.Frequency == 1)});
clear dat_stim

ustimlabel = 1;
for i = 2 : height(markers_50hz)
    if ~strcmp(markers_50hz.contact1{i}(1:end-2),markers_50hz.contact1{i-1}(1:end-2))
        ustimlabel(i) = ustimlabel(i-1) + 1;
    else
        ustimlabel(i) = ustimlabel(i-1);
    end
end

ureclabel = 1;
for i = 2 : size(dat_50hz.label,1)
    if ~strcmp(dat_50hz.label{i}(1:end-2),dat_50hz.label{i-1}(1:end-2))
        ureclabel(i) = ureclabel(i-1) + 1;
    else
        ureclabel(i) = ureclabel(i-1);
    end
end

ureclabelmax        = max(ureclabel);
ustimlabelmax       = max(ustimlabel);

% filter
cfgtemp             = [];
cfgtemp.hpfilter    = 'yes';
cfgtemp.hpfreq      = 100;
dat_50hz            = ft_preprocessing(cfgtemp,dat_50hz);

% loop over recorded micro bundels
for irec = 1 : ureclabelmax

    % loop over stimulating electrodes
    for istimlabel = 1 : ustimlabelmax
       

        % loop over stimulating electode contacts (looping over all
        % electrodes in the data)
        for istimelec = find(ustimlabel == istimlabel)
            
            fig = figure; hold;
            h_sum = 0;
            h_max = max(max(abs(dat_50hz.trial{istimelec}(find(ureclabel == irec),:))));
            
            % loop over recorded micro wires
            for irecwire = find(ureclabel == irec)

                h_current = max(abs(dat_50hz.trial{istimelec}(irecwire,:))) * 1.2;
                h_current = h_max;
                h_sum = h_sum + h_current;
                
                plot(dat_50hz.time{istimelec},      dat_50hz.trial{istimelec}(irecwire,:) + h_sum,'k');  
                text(dat_50hz.time{istimelec}(1),   h_sum,dat_50hz.label{irecwire},'HorizontalAlignment','right','fontsize',8);

%                 text(dat_50hz.time{istimelec}(end), h_sum,[num2str(markers_50hz.Current(istimelec)),'mA'],'fontsize',8);
%                 text(dat_50hz.time{istimelec}(1),   h_sum,markers_50hz.Label(istimelec),'HorizontalAlignment','right','fontsize',8);
            end
            
            % write to PDF
            titlestr1 = [dat_50hz.label{irecwire}(1:end-2),': 50Hz ',cell2mat(markers_50hz.Label(istimelec,:)),' ',num2str(markers_50hz.Current(istimelec,:)),'mA'];
            temp = strsplit(markers_50hz.DataFileFolder_NL(istimelec,:),'/');
            titlestr2 = sprintf('%s @ %.0f seconds',temp{size(temp,2)},markers_50hz.Seconds_NL(istimelec));

            title(sprintf('%s\n%s',titlestr1,titlestr2),'Interpreter', 'none');
            axis tight
            fig.Renderer='Painters';
            set(fig,'PaperOrientation','portrait');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            set(gca,'yticklabel','');
            fname = fullfile(cfg.imagesavedir,['Stim_50Hz_',dat_50hz.label{irecwire}(1:end-2),'_',cell2mat(markers_50hz.Label(istimelec,:)),'_',num2str(markers_50hz.Current(istimelec,:)),'mA.pdf']);
            fprintf('\n%s\n',fname);
            print(fig, '-dpdf', fname,'-r600');
%             print(fig, '-dpdf', fullfile(cfg.imagesavedir,['Stim_50Hz_',num2str(round(rand(1)*1000)),'.pdf']),'-r300');
            close all
        end
        
    end
    
end


% 
% % plot 50 hz
% h_max = 0;
% for itrial = 1 : size(dat_50_avg,2)
%     if max(abs(dat_50_avg{itrial}.trial{1})) > h_max
%         h_max = max(abs(dat_50_avg{itrial}.trial{1}));
%     end
% end
% 
% fig = figure; hold;
% h_sum = 0;
% for itrial = 1 : size(dat_50_avg,2)
%     h_current = max(abs(dat_50_avg{itrial}.trial{1}));
%     h_current = h_max;
%     h_sum = h_sum + h_current;
%     plot(dat_50_avg{itrial}.time{1},dat_50_avg{itrial}.trial{1} + h_sum,'k');
%     text(dat_50_avg{itrial}.time{1}(end),h_sum,num2str(trialinfo_50hz.Current(itrial)),'fontsize',8);
%     text(dat_50_avg{itrial}.time{1}(1),h_sum,trialinfo_50hz.Label(itrial),'HorizontalAlignment','right','fontsize',8);
% end
% 
% % write to PDF
% title(cfg.channel{1}(1:end-2));
% axis tight
% fig.Renderer='Painters';
% set(fig,'PaperOrientation','portrait');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0 1 1]);
% set(gca,'yticklabel','');
% print(fig, '-dpdf', fullfile(cfg.imagesavedir,['trials_50Hz_',cfg.channel{1}(1:end-2),'.pdf']),'-r600');
% 
% % plot 1hz
% h_max = 0;
% for itrial = 1 : size(dat_1_avg,2)
%     if max(abs(dat_1_avg{itrial}.trial{1})) > h_max
%         h_max = max(abs(dat_1_avg{itrial}.trial{1}));
%     end
% end
% 
% fig = figure; hold;
% h_sum = 0;
% for itrial = 1 : size(dat_1_avg,2)
%     h_current = max(abs(dat_1_avg{itrial}.trial{1}));
%     h_current = h_max;
%     h_sum = h_sum + h_current;
%     plot(dat_1_avg{itrial}.time{1},dat_1_avg{itrial}.trial{1} + h_sum,'k');
%     text(dat_1_avg{itrial}.time{1}(end),h_sum,num2str(trialinfo_1hz.Current(itrial)),'fontsize',8);
%     text(dat_1_avg{itrial}.time{1}(1),h_sum,trialinfo_1hz.Label(itrial),'HorizontalAlignment','right','fontsize',8);
% end
% 
% % write to PDF
% title(cfg.channel{1}(1:end-2));
% axis tight
% fig.Renderer='Painters';
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0 1 1]);
% set(gca,'yticklabel','');
% print(fig, '-dpdf', fullfile(cfg.imagesavedir,['trials_1Hz_',cfg.channel{1}(1:end-2),'.pdf']),'-r600');
% 
% 
% 
% 
% 
% %
% %
% %
% %
% %
% %
% %
% % % calculate correlation between micro and macro
% %
% % % first rereference micro and macro to most superficial macro
% %
% % chanindx = 0;
% % for ichan = 1 : size(dat_micro.label,1)
% %     if strcmp(cfg.channel,dat_micro.label{ichan})
% %         chanindx = ichan;
% %     end
% % end
% % if chanindx == 0
% %     fprintf('Can not find channel %s! \n',cfg.channel);
% % end
% %
% % % corrs = zeros(size(dat_micro.trial,2),size(dat_macro.label,1));
% % corrs = [];
% % for itrial = 1 : size(dat_micro.trial,2)
% %     fprintf('Correlating trial %d of %d \n',itrial, size(dat_micro.trial,2));
% %     t = 1 : (size(dat_micro.trial{itrial}(chanindx,:),2) - 2); % to allow for rounding issues
% %     A = dat_micro.trial{itrial}(chanindx,t)';
% %     B = dat_macro.trial{itrial}(:,t)';
% %     Br = B(:,1:7) - B(:,2:8);
% %     Ar = A;
% % %     Br = B - B(:,8);
% % %     Ar = A - B(:,8);
% %     corrs(itrial,:) = corr(Ar, Br, 'type','pearson');
% % end
% %
% % corravg = mean(corrs,1);
% % corrsem = std(corrs,1) ./ sqrt(size(dat_macro.trial,2));
% % [~, corrmax] = max(abs(corravg));
% %
% % fname = fullfile(cfg.datasavedir,['freqanalysis_',cfg.label,'.mat']);
% % if exist(fname,'file') && cfg.force == false
% %     load(fname,'FFT_micro_trials','TFR_micro_trials','TFR_macro_trials');
% % else
% %
% %     % frequency analysis of whole trials peaks
% %     cfgtemp                 = [];
% %     cfgtemp.method          = 'mtmfft';
% %     cfgtemp.output          = 'pow';
% %     cfgtemp.taper           = 'hanning';
% %     cfgtemp.foilim          = [1,15];
% %     cfgtemp.pad             = 'nextpow2';
% %
% %     FFT_micro_trials        = ft_freqanalysis(cfgtemp,dat_micro);
% %     % FFT_macro_trials        = ft_freqanalysis(cfgtemp,dat_macro);
% %
% %     % time frequency analysis around peaks
% %     cfgtemp                 = [];
% %     cfgtemp.method          = 'mtmconvol';
% %     cfgtemp.output          = 'pow';
% %     cfgtemp.taper           = 'hanning';
% %     cfgtemp.pad             = 'nextpow2';
% %     cfgtemp.keeptrials      = 'no';
% %     cfgtemp.foi             = 10:1:300;
% % %     cfgtemp.t_ftimwin       = ones(size(cfgtemp.foi))*0.5;
% %     cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
% %
% %     cfgtemp.toi             = -cfg.prestim:cfg.slidestep:cfg.poststim;
% %     TFR_micro_trials        = ft_freqanalysis(cfgtemp,dat_micro);
% %     TFR_macro_trials        = ft_freqanalysis(cfgtemp,dat_macro);
% %
% %     save(fname,'FFT_micro_trials','TFR_micro_trials','TFR_macro_trials','-v7.3');
% %
% % end
% %
% % %%%%%%%%%%%%%%%%%%%%%
% % %%% plot overview %%%
% % %%%%%%%%%%%%%%%%%%%%%
% %
% % close all
% % fig = figure;
% %
% % subplot(4,2,1);
% % hold;
% % h = 150;
% % n = 1;
% % for itrial = cfg.representtrials
% %     i1 = find(dat_micro_chansel.time{itrial} >= -cfg.prestim,1,'first');
% %     i2 = find(dat_micro_chansel.time{itrial} >= cfg.poststim,1,'first');
% %     plot(dat_micro_chansel.time{itrial}(i1:i2),dat_micro_chansel.trial{itrial}(i1:i2) + n*h,'color','k');
% %     n = n + 1;
% % end
% % ylabel('Trials');
% % xlabel('Time (s)');
% % axis tight
% %
% % subplot(4,2,2);
% % hold;
% % h = 150;
% % n = 1;
% % for itrial = cfg.representtrials
% %     i1 = find(dat_macro_chansel.time{itrial} >= -cfg.prestim,1,'first');
% %     i2 = find(dat_macro_chansel.time{itrial} >= cfg.poststim,1,'first');
% %     plot(dat_macro_chansel.time{itrial}(i1:i2),dat_macro_chansel.trial{itrial}(i1:i2) + n*h,'color','k');
% %     n = n + 1;
% % end
% % ylabel('Trials');
% % xlabel('Time (s)');
% % axis tight
% %
% % % plot grid of trial voltages (micro)
% % subplot(4,2,3);
% % image(trialgrid_micro*255);
% % colormap hot(255)
% % set(gca, 'XTickLabel', '');
% % % title('Raw trial amplitudes over time');
% % ylabel('Trials');
% % axis tight
% % ax = axis;
% %
% % % plot grid of trial voltages (micro)
% % subplot(4,2,4);
% % image(trialgrid_macro*255);
% % colormap hot(255)
% % % title('Raw trial amplitudes over time');
% % ylabel('Trials');
% % axis tight
% % set(gca, 'XTickLabel', '');
% %
% % % plot TFR micro
% % subplot(4,2,5);
% % cfgtemp = [];
% % cfgtemp.baseline = 'yes';
% % cfgtemp.baselinetype = 'relative';
% % cfgtemp.xlim = [-cfg.prestim cfg.poststim];
% % cfgtemp.colorbar = 'no';
% % cfgtemp.title = ' ';
% % ft_singleplotTFR(cfgtemp,TFR_micro_trials);
% % % title('Time-frequency-representation');
% % ylabel('Freq');
% % set(gca, 'XTickLabel', '');
% %
% % % plot TFR macro
% % subplot(4,2,6);
% % cfgtemp = [];
% % cfgtemp.channel = corrmax;
% % cfgtemp.baseline = 'yes';
% % cfgtemp.baselinetype = 'relative';
% % cfgtemp.xlim = [-cfg.prestim cfg.poststim];
% % cfgtemp.colorbar = 'no';
% % cfgtemp.title = ' ';
% % ft_singleplotTFR(cfgtemp,TFR_macro_trials);
% % % title('Time-frequency-representation');
% % ylabel('Freq');
% % set(gca, 'XTickLabel', '');
% %
% % % plot average peak micro
% % subplot(4,2,7); hold;
% % [~,~,W,P] = findpeaks(dat_micro_rptavg.avg,dat_micro_rptavg.time,'MinPeakProminence',40,'Annotate','extents');
% % findpeaks(dat_micro_rptavg.avg,dat_micro_rptavg.time,'MinPeakProminence',40,'Annotate','extents');
% % legend('off');
% % ylabel('Microvolts');
% % [~, ci] = max(P);
% % title(sprintf('Peak analysis: W=%.0fms, P=%.0fmV',W(ci)*1000,P(ci)));
% % axis tight
% % set(gca, 'XTickLabel', '');
% %
% % % plot correlations with macro electrodes
% % subplot(4,2,8); hold;
% % boxplot(corrs,'orientation','horizontal','outliersize',1);
% %
% % % % plot FFT micro
% % % subplot(6,2,11); hold;
% % % plot(FFT_micro_trials.freq,mean(FFT_micro_trials.powspctrm,1),'k');
% % % [ymax,imax] = max(mean(FFT_micro_trials.powspctrm,1));
% % % line([FFT_micro_trials.freq(imax),FFT_micro_trials.freq(imax)],[0,ymax],'color','r','linewidth',2);
% % % txt = ['\leftarrow ', num2str(FFT_micro_trials.freq(imax),3),'Hz'];
% % % text(FFT_micro_trials.freq(imax),ymax,txt,'color','r');
% % % % title('Spectral analysis');
% % % xlabel('Hz');
% % % ylabel('Power');
% % %
% % % % histogram of durations
% % % subplot(6,2,12); hold;
% % % triallength = zeros(size(dat_micro.trial,2),1);
% % % for itrial = 1 : size(dat_micro.trial,2)
% % %     triallength(itrial) = dat_micro.time{itrial}(end) - cfg.poststim;
% % % end
% % % [histcount,edges,bin] = histcounts(triallength,0:cfg.binsize:max(triallength)+1);
% % % bar(edges(2:end)-cfg.binsize/2,histcount,1,'facecolor','k')
% % % % title('Distribution of durations');
% % % ylabel('Nr. of observations');
% % % xlabel('Duration (s)');
% %
% % %
% % % [ymax,imax] = max(histcount);
% % %
% % % [m] = median(triallength);
% % %
% % % line([m, m],[0,ymax],'color','b','linewidth',2);
% % % txt = [num2str(edges(imax)+cfg.binsize/2,3),'s'];
% % % text(edges(imax)+cfg.binsize,ymax*0.8,txt,'color','b')
% % %
% % % line([edges(imax)+cfg.binsize/2,edges(imax)+cfg.binsize/2],[0,ymax],'color','r','linewidth',2);
% % % txt = [num2str(edges(imax)+cfg.binsize/2,3),'s'];
% % % text(edges(imax)+cfg.binsize,ymax,txt,'color','r')
% % %
% % % [yleft,ileft] = find(histcount>0,1,'first');
% % % txt = ['\downarrow ', num2str(edges(ileft)+cfg.binsize/2) ,'s'];
% % % text(edges(ileft)+cfg.binsize,yleft+20,txt,'color','r')
% % %
% % % [yright,iright] = find(histcount>0,1,'last');
% % % txt = ['\downarrow ', num2str(edges(iright)+cfg.binsize/2) ,'s'];
% % % text(edges(iright)+cfg.binsize/2,yright+20,txt,'color','r')
% %
% % % print to file
% % set(fig,'PaperOrientation','landscape');
% % set(fig,'PaperUnits','normalized');
% % set(fig,'PaperPosition', [0 0 1 1]);
% % print(fig, '-dpdf', fullfile(cfg.imagesavedir,['overview_analyses_',cfg.label,'.pdf']),'-r600');
% %
% %
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% plot single event examples %%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % fig = figure;
% % n = 1;
% % for itrial = cfg.representtrials(1:10)
% %
% %     cfgtemp                 = [];
% %     cfgtemp.trials          = itrial;
% %     cfgtemp.method          = 'mtmconvol';
% %     cfgtemp.output          = 'pow';
% %     cfgtemp.taper           = 'hanning';
% %     cfgtemp.pad             = 'nextpow2';
% %     cfgtemp.foi             = 10:1:300;
% %     cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
% %     cfgtemp.toi             = -cfg.prestim:cfg.slidestep:cfg.poststim;
% %     cfgtemp.title           = ['Micro trialnr. ' num2str(itrial)];
% %     cfgtemp.channel         = cfg.channel;
% %     TFR_micro_single        = ft_freqanalysis(cfgtemp,dat_micro);
% %
% %     cfgtemp.channel         = 1;
% %     cfgtemp.title           = ['Micro trialnr. ' num2str(itrial)];
% %     TFR_macro_single        = ft_freqanalysis(cfgtemp,dat_macro);
% %
% %     subplot(size(cfg.representtrials,2)/2,2,n); hold;
% %     cfgtemp = [];
% %     cfgtemp.channel = cfg.channel;
% %     cfgtemp.baseline = 'yes';
% %     cfgtemp.baselinetype = 'relative';
% %     ft_singleplotTFR(cfgtemp,TFR_micro_single);
% %
% %     % add scaled LFP line
% %     hold;
% %     ax = axis;
% %     scaled = (dat_micro.trial{itrial}(chanindx,:) + abs(min(dat_micro.trial{itrial}(chanindx,:))));
% %     scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
% %     i1 = find(dat_micro.time{itrial} > TFR_micro_single.time(1),1,'first');
% %     i2 = find(dat_micro.time{itrial} < TFR_micro_single.time(end),1,'last');
% %     plot(dat_micro.time{itrial}(i1:i2),scaled(i1:i2),'linewidth',1,'color',[1 1 1]);
% %
% %     axis tight
% %
% %     subplot(size(cfg.representtrials,2)/2,2,n+1);
% %     cfgtemp = [];
% %     cfgtemp.channel = 'all';
% %     cfgtemp.baseline = 'yes';
% %     cfgtemp.baselinetype = 'relative';
% %     ft_singleplotTFR(cfgtemp,TFR_macro_single);
% %
% %     % add scaled LFP line
% %     hold;
% %     ax = axis;
% %     scaled = (dat_macro.trial{itrial}(chanindx,:) + abs(min(dat_macro.trial{itrial}(chanindx,:))));
% %     scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
% %     i1 = find(dat_macro.time{itrial} > TFR_macro_single.time(1),1,'first');
% %     i2 = find(dat_macro.time{itrial} < TFR_macro_single.time(end),1,'last');
% %     plot(dat_macro.time{itrial}(i1:i2),scaled(i1:i2),'linewidth',1,'color',[1 1 1]);
% %
% %     n = n + 2;
% %
% %     colormap hot
% % end
% %
% % % print to file
% % set(fig,'PaperOrientation','landscape');
% % set(fig,'PaperUnits','normalized');
% % set(fig,'PaperPosition', [0 0 1 1]);
% % print(fig, '-dpdf', fullfile(cfg.imagesavedir,['single_trials_',cfg.label,'.pdf']),'-r600');
% %
% %
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% calculate correlation LFP x power %%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % fname = fullfile(cfg.datasavedir,['freqanalysis_singletrial_micro_',cfg.label,'.mat']);
% % if exist(fname,'file') && cfg.force == false
% %     fprintf('loading %s, which can take a while \n',fname);
% %     load(fname,'TFR_micro_singletrial','lfp_corr');
% % else
% %
% %     % time frequency analysis around peaks
% %     cfgtemp                 = [];
% %     cfgtemp.method          = 'mtmconvol';
% %     cfgtemp.output          = 'pow';
% %     cfgtemp.taper           = 'hanning';
% %     cfgtemp.pad             = 'nextpow2';
% %     cfgtemp.keeptrials      = 'yes';
% %     cfgtemp.foi             = 10:1:300;
% %     cfgtemp.channel         = chanindx;
% %     cfgtemp.t_ftimwin       = ones(size(cfgtemp.foi))*0.5;
% % %     cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
% %     cfgtemp.toi             = -cfg.prestim:cfg.slidestep:cfg.poststim;
% %     TFR_micro_singletrial   = ft_freqanalysis(cfgtemp,dat_micro);
% %
%     clear lfp_resampled lfp_corr
%     for itrial = 1 : size(dat_micro.trial,2)
%         fprintf('Trial %d of %d \n',itrial, size(dat_micro.trial,2));
%
%         cfgtemp = [];
%         cfgtemp.time{1}     = TFR_micro_singletrial.time;
%         cfgtemp.trials      = itrial;
%         cfgtemp.method      = 'pchip';
%         temp                = ft_resampledata(cfgtemp,dat_micro);
%         lfp_corr(itrial,:)  = corr(permute(TFR_micro_singletrial.powspctrm(itrial,1,:,:),[4,3,2,1]),temp.trial{1}(chanindx,:)','rows','pairwise');
%
%     end
%     save(fname,'TFR_micro_singletrial','lfp_corr','-v7.3');
% end
%
% lfp_corr_avg = mean(lfp_corr,1);
% lfp_corr_std = std(lfp_corr,1);
%
% % plot average peak micro
% fig = figure;
%
% subplot(3,1,1);
% [PKS,LOCS,W,P] = findpeaks(lfp_corr_avg,TFR_micro_singletrial.freq,'MinPeakProminence',0.05,'Annotate','extents');
% findpeaks(lfp_corr_avg,TFR_micro_singletrial.freq,'MinPeakProminence',0.05,'Annotate','extents');
% legend('off');
% ylabel('Correlation rho');
% [~, ci] = max(P);
% title(sprintf('%.1fHz ',LOCS(ci)));
% axis tight
% set(gca, 'XTickLabel', '');
%
% subplot(3,1,2); hold;
% patch([TFR_micro_singletrial.freq,TFR_micro_singletrial.freq(end:-1:1)],[lfp_corr_avg+lfp_corr_std,lfp_corr_avg(end:-1:1)-lfp_corr_std(end:-1:1)],[0.8,0.8,0.8],'edgecolor','none');
% plot(TFR_macro_single.freq,mean(lfp_corr,1),'color',[0 0 0]);
% axis tight
% ylabel('Power');
% ylabel('Correlation rho');
% set(gca, 'XTickLabel', '');
%
% FFT_avg = squeeze(nanmean(nanmean(TFR_micro_singletrial.powspctrm,4),1))';
% FFT_std = squeeze(std(nanmean(TFR_micro_singletrial.powspctrm,4),1))';
%
% subplot(3,1,3); hold;
%
% % dont know why it fails one way or another
% try
%     patch([TFR_micro_singletrial.freq,TFR_micro_singletrial.freq(end:-1:1)],[FFT_avg+FFT_std,FFT_avg(end:-1:1)-FFT_std(end:-1:1)],[0.8,0.8,0.8],'edgecolor','none');
% catch
% end
% try
%     patch([TFR_micro_singletrial.freq,TFR_micro_singletrial.freq(end:-1:1)],[(FFT_avg(:,chanindx)+FFT_std(:,chanindx))',(FFT_avg(end:-1:1,chanindx)-FFT_std(end:-1:1,chanindx))'],[0.8,0.8,0.8],'edgecolor','none');
% catch
% end
%
% plot(TFR_macro_single.freq,FFT_avg,'color',[0 0 0]);
% set(gca, 'YScale', 'log');
% ylabel('Power');
% axis tight
% xlabel('Frequency (Hz)');
%
% % print to file
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0 1 1]);
% print(fig, '-dpdf', fullfile(cfg.imagesavedir,['spectral_correlations_',cfg.label,'.pdf']),'-r600');
%
% close all
