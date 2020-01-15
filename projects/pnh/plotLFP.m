function [FFT_micro_trials,TFR_micro_trials,TFR_macro_trials,stat_TFR_micro, corrs] = plotLFP(cfg,dat_micro,dat_macro,force)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'_TFR_LFP.mat']);
if exist(fname,'file') && force == false
    fprintf('****************************************\n');
    fprintf('*** Loading precomputed TFR-LFP data ***\n');
    fprintf('****************************************\n\n');
    
    load(fname,'FFT_micro_trials','TFR_micro_trials','TFR_macro_trials','stat_TFR_micro', 'corrs');
else
    
    for imarker = 1 : size(cfg.name,2)
        
        % find datafilename corresponding to channel
        channelindx = [];
        for ilabel = 1 : size(dat_micro{imarker}.label,1)
            if strfind(dat_micro{imarker}.label{ilabel},cfg.align.channel{imarker})
                channelindx = ilabel;
            end
        end
        
        % select channel for plotting
        cfgtemp                 = [];
        cfgtemp.channel         = channelindx;
        dat_micro_chansel       = ft_selectdata(cfgtemp,dat_micro{imarker});
        cfgtemp.channel         = 1; % closes to micro
        dat_macro_chansel       = ft_selectdata(cfgtemp,dat_macro{imarker});
        
        % average over trials for plotting
        cfgtemp                 = [];
        cfgtemp.vartrllength    = 2;
        dat_micro_rptavg        = ft_timelockanalysis(cfgtemp,dat_micro_chansel);
        cfgtemp.channel         = 1;
        dat_macro_rptavg        = ft_timelockanalysis(cfgtemp,dat_macro_chansel);
        
        % select time interval
        cfgtemp                 = [];
        cfgtemp.latency         = [cfg.epoch.toi{imarker}(1) cfg.epoch.toi{imarker}(2)];
        dat_micro_rptavg        = ft_selectdata(cfgtemp,dat_micro_rptavg);
        dat_macro_rptavg        = ft_selectdata(cfgtemp,dat_macro_rptavg);
        
        % map all trials onto time-independant grid (average over channels)
        maxsamples              = -cfg.epoch.toi{imarker}(1)*dat_micro_chansel.fsample + cfg.epoch.toi{imarker}(2)*dat_micro_chansel.fsample;
        trialgrid_micro         = nan(size(dat_micro_chansel.trial,2),maxsamples);
        trialgrid_macro         = nan(size(dat_macro_chansel.trial,2),maxsamples);
        
        for itrial = 1 : size(dat_micro{imarker}.trial,2)
            
            i1 = find(dat_micro_chansel.time{itrial} >= cfg.epoch.toi{imarker}(1),1,'first');
            i2 = find(dat_micro_chansel.time{itrial} >= cfg.epoch.toi{imarker}(2),1,'first');
            
            n = (dat_micro_chansel.trial{itrial}(i1:i2) - mean(dat_micro_chansel.trial{itrial}(i1:i2))) / max(abs(dat_micro_chansel.trial{itrial}(i1:i2)));
            %     n = dat_micro_chansel.trial{itrial}(i1:i2);
            
            trialgrid_micro(itrial,1:size(n,2)) = n;
            
            n = (dat_macro_chansel.trial{itrial}(i1:i2) - mean(dat_macro_chansel.trial{itrial}(i1:i2))) / max(abs(dat_macro_chansel.trial{itrial}(i1:i2)));
            %     n = dat_macro_chansel.trial{itrial}(i1:i2);
            
            trialgrid_macro(itrial,1:size(n,2)) = n;
        end
        
        
        %% calculate correlation between first rereference micro and macro to most superficial macro
        
        % chanindx = 0;
        % for ichan = 1 : size(dat_micro{imarker}.label,1)
        %     if strcmp(cfg.channel,dat_micro{imarker}.label{ichan})
        %         chanindx = ichan;
        %     end
        % end
        % if chanindx == 0
        %     fprintf('Can not find channel %s! \n',cfg.channel);
        % end
        
        % corrs = zeros(size(dat_micro.trial,2),size(dat_macro.label,1));
        
        corrs = [];
        for itrial = 1 : size(dat_micro{imarker}.trial,2)
            fprintf('Correlating trial %d of %d \n',itrial, size(dat_micro{imarker}.trial,2));
            t = 1 : (size(dat_micro{imarker}.trial{itrial}(channelindx,:),2) - 2); % to allow for rounding issues
            A = dat_micro{imarker}.trial{itrial}(channelindx,t)';
            B = dat_macro{imarker}.trial{itrial}(:,t)';
            Br = B(:,1:7) - B(:,2:8);
            Ar = A;
            %     Br = B - B(:,8);
            %     Ar = A - B(:,8);
            [corrs{imarker}.c(itrial,:), corrs{imarker}.ptrial(itrial,:)] = corr(Ar, Br, 'type','pearson');
        end
        
        corrs{imarker}.avg = mean(corrs{imarker}.c,1);
        corrs{imarker}.sem = std(corrs{imarker}.c,1) ./ sqrt(size(dat_macro{imarker}.trial,2));
        [~, corrs{imarker}.max] = max(abs(corrs{imarker}.avg ));
        [corrs{imarker}.h,corrs{imarker}.p,corrs{imarker}.ci,corrs{imarker}.stats] = ttest(corrs{imarker}.c);
        
        % time frequency analysis around peaks
        cfgtemp                         = [];
        cfgtemp.channel                 = channelindx;
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'hanning';
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.foi                     = 10:1:200;
        cfgtemp.t_ftimwin               = 7./cfgtemp.foi;
        cfgtemp.toi                     = cfg.epoch.toi{imarker}(1):cfg.LFP.slidestep:cfg.epoch.toi{imarker}(2);
        TFR_micro_trials{imarker}       = ft_freqanalysis(cfgtemp,dat_micro{imarker});
        
        cfgtemp.channel                 = corrs{imarker}.max;
        cfgtemp.keeptrials              = 'no';
        TFR_macro_trials{imarker}       = ft_freqanalysis(cfgtemp,dat_macro{imarker});
        
        %% cluster stats on TFR
        cfgtemp                         = [];
        cfgtemp.baseline                = cfg.LFP.baselinewindow{imarker};
        cfgtemp.baselinetype            = 'relchange';
        TFR_micro_trials{imarker}       = ft_freqbaseline(cfgtemp,TFR_micro_trials{imarker});
        
        cfgtemp                         = [];
        cfgtemp.baseline                = cfg.LFP.baselinewindow{imarker};
        cfgtemp.baselinetype            = 'relchange';
        TFR_macro_trials{imarker}       = ft_freqbaseline(cfgtemp,TFR_macro_trials{imarker});
        
        cfgtemp                         = [];
        cfgtemp.latency                 = [cfg.LFP.baselinewindow{imarker}(2) cfg.epoch.toi{imarker}(2)];
        cfgtemp.frequency               = 'all';
        cfgtemp.method                  = 'montecarlo';
        cfgtemp.statistic               = 'ft_statfun_indepsamplesT';
        cfgtemp.correctm                = 'cluster';
        cfgtemp.clusteralpha            = 0.01;
        cfgtemp.clusterstatistic        = 'maxsum';
        cfgtemp.clustertail             = 0;
        cfgtemp.alpha                   = 0.025;
        cfgtemp.numrandomization        = 500;
        cfgtemp.ivar                    = 1;
        cfgtemp.design                  = [ones(1,size(TFR_micro_trials{imarker}.trialinfo,1)) ones(1,size(TFR_micro_trials{imarker}.trialinfo,1))*2];
        
        dummy                           = TFR_micro_trials{imarker};
        dummy.powspctrm                 = zeros(size(dummy.powspctrm));
        
        stat_TFR_micro{imarker}         = ft_freqstatistics(cfgtemp,TFR_micro_trials{imarker}, dummy);
        stat_TFR_micro{imarker}.corrs   = corrs{imarker};
        
        posclust = zeros(size(stat_TFR_micro{imarker}.mask));
        for ipos = 1 : size(stat_TFR_micro{imarker}.posclusters,2)
            if stat_TFR_micro{imarker}.posclusters(ipos).prob < 0.025
                posclust(stat_TFR_micro{imarker}.posclusterslabelmat == ipos) = 1;
            end
        end
        
        negclust = zeros(size(stat_TFR_micro{imarker}.mask));
        for ineg = 1 : size(stat_TFR_micro{imarker}.negclusters,2)
            if stat_TFR_micro{imarker}.negclusters(ineg).prob < 0.025
                negclust(stat_TFR_micro{imarker}.negclusterslabelmat == ineg) = 1;
            end
        end
        
        d = size(TFR_micro_trials{imarker}.powspctrm,4) - size(stat_TFR_micro{imarker}.mask,3);
        
        TFR_micro_trials{imarker}.mask      = cat(3, zeros(size(stat_TFR_micro{imarker}.mask,1),size(stat_TFR_micro{imarker}.mask,2),d), stat_TFR_micro{imarker}.mask);
        TFR_micro_trials{imarker}.mask      = repmat(TFR_micro_trials{imarker}.mask, size(TFR_micro_trials{imarker}.powspctrm,1),1,1,1);
        TFR_micro_trials{imarker}.mask      = reshape(TFR_micro_trials{imarker}.mask, size(TFR_micro_trials{imarker}.powspctrm));
        
        TFR_micro_trials{imarker}.maskpos   = cat(3, zeros(size(stat_TFR_micro{imarker}.mask,1),size(stat_TFR_micro{imarker}.mask,2),d), posclust);
        TFR_micro_trials{imarker}.maskpos   = repmat(TFR_micro_trials{imarker}.maskpos, size(TFR_micro_trials{imarker}.powspctrm,1),1,1,1);
        TFR_micro_trials{imarker}.maskpos   = reshape(TFR_micro_trials{imarker}.maskpos, size(TFR_micro_trials{imarker}.powspctrm));
        
        TFR_micro_trials{imarker}.maskneg   = cat(3, zeros(size(stat_TFR_micro{imarker}.mask,1),size(stat_TFR_micro{imarker}.mask,2),d), negclust);
        TFR_micro_trials{imarker}.maskneg   = repmat(TFR_micro_trials{imarker}.maskneg, size(TFR_micro_trials{imarker}.powspctrm,1),1,1,1);
        TFR_micro_trials{imarker}.maskneg   = reshape(TFR_micro_trials{imarker}.maskneg, size(TFR_micro_trials{imarker}.powspctrm));
        
        TFR_micro_trials{imarker}.masked    = TFR_micro_trials{imarker}.powspctrm;
        TFR_micro_trials{imarker}.masked(~TFR_micro_trials{imarker}.mask) = nan;
        
        TFR_micro_trials{imarker}.maskedpos = TFR_micro_trials{imarker}.powspctrm;
        TFR_micro_trials{imarker}.maskedpos(~TFR_micro_trials{imarker}.maskpos) = nan;
        
        TFR_micro_trials{imarker}.maskedneg = TFR_micro_trials{imarker}.powspctrm;
        TFR_micro_trials{imarker}.maskedneg(~TFR_micro_trials{imarker}.maskneg) = nan;
        
        powpos = nanmean(TFR_micro_trials{imarker}.maskedpos,1);
        powpos = nanmean(powpos,4);
        powpos = squeeze(powpos);
        
        powneg = nanmean(TFR_micro_trials{imarker}.maskedneg,1);
        powneg = nanmean(powneg,4);
        powneg = squeeze(powneg);
        
        powall = nanmean(TFR_micro_trials{imarker}.powspctrm,1);
        powall = nanmean(powall,4);
        powall = squeeze(powall);
        
        %% Figure
        fig = figure; hold;
        
        plot(TFR_micro_trials{imarker}.freq,powall)
        plot(TFR_micro_trials{imarker}.freq,powpos)
        plot(TFR_micro_trials{imarker}.freq,powneg)
        
        [PKS,LOCS] = findpeaks(powall,'NPeaks',10);
        TFR_micro_trials{imarker}.maxfreq = TFR_micro_trials{imarker}.freq(LOCS);
        TFR_micro_trials{imarker}.maxpow  = PKS;
        
        % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'FFT_',cfg.name{imarker},'.pdf']),'-r600');
        
        cfgtemp                 = [];
        cfgtemp.method          = 'mtmfft';
        cfgtemp.output          = 'pow';
        cfgtemp.taper           = 'hanning';
        cfgtemp.pad             = 'nextpow2';
        cfgtemp.keeptrials      = 'no';
        cfgtemp.foi             = 1:200;
        %     cfgtemp.t_ftimwin       = ones(size(cfgtemp.foi))*0.5;
        cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
        
        % cfgtemp.toi             = cfg.epoch.toi{imarker}(1):cfg.LFP.slidestep:cfg.epoch.toi{imarker}(2);
        FFT_micro_trials{imarker} = ft_freqanalysis(cfgtemp,dat_micro{imarker});
        
        save(fname,'FFT_micro_trials','TFR_micro_trials','TFR_macro_trials','stat_TFR_micro','corrs','-v7.3');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%
        
        % Find representative trials
        
        cfgtemp                 = [];
        cfgtemp.vartrllength    = 2;
        avg                     = ft_timelockanalysis(cfgtemp,dat_micro{imarker});
        clear c
        for itrial = 1 : size(dat_micro{imarker}.trial,2)
            fprintf('correlating trial %d out of %d \n',itrial,size(dat_micro{imarker}.trial,2));
            i1          = find(dat_micro{imarker}.time{itrial} >= cfg.epoch.toi{imarker}(1),1,'first');
            i2          = find(dat_micro{imarker}.time{itrial} >= cfg.epoch.toi{imarker}(2),1,'first');
            c(itrial)   = corr(avg.avg(channelindx,i1:i2)',dat_micro{imarker}.trial{itrial}(channelindx,i1:i2)');
        end
        
        %% FIGURE
        
        fig = figure;
        subplot(5,2,1);
        hold;
        h = 150;
        n = 1;
        [Y, I] = sort(c,'descend');
        for itrial = I(1:10)
            i1 = find(dat_micro{imarker}.time{itrial} >= cfg.epoch.toi{imarker}(1),1,'first');
            i2 = find(dat_micro{imarker}.time{itrial} >= cfg.epoch.toi{imarker}(2),1,'first');
            plot(dat_micro{imarker}.time{itrial}(i1:i2),dat_micro{imarker}.trial{itrial}(channelindx,i1:i2) + n*h,'color','k');
            n = n + 1;
        end
        ylabel('Trials');
        xlabel('Time (s)');
        axis tight
        
        subplot(5,2,2);
        hold;
        h = 150;
        n = 1;
        for itrial = I(1:10)
            i1 = find(dat_macro_chansel.time{itrial} >= cfg.epoch.toi{imarker}(1),1,'first');
            i2 = find(dat_macro_chansel.time{itrial} >= cfg.epoch.toi{imarker}(2),1,'first');
            plot(dat_macro{imarker}.time{itrial}(i1:i2),dat_macro{imarker}.trial{itrial}(1,i1:i2) + n*h,'color','k');
            n = n + 1;
        end
        ylabel('Trials');
        xlabel('Time (s)');
        axis tight
        
        % plot grid of trial voltages (micro)
        subplot(5,2,3);
        image(trialgrid_micro*255);
        colormap hot(255)
        set(gca, 'XTickLabel', '');
        % title('Raw trial amplitudes over time');
        ylabel('Trials');
        axis tight
        ax = axis;
        
        % plot grid of trial voltages (micro)
        subplot(5,2,4);
        image(trialgrid_macro*255);
        %     colormap hot(255)
        colormap parula(255^2)
        
        % title('Raw trial amplitudes over time');
        ylabel('Trials');
        axis tight
        set(gca, 'XTickLabel', '');
        
        % plot TFR micro
        subplot(5,2,5);
        cfgtemp                 = [];
        %     cfgtemp.baseline        = cfg.LFP.baselinewindow{imarker};
        %     cfgtemp.baselinetype    = 'relchange';
        cfgtemp.xlim            = cfg.epoch.toi{imarker};
        cfgtemp.colorbar        = 'no';
        cfgtemp.colorbar        = 'yes';
        cfgtemp.zlim            = 'maxabs';
        cfgtemp.title           = ' ';
        cfgtemp.parameter       = 'powspctrm';
        cfgtemp.colormap        = parula(5000);
        cfgtemp.renderer        = 'painters';
        ft_singleplotTFR(cfgtemp,TFR_micro_trials{imarker});
        
        % title('Time-frequency-representation');
        ylabel('Freq');
        set(gca, 'XTickLabel', '');
        
        % plot TFR macro
        subplot(5,2,6);
        cfgtemp = [];
        %     cfgtemp.baseline = 'yes';
        %     cfgtemp.baselinetype = 'relative';
%         cfgtemp.xlim = cfg.epoch.toi{imarker};
cfgtemp.xlim = [-5 5];
        cfgtemp.colorbar = 'no';
        cfgtemp.colorbar = 'yes';
        cfgtemp.zlim     = 'maxabs';
        cfgtemp.title = ' ';
        ft_singleplotTFR(cfgtemp,TFR_macro_trials{imarker});
        % title('Time-frequency-representation');
        ylabel('Freq');
        set(gca, 'XTickLabel', '');
        
        
        % plot TFR micro MASKED
        subplot(5,2,7);
        cfgtemp                 = [];
        %     cfgtemp.baseline        = cfg.LFP.baselinewindow{imarker};
        %     cfgtemp.baselinetype    = 'relchange';
        cfgtemp.xlim            = cfg.epoch.toi{imarker};
        cfgtemp.colorbar        = 'no';
        cfgtemp.colorbar        = 'yes';
        cfgtemp.zlim            = 'maxabs';
        cfgtemp.title           = ' ';
        cfgtemp.maskstyle       = 'opacity';
        cfgtemp.maskparameter   = 'mask';
        cfgtemp.parameter       = 'powspctrm';
        cfgtemp.colormap        = parula(5000);
        cfgtemp.renderer        = 'painters';
        ft_singleplotTFR(cfgtemp,TFR_micro_trials{imarker});
        % title('Time-frequency-representation');
        ylabel('Freq');
        set(gca, 'XTickLabel', '');
        
        % plot average peak micro
        subplot(5,2,9); hold;
        [~,~,W,P] = findpeaks(dat_micro_rptavg.avg,dat_micro_rptavg.time,'MinPeakProminence',40,'Annotate','extents');
        findpeaks(dat_micro_rptavg.avg,dat_micro_rptavg.time,'MinPeakProminence',40,'Annotate','extents');
        legend('off');
        ylabel('Microvolts');
        [~, ci] = max(P);
        title(sprintf('Peak analysis: W=%.0fms, P=%.0fmV',W(ci)*1000,P(ci)));
        axis tight
        set(gca, 'XTickLabel', '');
        
        % plot correlations with macro electrodes
        subplot(5,2,10); hold;
        boxplot(corrs.c,'orientation','horizontal','outliersize',1);
        
        % % plot FFT micro
        % subplot(6,2,11); hold;
        % plot(FFT_micro_trials.freq,mean(FFT_micro_trials.powspctrm,1),'k');
        % [ymax,imax] = max(mean(FFT_micro_trials.powspctrm,1));
        % line([FFT_micro_trials.freq(imax),FFT_micro_trials.freq(imax)],[0,ymax],'color','r','linewidth',2);
        % txt = ['\leftarrow ', num2str(FFT_micro_trials.freq(imax),3),'Hz'];
        % text(FFT_micro_trials.freq(imax),ymax,txt,'color','r');
        % % title('Spectral analysis');
        % xlabel('Hz');
        % ylabel('Power');
        %
        % % histogram of durations
        % subplot(6,2,12); hold;
        % triallength = zeros(size(dat_micro.trial,2),1);
        % for itrial = 1 : size(dat_micro.trial,2)
        %     triallength(itrial) = dat_micro.time{itrial}(end) - cfg.poststim;
        % end
        % [histcount,edges,bin] = histcounts(triallength,0:cfg.binsize:max(triallength)+1);
        % bar(edges(2:end)-cfg.binsize/2,histcount,1,'facecolor','k')
        % % title('Distribution of durations');
        % ylabel('Nr. of observations');
        % xlabel('Duration (s)');
        
        %
        % [ymax,imax] = max(histcount);
        %
        % [m] = median(triallength);
        %
        % line([m, m],[0,ymax],'color','b','linewidth',2);
        % txt = [num2str(edges(imax)+cfg.binsize/2,3),'s'];
        % text(edges(imax)+cfg.binsize,ymax*0.8,txt,'color','b')
        %
        % line([edges(imax)+cfg.binsize/2,edges(imax)+cfg.binsize/2],[0,ymax],'color','r','linewidth',2);
        % txt = [num2str(edges(imax)+cfg.binsize/2,3),'s'];
        % text(edges(imax)+cfg.binsize,ymax,txt,'color','r')
        %
        % [yleft,ileft] = find(histcount>0,1,'first');
        % txt = ['\downarrow ', num2str(edges(ileft)+cfg.binsize/2) ,'s'];
        % text(edges(ileft)+cfg.binsize,yleft+20,txt,'color','r')
        %
        % [yright,iright] = find(histcount>0,1,'last');
        % txt = ['\downarrow ', num2str(edges(iright)+cfg.binsize/2) ,'s'];
        % text(edges(iright)+cfg.binsize/2,yright+20,txt,'color','r')
        
        % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'overview_analyses_',cfg.name{imarker},'.pdf']),'-r600');
        
    end
    
    
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot single event examples %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % fig = figure;
    % n = 1;
    % for itrial = cfg.representtrials(1:10)
    %
    %     cfgtemp                 = [];
    %     cfgtemp.trials          = itrial;
    %     cfgtemp.method          = 'mtmconvol';
    %     cfgtemp.output          = 'pow';
    %     cfgtemp.taper           = 'hanning';
    %     cfgtemp.pad             = 'nextpow2';
    %     cfgtemp.foi             = 10:1:300;
    %     cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
    %     cfgtemp.toi             = -cfg.prestim:cfg.slidestep:cfg.poststim;
    %     cfgtemp.title           = ['Micro trialnr. ' num2str(itrial)];
    %     cfgtemp.channel         = cfg.channel;
    %     TFR_micro_single        = ft_freqanalysis(cfgtemp,dat_micro);
    %
    %     cfgtemp.channel         = 1;
    %     cfgtemp.title           = ['Micro trialnr. ' num2str(itrial)];
    %     TFR_macro_single        = ft_freqanalysis(cfgtemp,dat_macro);
    %
    %     subplot(size(cfg.representtrials,2)/2,2,n); hold;
    %     cfgtemp = [];
    %     cfgtemp.channel = cfg.channel;
    %     cfgtemp.baseline = 'yes';
    %     cfgtemp.baselinetype = 'relative';
    %     ft_singleplotTFR(cfgtemp,TFR_micro_single);
    %
    %     % add scaled LFP line
    %     hold;
    %     ax = axis;
    %     scaled = (dat_micro.trial{itrial}(chanindx,:) + abs(min(dat_micro.trial{itrial}(chanindx,:))));
    %     scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
    %     i1 = find(dat_micro.time{itrial} > TFR_micro_single.time(1),1,'first');
    %     i2 = find(dat_micro.time{itrial} < TFR_micro_single.time(end),1,'last');
    %     plot(dat_micro.time{itrial}(i1:i2),scaled(i1:i2),'linewidth',1,'color',[1 1 1]);
    %
    %     axis tight
    %
    %     subplot(size(cfg.representtrials,2)/2,2,n+1);
    %     cfgtemp = [];
    %     cfgtemp.channel = 'all';
    %     cfgtemp.baseline = 'yes';
    %     cfgtemp.baselinetype = 'relative';
    %     ft_singleplotTFR(cfgtemp,TFR_macro_single);
    %
    %     % add scaled LFP line
    %     hold;
    %     ax = axis;
    %     scaled = (dat_macro.trial{itrial}(chanindx,:) + abs(min(dat_macro.trial{itrial}(chanindx,:))));
    %     scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
    %     i1 = find(dat_macro.time{itrial} > TFR_macro_single.time(1),1,'first');
    %     i2 = find(dat_macro.time{itrial} < TFR_macro_single.time(end),1,'last');
    %     plot(dat_macro.time{itrial}(i1:i2),scaled(i1:i2),'linewidth',1,'color',[1 1 1]);
    %
    %     n = n + 2;
    %
    %     colormap hot
    % end
    %
    % % print to file
    % set(fig,'PaperOrientation','landscape');
    % set(fig,'PaperUnits','normalized');
    % set(fig,'PaperPosition', [0 0 1 1]);
    % print(fig, '-dpdf', fullfile(cfg.imagesavedir,['single_trials_',cfg.label,'.pdf']),'-r600');
    %
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% calculate correlation LFP x power %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % fname = fullfile(cfg.datasavedir,['freqanalysis_singletrial_micro_',cfg.label,'.mat']);
    % if exist(fname,'file') && cfg.force == false
    %     fprintf('loading %s, which can take a while \n',fname);
    %     load(fname,'TFR_micro_singletrial','lfp_corr');
    % else
    %
    %     % time frequency analysis around peaks
    %     cfgtemp                 = [];
    %     cfgtemp.method          = 'mtmconvol';
    %     cfgtemp.output          = 'pow';
    %     cfgtemp.taper           = 'hanning';
    %     cfgtemp.pad             = 'nextpow2';
    %     cfgtemp.keeptrials      = 'yes';
    %     cfgtemp.foi             = 10:1:300;
    %     cfgtemp.channel         = chanindx;
    %     cfgtemp.t_ftimwin       = ones(size(cfgtemp.foi))*0.5;
    %     %     cfgtemp.t_ftimwin       = 7./cfgtemp.foi;
    %     cfgtemp.toi             = -cfg.prestim:cfg.slidestep:cfg.poststim;
    %     TFR_micro_singletrial   = ft_freqanalysis(cfgtemp,dat_micro);
    %
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
