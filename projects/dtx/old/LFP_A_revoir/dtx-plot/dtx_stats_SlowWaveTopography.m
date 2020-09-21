function begin = dtx_stats_SlowWaveTopography(cfg,isPatient,data,ipart,imarker,saveplot)
%plot topography of event related to one or two markers, according to 1020 human EEG layout.
%One figure per marker

%plots eeg
%plot topoplot

% no baseline correction : already done for each trial before calling this function

%rename prefix in case of "merge" data
if isfield(cfg, 'merge')
    if cfg.merge == true
        if ipart > 1 && ipart == length(cfg.directorylist) %last part = merge (except if only one part, nothing to merge)
            cfg.prefix = [cfg.prefix, 'MERGED-'];
        else
            cfg.prefix = [cfg.prefix, cfg.directorylist{ipart}{:}, '-'];
        end
    end
end

%% get data

data = data{ipart}{imarker};

%remove EMG channels
cfgtemp                     = [];
cfgtemp.channel             = {'all','-*EMG*'};% {'eeg1020'};
data                        = ft_selectdata(cfgtemp,data);

%bl correctop,
cfgtemp = [];
cfg.demean = 'yes';
cfg.baselinewindow = cfg.stats.toibaseline{imarker};
data = ft_preprocessing(cfgtemp,data);
% avg
cfgtemp = [];
if isPatient; cfgtemp.channel = {'EEG'}; end
cfgtemp.keeptrials = 'no';
dat_EEG_avg = ft_timelockanalysis(cfgtemp,data);

% structure for timelock analysis
cfgtemp = [];
if isPatient; cfgtemp.channel = {'EEG'}; end
cfgtemp.keeptrials = 'yes';
dat_EEG_pour_stats = ft_timelockanalysis(cfgtemp,data);

% create baseline trials
t1_bl_indx                  = find(dat_EEG_avg.time > cfg.stats.toibaseline{imarker}(1),1,'first');
t2_bl_indx                  = find(dat_EEG_avg.time < cfg.stats.toibaseline{imarker}(2),1,'last');
dat_EEG_bl                  = dat_EEG_pour_stats;
dat_EEG_bl.trial            = [];
dat_EEG_bl.trial = ones(size(dat_EEG_pour_stats.trial)) .* nanmean(dat_EEG_pour_stats.trial(:,:,t1_bl_indx:t2_bl_indx),3);%replace baseline by mean, all trials, all channels, bl period

%create trials with real bl (no mean). Dupplicate baseline until it as the
%same size as the whole data.
% size_bl                                         = size(t1_bl_indx:t2_bl_indx,2);
% size_data                                       = size(dat_EEG_pour_stats.trial,3);
% dat_EEG_bl.trial                                = dat_EEG_pour_stats.trial(:,:,t1_bl_indx:t2_bl_indx);
% while size(dat_EEG_bl.trial,3) + size_bl <= size_data
%     dat_EEG_bl.trial(:,:,end+1:end+size_bl)     = dat_EEG_pour_stats.trial(:,:,t1_bl_indx:t2_bl_indx);
% end
% missing_points                                  = size_data - size(dat_EEG_bl.trial,3);
% dat_EEG_bl.trial(:,:,end+1:end+missing_points)  = dat_EEG_pour_stats.trial(:,:,t1_bl_indx:t1_bl_indx+missing_points-1);

imagesavedir = cfg.imagesavedir;
for itest_stat = 1:2
    
    if itest_stat == 1
        latency = 'all';
        idx_begin_stat = 1;
        cfg.imagesavedir = fullfile(imagesavedir, 'stats_all_period');
    elseif itest_stat == 2
        latency = cfg.align.toiactive{1};
        idx_begin_stat = find(dat_EEG_avg.time > cfg.align.toiactive{1}(1),1, 'first');
        cfg.imagesavedir = fullfile(imagesavedir, 'stats_active_period');
    end
    %% compute stats
    %     load('elec1020_neighb.mat','neighbours'); %neighbours file from fieldtrip
    % Found here : https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
    
    %creating neighbours just for fill a mandatory field of timelockstats,
    %set each electrode as isolated (neigbour : only the electrode itself)
    for ichan = 1:length(data.label)
        neighbours(ichan).label         = data.label{ichan};
        neighbours(ichan).neighblabel   = data.label(ichan);
    end
    
    % statistics
    cfgtemp = [];
    cfgtemp.statistic                       = 'ft_statfun_depsamplesT';
    cfgtemp.clusterstatistic                = 'maxsum';
    cfgtemp.alpha                           = cfg.stats.alpha;
    cfgtemp.clusteralpha                    = 0.01;
    cfgtemp.minnbchan                       = 0; %no analysis on neighbours, but the cfg.neighbours field is still mandatory
    cfgtemp.neighbours                      = neighbours;
    cfgtemp.tail                            = 0; %two-sided test
    cfgtemp.method                          = 'montecarlo';
    cfgtemp.computestat                     = 'yes';
    cfgtemp.correctm                        = 'cluster';
    cfgtemp.latency                         = latency;%'all';%cfg.align.toiactive{1}; % active period
    cfgtemp.design(1,:)                     = [ones(1,size(dat_EEG_pour_stats.trial,1)) ones(1,size(dat_EEG_bl.trial,1)) *2];
    cfgtemp.design(2,:)                     = [1 : size(dat_EEG_pour_stats.trial,1) 1 : size(dat_EEG_bl.trial,1)];
    cfgtemp.ivar                            = 1; %1 = active, 2 = bl
    cfgtemp.uvar                            = 2; %trials nr, for each condition (active, bl)
    cfgtemp.numrandomization                = 50;
    cfgtemp.layout                          = 'EEG1020';
    cfgtemp.channel                         = 'all';
    
    stat                             = ft_timelockstatistics(cfgtemp,dat_EEG_pour_stats,dat_EEG_bl);
    
    
    %% plot EEG chan by chan with positive and negative clusters
    %h automatic setting :
    cfgtemp = [];
    cfgtemp.channel = cfg.LFP.electrodetoplot{imarker};
    data_h = ft_selectdata(cfgtemp,data);
    for itrial = 1 : size(data_h.trial,2)
        t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data.fsample; % offset for which t = 0;
        h_temp(itrial) = max(data_h.trial{itrial}(1,round(-0.5*data.fsample)+t_0: round(0.5*data.fsample)+t_0)); %amplitude of peak vs baseline
        bl_period_inf = data_h.time{itrial}>cfg.stats.toibaseline{imarker}(1);
        bl_period_sup = data_h.time{itrial}<cfg.stats.toibaseline{imarker}(2);
        bl_period = logical(bl_period_inf .* bl_period_sup);
        h_temp(itrial) = max(data_h.trial{itrial}(1, data_h.time{itrial}>-1 & data.time{itrial}<1)) - nanmean(data_h.trial{itrial}(1, bl_period));
    end
    h = mean(h_temp);
    
    fig1 = figure;
    hold;
    
    nb_channels = length(data.label);
    % re-order channels according to config, for plot
    [ichan_config,ichan_data] = match_str(cfg.labels.macro,dat_EEG_avg.label);
    
    h_amplitude = 0;
    
    for ichan = ichan_data'
        h_amplitude = h_amplitude +1;
        
        %plot avg
        plot(dat_EEG_avg.time,dat_EEG_avg.avg(ichan,:)+(nb_channels+1)*h-h*h_amplitude,'k');
        
        %select positive clusters and plot
        if isfield(stat,'posclusters')
            for ipos = 1 : size(stat.posclusters,2)
                if stat.posclusters(ipos).prob < cfg.stats.alpha
                    sel = [];
                    sel = find(stat.posclusterslabelmat(ichan,:) == ipos);
                    plot(stat.time(sel),dat_EEG_avg.avg(ichan,round(sel + idx_begin_stat - 2))+(nb_channels+1)*h-h*h_amplitude,'.g','linewidth',2); % +bl index to begin in the active period if any (here is 'all' so no need)
                end
            end
        end
        
        %select negative clusters and plot
        if isfield(stat,'negclusters')
            for ineg = 1 : size(stat.negclusters,2)
                if stat.negclusters(ineg).prob < cfg.stats.alpha
                    sel = [];
                    sel = find(stat.negclusterslabelmat(ichan,:) == ineg);
                    plot(stat.time(sel),dat_EEG_avg.avg(ichan,round(sel + idx_begin_stat - 2))+(nb_channels+1)*h-h*h_amplitude,'.r','linewidth',2);
                end
            end
        end
    end %ichan
    
    
    plot([0 0],[0 (nb_channels+1)*h], '--', 'Linewidth', 0.5, 'Color', [0.6 0.6 0.6]);
    
    axis tight
    xlim(cfg.stats.toiploteeg);
    %ylim([0 (nb_channels+1)*h]);
    xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
    ylabel('Channel name', 'Fontsize',15);
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    title(sprintf('%s : average of %d trials',cfg.LFP.name{imarker}, length(data.trial)),'Interpreter','none','Fontsize',20);
    set(gca,'TickDir','out');
    tick = h;
    yticks(h : tick : nb_channels*h);
    set(gca, 'YTickLabel',[flip(cfg.labels.macro')]);
    
    
    
    %% Topoplot
    %for now, only for positive clusters
    %a star on the topoplot : all the samples af the time interval are significant
    
    
    
    %find significant clusters
    pos_cluster_pvals   = [stat.posclusters(:).prob];
    pos_signif_clust    = find(pos_cluster_pvals < cfg.stats.alpha);
    %Boolean matrix with, for each sample, true/false if significant at this sample
    pos_idx             = ismember(stat.posclusterslabelmat, pos_signif_clust);
    
    %find timings and indexes for the different subplots
    timestep                    = cfg.topoplot.timestep;  % timestep between time windows for each subplot (in seconds)
    toi                         = cfg.topoplot.toi;       % begin and end of topoplots
    subplot_lim                 = [toi(1):timestep:toi(2)];
    t1_toi_indx                 = find(dat_EEG_avg.time > toi(1),1,'first');
    t2_toi_indx                 = find(dat_EEG_avg.time < toi(2),1,'last');
    subplot_lim_idx             = [t1_toi_indx:timestep*data.fsample:t2_toi_indx];
    
    [ichan_dataavg,ichan_stat] = match_str(dat_EEG_avg.label, stat.label); %be sure labels are in the same order in avg and in stat
    
%     if isPatient
%         %Find zlim according to values at t=0
%         figtemp                 = figure;
%         cfgtemp                 = [];
%         cfgtemp.layout          = 'EEG1020';
%         cfgtemp.zlim            = 'maxabs';
%         cfgtemp.xlim            = [-0.125 0.125];
%         ft_topoplotER(cfgtemp,dat_EEG_avg);
%         zaxis_t0                = get(gca,'CLim');
%         close(figtemp);
%         
%         fig2 = figure;
%         
%         for isubplot = 1:length(subplot_lim_idx)-1
%             subplot(8,ceil(length(subplot_lim_idx)/8),isubplot)
%             
%             %re-initialize for each loop
%             positive_chans = zeros(numel(dat_EEG_avg.label),1);
%             positive_chans(ichan_dataavg) = all(pos_idx(ichan_stat, subplot_lim_idx(isubplot):subplot_lim_idx(isubplot+1)), 2);
%             %     positive_chans(ichan1) = any(pos_idx(ichan2, subplot_lim_idx(isubplot):subplot_lim_idx(isubplot+1)), 2);
%             
%             cfgtemp                     = [];
%             cfgtemp.layout              = 'EEG1020';
%             cfgtemp.colorbar            = 'no';
%             cfgtemp.zlim                = zaxis_t0;
%             cfgtemp.xlim                = [subplot_lim(isubplot) subplot_lim(isubplot+1)];
%             cfgtemp.fontsize            = 8;
%             cfgtemp.interactive         = 'no';
%             cfgtemp.renderer            = 'painters';
%             cfgtemp.comment             = 'xlim';
%             cfgtemp.commentpos          = 'title';
%             cfgtemp.highlight           = 'on';
%             cfgtemp.highlightchannel    = find(positive_chans);
%             cfgtemp.highlightcolor      = 'k';
%             ft_topoplotER(cfgtemp,dat_EEG_avg);
%         end
%         
%     end
    %% plot overdraw of each positive electrode
    
    %be sure than channels have the same index
    [ichan_dataavg,ichan_stat] = match_str(dat_EEG_avg.label, stat.label); %be sure labels are in the same order in avg and in stat
    %find active period indexes
    %remove bl indexes because stats indexes are only the active indexes
    t1_activeperiod_idx             = find(dat_EEG_avg.time > cfg.stats.toilatency(1),1,'first') - idx_begin_stat + 1;
    t2_activeperiod_idx             = find(dat_EEG_avg.time < cfg.stats.toilatency(2),1,'last') - idx_begin_stat + 1;
    %search chans with at least one positive sample in this period
    positive_chans = zeros(numel(dat_EEG_avg.label),1);
    positive_chans(ichan_dataavg) = any(pos_idx(ichan_stat, t1_activeperiod_idx:t2_activeperiod_idx), 2);
    
    fig3 = figure;
    hold;
    
    if sum(positive_chans) >= 1
        
        plot(dat_EEG_avg.time,dat_EEG_avg.avg(logical(positive_chans),:), 'LineWidth', 2);
        
        axis tight
        xlim([-1 1]);
        ax = axis;
        plot([0 0], [ax(3) ax(4)], '--', 'Linewidth', 1, 'Color', [0.6 0.6 0.6]);
        xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('µV', 'Fontsize',15);
        set(gca, 'FontWeight','bold', 'Fontsize',15);
        title(sprintf('%s : average of %d trials',cfg.LFP.name{imarker}, length(data.trial)),'Interpreter','none','Fontsize',20);
        set(gca,'TickDir','out');
        legend(dat_EEG_avg.label{logical(positive_chans)},'location','northeast','fontsize',10);
        
    end
    
    fig4 = figure;
    hold;
    
    if sum(positive_chans) >= 1
%         norm_values = max(dat_EEG_avg.avg(logical(positive_chans),:),[],2);
        toi = [-1 1];
        t = dat_EEG_avg.time;
        plot(t(t>toi(1)&t<toi(2)),normalize(dat_EEG_avg.avg(logical(positive_chans),t>toi(1)&t<toi(2)),2,'range'), 'LineWidth', 2);
        
        axis tight
        xlim(toi);
        ylim([0 1]);
        ax = axis;
        plot([0 0], [ax(3) ax(4)], '--', 'Linewidth', 1, 'Color', [0.6 0.6 0.6]);
        xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('normalized', 'Fontsize',15);
        set(gca, 'FontWeight','bold', 'Fontsize',15);
        title(sprintf('%s : average of %d trials - NORMALIZED',cfg.LFP.name{imarker}, length(data.trial)),'Interpreter','none','Fontsize',20);
        set(gca,'TickDir','out');
        legend(dat_EEG_avg.label{logical(positive_chans)},'location','northeast','fontsize',10);
%         normalize
    end
    
    %     %% movie
    %     fig2 = figure;
    %     cfgtemp = [];
    %     cfgtemp.layout          = 'EEG1020';
    %     cfgtemp.framespersec    = 1;
    %     cfgtemp.samperframe     = 256;
    %     cfgtemp.colorbar        = 'yes';
    %     cfgtemp.zlim            = zaxis_t0;
    %     cfgtemp.xlim            = [-3 10];
    %     %cfgtemp.samperframe     = 100;
    %     ft_movieplotER(cfgtemp,dat_EEG_avg);
    
    
    
    
    
    %% print to file
    if saveplot
        
        if ~(exist(cfg.imagesavedir)==7)
            mkdir(cfg.imagesavedir);
            fprintf('Create folder %s\n',cfg.imagesavedir);
        end
        
        set(fig1,'PaperOrientation','landscape');
        set(fig1,'PaperUnits','normalized');
        set(fig1,'PaperPosition', [0 0 1 1]);
        set(fig1,'Renderer','Painters');
        print(fig1, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_avg_allchans_stats']),'-r600');
        print(fig1, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_avg_allchans_stats']),'-r600');
        
%         if isPatient
%             set(fig2,'PaperOrientation','landscape');
%             set(fig2,'PaperUnits','normalized');
%             set(fig2,'PaperPosition', [0 0 1 1]);
%             set(fig2,'Renderer','Painters');
%             print(fig2, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_topography_timecourse_stats']),'-r600');
%             print(fig2, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_topography_timecourse_stats']),'-r600');
%         end
        set(fig3,'PaperOrientation','landscape');
        set(fig3,'PaperUnits','normalized');
        set(fig3,'PaperPosition', [0 0 1 1]);
        set(fig3,'Renderer','Painters');
        print(fig3, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_positive_chans']),'-r600');
        print(fig3, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_positive_chans']),'-r600');
        
        set(fig4,'PaperOrientation','landscape');
        set(fig4,'PaperUnits','normalized');
        set(fig4,'PaperPosition', [0 0 1 1]);
        set(fig4,'Renderer','Painters');
        print(fig4, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_positive_chans_NORMALIZED']),'-r600');
        print(fig4, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_positive_chans_NORMALIZED']),'-r600');
        
      
        close all
    end
end %itest_stats
end


