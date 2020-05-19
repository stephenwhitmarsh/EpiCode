function [sdf, clusterstats] = plot_spikedensity_stats(cfg,SpikeTrials, ipart, ilabel, saveplot)
% En tête a faire

%% Compute spike density 
%trialinfo: column 3 StartSample 4 EndSample

cfgtemp                         = [];
cfgtemp.fsample                 = cfg.spike.resamplefs{ilabel};   % sample at 1000 hz
cfgtemp.keeptrials              = 'yes';
cfgtemp.latency                 = 'maxperiod';%[cfg.epoch.toi{ilabel}(1), cfg.epoch.toi{ilabel}(2)];
cfgtemp.timwin                  = [-1/cfg.spike.resamplefs{ilabel}/2 1/cfg.spike.resamplefs{ilabel}/2];%timewin same length as period
[sdf_event]                     = ft_spikedensity(cfgtemp,SpikeTrials{ipart}{ilabel});

%% Compute cluster stats per unit for events

% prepare dummy data with baseline value per trial for stats
slim(1)                         = find(sdf_event.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
slim(2)                         = find(sdf_event.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
sdf_bl                          = sdf_event;
sdf_bl.trial                    = ones(size(sdf_event.trial)) .* nanmean(sdf_event.trial(:,:,slim(1):slim(2)),3); % replace with mean


for i_unit = 1 : size(SpikeTrials{ipart}{ilabel}.label, 2)
    
    % statistics
    cfgtemp = [];
    cfgtemp.statistic                       = 'ft_statfun_depsamplesT';
    cfgtemp.alpha                           = cfg.stats.alpha;
    cfgtemp.clusteralpha                    = 0.01;
    cfgtemp.method                          = 'montecarlo';
    cfgtemp.computestat                     = 'yes';
    cfgtemp.correctm                        = 'cluster';
    cfgtemp.latency                         = [cfg.stats.bltoi{ilabel}(2) sdf_event.time(end)]; % active period starts after baseline
    cfgtemp.ivar                            = 1;
    cfgtemp.uvar                            = 2;
    cfgtemp.design(1,:)                     = [ones(1,size(sdf_event.trial,1)) ones(1,size(sdf_bl.trial,1)) *2];
    cfgtemp.design(2,:)                     = [1 : size(sdf_event.trial,1) 1 : size(sdf_bl.trial,1)];
    cfgtemp.numrandomization                = cfg.stats.numrandomization;
    cfgtemp.channel                         = i_unit;
    stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit} = ft_timelockstatistics(cfgtemp,sdf_smooth_event,sdf_bl);
    
    % calculate baseline
    slim(1) = find(sdf_event.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
    slim(2) = find(sdf_event.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
    stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.avg        = nanmean(sdf_event.avg(i_unit,slim(1):slim(2)),2);
    stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.var        = nanmean(sdf_event.var(i_unit,slim(1):slim(2)),2);
    stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.dof        = nanmean(sdf_event.dof(i_unit,slim(1):slim(2)),2);
    stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.trialavg   = nanmean(sdf_event.trial(:,i_unit,slim(1):slim(2)),3);
end

bar(sdf_event.time,sdf_event.avg(i_unit,:),1,'edgecolor','none','facecolor',[0 0 0]);
plot(sdf_event.time,sqrt(sdf_event.var(i_unit,:))+sdf_event.avg(i_unit,:),'-','LineWidth',2,'Color',[0.6 0.6 0.6]);

ylabel('Spikerate (Hz)');
%         xlabel(sprintf('Time from %s (s)', cfg.name{ilabel}));
axis tight
ax = axis;

baseline = stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.avg;
idx_begin_stat = find(sdf_event.time > cfg.stats.bltoi{ilabel}(2), 1, 'first');

% plot positive clusters
if isfield(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit},'posclusters')
    for ipos = 1 : size(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.posclusters,2)
        if stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.posclusters(ipos).prob < cfg.stats.alpha
            sel = [];
            sel = find(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.posclusterslabelmat == ipos);
            bar(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),1,'facecolor','g','edgecolor','g');
            
            % plot percentage
            [~,max_idx_temp] = max(sdf_event.avg(i_unit,sel+idx_begin_stat-2));
            max_idx = sel(max_idx_temp) + idx_begin_stat - 2;
            x_max = sdf_event.time(max_idx);
            y_max = sdf_event.avg(i_unit,max_idx);
            percent_increase = (sdf_event.avg(i_unit,max_idx)-baseline) / baseline * 100;
            y = ax(4)*0.9;
            text(x_max,y,sprintf('max = %.1fHz (+%.1f%%)',y_max,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',7);
            plot(x_max,y_max,'*b');
            
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.perc{ipos} = percent_increase;
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.x{ipos} = x_max;
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.y{ipos} = y_max;
        end
    end
end

% plot negative clusters
if isfield(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit},'negclusters')
    for ineg = 1 : size(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.negclusters,2)
        if stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.negclusters(ineg).prob < cfg.stats.alpha
            sel = [];
            sel = find(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.negclusterslabelmat == ineg);
            bar(stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),1,'facecolor','r','edgecolor','r');
            
            % plot percentage
            [~,min_idx_temp] = min(sdf_event.avg(i_unit,sel+idx_begin_stat-2));
            min_idx = sel(min_idx_temp) + idx_begin_stat - 2;
            x_min = sdf_event.time(min_idx);
            y_min = sdf_event.avg(i_unit,min_idx);
            percent_decrease = (sdf_event.avg(i_unit,min_idx)-baseline) / baseline * 100;
            %                     text(x_min,y_min+y_min/100,sprintf('%.1fHz (%.1f%%)\n',y_min,percent_decrease),'HorizontalAlignment','center','VerticalAlignment','middle');
            %                     plot(x_min,y_min,'*b');
            
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.perc{ipos} = percent_decrease;
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.x{ipos} = x_min;
            stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.maxcluster.y{ipos} = y_min;
        end
    end
end

% plot baseline patch
x = [cfg.stats.bltoi{ilabel}(1) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(1)];
y = [ax(3) ax(3) ax(4) ax(4)];
patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);

% plot baseline
x = [];
x = [ax(1) ax(2)];
plot(x,[baseline, baseline],':k');

% plot baseline text
x = (cfg.stats.bltoi{ilabel}(1) + cfg.stats.bltoi{ilabel}(2))/2;
y = ax(4)*0.9;
d = stats{ipart}.(cfg.name{ilabel}).clusterstat{i_unit}.bl.avg;
text(x,y,sprintf('baseline = %.1f Hz',d),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',7);

legend('Average','SD');
set(gca, 'FontWeight', 'bold');

%% sava data
if saveplot
    
    set(gca,'Fontsize',15);
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprintf('Create forlder %s',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'part',num2str(ipart),'_',cfg.name{ilabel},'_',channame,'_Spikedensity_stats.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'part',num2str(ipart),'_',cfg.name{ilabel},'_',channame,'_Spikedensity_stats.png']),'-r600');
    close all
end


end

