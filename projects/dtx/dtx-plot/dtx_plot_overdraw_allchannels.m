function dtx_plot_overdraw_allchannels(cfg,data,ipart,imarker,toi,saveplot)
%plot overdraw and avg of all channels indicates in cfg.labels.macro and
%cfg.LFP.emg.

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

data = data{ipart}{imarker};

nb_channels = length(data.label);
fig = figure;
%subplot(1,2,1)
hold;

%h automatic setting :
cfgtemp = [];
cfgtemp.channel = cfg.LFP.electrodetoplot{imarker};
data_h = ft_selectdata(cfgtemp,data);

%h automatic setting :
for itrial = 1 : length(data.trial)
%     h_temp_max = max(data.trial{itrial}(1,:));
%     h_temp_min = min(data.trial{itrial}(1,:));
%     h_temp(itrial) = h_temp_max - h_temp_min;
t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data.fsample; % offset for which t = 0;
h_temp(itrial) = max(data_h.trial{itrial}(1,round(-0.5*data.fsample)+t_0: round(0.5*data.fsample)+t_0)); %amplitude of peak vs baseline
bl_period_inf = data_h.time{itrial}>cfg.align.toibaseline{imarker}(1);
bl_period_sup = data_h.time{itrial}<cfg.align.toibaseline{imarker}(2);
bl_period = logical(bl_period_inf .* bl_period_sup);
h_temp(itrial) = max(data_h.trial{itrial}(1, data_h.time{itrial}>-1 & data.time{itrial}<1)) - nanmean(data_h.trial{itrial}(1, bl_period));
end

h = mean(h_temp);


for ichan = 1:nb_channels
    
    %select channel
    cfgtemp = [];
    cfgtemp.channel = data.label{ichan};
    data_1chan = ft_selectdata(cfgtemp,data);
    
    %plot all trials
    for itrial = 1:length(data_1chan.trial)
        plot(data_1chan.time{itrial},data_1chan.trial{itrial}+(nb_channels+1)*h-h*ichan,'color',[0.6 0.6 0.6]); %first on top
    end

    
    cfgtemp                  = [];
    data_temp_avg            = ft_timelockanalysis(cfgtemp,data_1chan);
    dat_avg{ichan}           = data_temp_avg;
    
end

for ichan = 1: nb_channels
    plot(dat_avg{ichan}.time,dat_avg{ichan}.avg+(nb_channels+1)*h-h*ichan,'k','LineWidth',2);
end

plot([0 0],[0 (nb_channels+1)*h], '--r', 'Linewidth', 1);

axis tight
xlim(toi);
ylim([0 (nb_channels+1)*h]);
xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Channel name', 'Fontsize',15);
title(sprintf('%s : overdraw of %d trials',cfg.LFP.name{imarker}, length(data.trial)),'Interpreter','none','Fontsize',18);
set(gca, 'FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(h : tick : nb_channels*h);
set(gca,'TickDir','out');
set(gca, 'YTickLabel',[flip(data.label')]);


%% print to file
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end

    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_allchannels_eeg','_scale[',num2str(toi),']']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_overdraw_allchannels_eeg','_scale[',num2str(toi),']']),'-r600');
    close all
end

end

