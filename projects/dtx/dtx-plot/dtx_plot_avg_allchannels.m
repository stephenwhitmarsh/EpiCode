function dat_avg = dtx_plot_avg_allchannels(cfg,data,imarker,saveplot)
%plot avg of all channels indicates in cfg.labels.macro and
%cfg.LFP.emg.


data = data{imarker};
nb_channels = length(cfg.labels.macro);
fig = figure;
hold;

%h automatic setting :
cfgtemp = [];
cfgtemp.channel = cfg.align.channel{imarker};
data_h = ft_selectdata(cfgtemp,data);

for itrial = 1 : size(data_h.trial,2)
%     h_temp_max = max(data_h.trial{itrial}(1,:));
%     h_temp_min = min(data_h.trial{itrial}(1,:));
%     h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
% end
% h = mean(h_temp_amplitude)*2;
t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data_h.fsample; % offset for which t = 0;
h_temp(itrial) = max(data_h.trial{itrial}(1,round(-0.5*data_h.fsample)+t_0: round(0.5*data_h.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
end
h = mean(h_temp);


for ichan = 1:nb_channels
    
    %select channel
    cfgtemp = [];
    cfgtemp.channel = cfg.labels.macro{ichan};
    data_1chan = ft_selectdata(cfgtemp,data);
        
    cfgtemp                  = [];
    data_temp_avg            = ft_timelockanalysis(cfgtemp,data_1chan);
    dat_avg{ichan}           = data_temp_avg;
    
    plot(dat_avg{ichan}.time,dat_avg{ichan}.avg+(nb_channels+1)*h-h*ichan,'k');

    
end

plot([0 0],[0 (nb_channels+1)*h], '--r', 'Linewidth', 1);

axis tight
xlim(cfg.epoch.toi{1});
%ylim([0 (nb_channels+1)*h]);
xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Channel name', 'Fontsize',15);
title(sprintf('%d seizures', length(data.trial)),'Interpreter','none','Fontsize',20);
set(gca, 'FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(h : tick : nb_channels*h);
set(gca,'TickDir','out');
set(gca, 'YTickLabel',[flip(cfg.labels.macro)]);


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
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'avg_allchannels_eeg',cfg.LFP.name{imarker}]),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'avg_allchannels_eeg',cfg.LFP.name{imarker}]),'-r600');
    close all
end

end

