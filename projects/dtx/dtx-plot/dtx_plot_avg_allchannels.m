function [EEG_avg_allchan, EEG_align] = dtx_plot_avg_allchannels(cfg,data,ipart,imarker,saveplot)
%plot avg of all channels indicates in cfg.labels.macro and
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
nb_channels = length(cfg.labels.macro);
fig = figure;
hold;

cfgtemp                  = [];
EEG_avg_allchan          = ft_timelockanalysis(cfgtemp,data);
EEG_avg_allchan.ID  = cfg.prefix(1:end-1);

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

    if strcmp(cfg.labels.macro{ichan},cfg.align.channel{imarker}) %if chan is align channel
        EEG_align = dat_avg{ichan};
        EEG_align.label = 'chan_SlowWave';
        EEG_align.ID = cfg.prefix(1:end-1);
    end
    
end

plot([0 0],[0 (nb_channels+1)*h], '--r', 'Linewidth', 1);

axis tight
xlim(cfg.epoch.toi{1});
%ylim([0 (nb_channels+1)*h]);
xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Channel name', 'Fontsize',15);
title(sprintf('%s : average of %d trials',cfg.LFP.name{imarker}, length(data.trial)),'Interpreter','none','Fontsize',20);
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
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_avg_allchannels_eeg']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_avg_allchannels_eeg']),'-r600');
    close all
end

end

