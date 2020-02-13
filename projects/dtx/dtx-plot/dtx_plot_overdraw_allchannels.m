function dtx_plot_overdraw_allchannels(cfg,data,iEEG,iEMG,separate_R_L,saveplot)
%plot overdraw and avg of all channels indicates in cfg.labels.macro and
%cfg.labels.emg.

%iEEG{imarker} is the index of the channel used for alignment. It is used as
%reference for setting h value.
%if iEMG{imarker} = 'no' : indicate absence of EMG
%if iEMG{imarker} = false : ignoring everything concerning EMG
%if separate_R_L = true, subplot right and left and indicate absence of one
%side of seizure if so

%ATTENTION : if there are more than 2 markers which define trials, this
%script is not yet adapted.


fig = figure;
fig.Renderer='Painters';


for imarker = 1 : length(data)
    if strcmp(iEMG{imarker},'no') || iEMG{imarker} == false
        n = (length(cfg.labels.macro));
    else
        n = (length(cfg.labels.macro)+length(cfg.labels.emg));
    end
    
    %if separate_R_L == false : no subplot and only plot 1 marker
    if separate_R_L == true
        subplot(1,2,imarker);
    end
    hold;
    
    %h automatic setting :
    for itrial = 1 : size(data{imarker}.trial,2)
        h_temp_max = max(data{imarker}.trial{itrial}(iEEG{imarker},:));
        h_temp_min = min(data{imarker}.trial{itrial}(iEEG{imarker},:));
        h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
    end
    h = mean(h_temp_amplitude)*2;

    for ichan = 1:n
        
        %select channel
        cfgtemp = [];
        channelisEMG=false;
        if ichan <= length(cfg.labels.macro)
            cfgtemp.channel = cfg.labels.macro{ichan};
            channelisEMG=false;
        elseif ~(strcmp(iEMG{imarker},'no') || iEMG{imarker} == false) %if there is EMG data
            ichan_emg = ichan-length(cfg.labels.macro);
            cfgtemp.channel = cfg.labels.emg{ichan_emg};
            channelisEMG=true;
        end
        dat_temp = ft_selectdata(cfgtemp,data{imarker});
        
        %plot all trials
        for itrial = 1:length(dat_temp.trial)
            if channelisEMG
                dat_temp.trial{itrial} = dat_temp.trial{itrial}/max(dat_temp.trial{itrial})*h/5; %to have an amplitude comparable to EEG
            end %à régler
            plot(dat_temp.time{itrial},dat_temp.trial{itrial}+(n+1)*h-h*ichan,'color',[0.6 0.6 0.6]); %first on top
        end
        
        %plot avg
        if channelisEMG %avg of abs
            for itrial = 1:length(dat_temp.trial)
                dat_temp.trial{itrial}                = abs(dat_temp.trial{itrial})*3; %a régler
            end
        end
        
        cfgtemp                  = [];
        dat_temp_avg             = ft_timelockanalysis(cfgtemp,dat_temp);
        dat_avg{ichan}           = dat_temp_avg;
        
        
    end
    for ichan = 1: n
        plot(dat_avg{ichan}.time,dat_avg{ichan}.avg+(n+1)*h-h*ichan,'k','LineWidth',2);
    end
    
    plot([0 0],[0 (n+1)*h], '--r', 'Linewidth', 1);
    
    axis tight
    xlim(cfg.epoch.toi{1});
    ylim([0 (n+1)*h]);
    xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
    ylabel('Channel name', 'Fontsize',15);
    title(sprintf('%d seizures', length(data{imarker}.trial)),'Interpreter','none','Fontsize',18);
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    tick = h;
    yticks(h : tick : n*h);
    set(gca,'TickDir','out');
    if strcmp(iEMG{imarker},'no') || iEMG{imarker} == false
        set(gca, 'YTickLabel',[flip(cfg.labels.macro)]);
    else
        set(gca, 'YTickLabel',[flip(cfg.labels.emg), flip(cfg.labels.macro)]); 
    end
    
    
end

% Case separate_R_L == true but there are no controlateral seizures
if  separate_R_L == true && length(data)==1 %script adapted only for 2 markers maximum
    subplot(1,2,2);
    set(gca,'TickLength',[0 0]);
    yticklabels([]);
    xticklabels([]);
    set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
    text(0.2,0.5,sprintf('No controlateral seizure for \n%s', cfg.prefix(1:end-1)),'Interpreter','none','Fontsize',15);
end




%% print to file
if saveplot
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'overdraw_allchannels',[cfg.prefix,'overdraw_allchannels']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,'overdraw_allchannels',[cfg.prefix,'overdraw_allchannels']),'-r600');
    close all
end

end

