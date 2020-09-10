function dtx_plot_overdraw_allchannels(cfg,data,iEEG,iEMG,separateMarkers)
%iEEG et iEMG are the indexes of the channels in data.
%if iEMG = 'no' : indicate absence of EMG
%if iEMG = false : ignoring everything concerning EMG
%if separateMarkers = true, 

%abscisse_scale = 10;%s

fig = figure;
fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap

%% EEG
if ~(iEMG == false)
    subplot(1,2,1);
end
hold;

%h automatic setting :
for itrial = 1 : size(data{imarker}.trial,2)
    h_temp_max = max(data{imarker}.trial{itrial}(iEEG,:));
    h_temp_min = min(data{imarker}.trial{itrial}(iEEG,:));
    h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
end
h = mean(h_temp_amplitude);

plot([0 0],[0 size(data{imarker}.trial,2)*h+h], 'r', 'Linewidth', 2);
%plot([0 0],[0 size(data{imarker}.trial,2)*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
for itrial = 1 : size(data{imarker}.trial,2)
    plot(data{imarker}.time{itrial},data{imarker}.trial{itrial}(iEEG,:)+ itrial*h,'k'); %first on top
end


xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('Aligned data from %s', data{imarker}.label{iEEG}),'Interpreter','none','Fontsize',18);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(tick : tick*10 : size(data{imarker}.trial,2)*h);
yticklabels(size(data{imarker}.trial,2) : -10 : 0);
set(gca,'TickDir','out');
axis tight
xlim(cfg.epoch.toi{1});

%% EMG
if ~(iEMG == false)
    
    subplot(1,2,2);
    hold;
    
    if ~(iEMG=='no')
        
        
        %h automatic setting :
        for itrial = 1 : size(data{imarker}.trial,2)
            h_temp_max = max(data{imarker}.trial{itrial}(iEMG,:));
            h_temp_min = min(data{imarker}.trial{itrial}(iEMG,:));
            h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
        end
        h = mean(h_temp_amplitude);
        
        
        
        %plot([0 0],[0 size(data{imarker}.trial,2)*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
        plot([0 0],[0 size(data{imarker}.trial,2)*h+h], 'r', 'Linewidth', 2);
        
        for itrial = 1 : size(data{imarker}.trial,2)
            plot(data{imarker}.time{itrial},data{imarker}.trial{itrial}(iEMG,:)+ itrial*h,'k'); %first on top
        end
        
        
        xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('Number of seizures', 'Fontsize',15);
        title(sprintf('Aligned data from %s', data{imarker}.label{iEMG}),'Interpreter','none','Fontsize',18);
        set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
        tick = h;
        yticks(tick : tick*10 : size(data{imarker}.trial,2)*h);
        yticklabels(size(data{imarker}.trial,2) : -10 : 0);
        set(gca,'TickDir','out');
        axis tight
        xlim(cfg.epoch.toi{1});
        
        
    else
        set(gca,'TickLength',[0 0]);
        yticklabels([]);
        xticklabels([]);
        set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
        text(0.2,0.5,sprintf('No EMG data \nassociated with %s', cfg.LFP.name{imarker}),'Interpreter','none','Fontsize',15);
    end
end

%% print to file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
set(fig,'Renderer','Painters');
print(fig, '-dpdf', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG},'_',data{imarker}.label{iEMG}]),'-r600');
print(fig, '-dpng', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG},'_',data{imarker}.label{iEMG}]),'-r600');
close all

end

