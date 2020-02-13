function dtx_plot_timecourse_eeg_emg(cfg,data,iEEG,iEMG,imarker,saveplot)
%iEEG{imarker} et iEMG{imarker} are the indexes of the channels in data.
% if iEMG{imarker} = 'no' : subplot indicating absence of EMG
%if iEMG{imarker} = false : ignoring everything concerning EMG

%abscisse_scale = 10;%s

n = size(data{imarker}.trial,2);
fig = figure;
fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap


%% EEG
subplot(1,2,1);

if ~strcmp(iEMG{imarker},'no') %if iEMG='no' : do subplot(1,2,1)
    if iEMG{imarker} == false %if iEMG=false : do not subplot
        subplot(1,1,1);
    end
end
hold;

%h automatic setting :
for itrial = 1 : n
    h_temp_max = max(data{imarker}.trial{itrial}(iEEG{imarker},:));
    h_temp_min = min(data{imarker}.trial{itrial}(iEEG{imarker},:));
    h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
end
h = mean(h_temp_amplitude);

plot([0 0],[0 n*h+h], 'r', 'Linewidth', 2);
%plot([0 0],[0 n*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
for itrial = 1 : n
    plot(data{imarker}.time{itrial},data{imarker}.trial{itrial}(iEEG{imarker},:)+ (n+1)*h - itrial*h,'k'); %first on top
end


xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('Aligned data from %s', data{imarker}.label{iEEG{imarker}}),'Interpreter','none','Fontsize',18);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(tick : tick*10 : n*h);
yticklabels(n : -10 : 0);
set(gca,'TickDir','out');
axis tight
xlim(cfg.epoch.toi{1});

%% EMG

if ~strcmp(iEMG{imarker},'no')
    if ~(iEMG{imarker} == false)
        
        subplot(1,2,2);
        hold;
        
        %h automatic setting :
        for itrial = 1 : n
            h_temp_max = max(data{imarker}.trial{itrial}(iEMG{imarker},:));
            h_temp_min = min(data{imarker}.trial{itrial}(iEMG{imarker},:));
            h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
        end
        h = mean(h_temp_amplitude);
        
        
        
        %plot([0 0],[0 n*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
        plot([0 0],[0 n*h+h], 'r', 'Linewidth', 2);
        
        for itrial = 1 : n
            t = data{imarker}.time{itrial};
            plot(t,data{imarker}.trial{itrial}(iEMG{imarker},:)+ (n+1)*h - itrial*h,'k'); %first on top
        end
        
        
        
        xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('Number of seizures', 'Fontsize',15);
        title(sprintf('Aligned data from %s', data{imarker}.label{iEMG{imarker}}),'Interpreter','none','Fontsize',18);
        set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
        tick = h;
        yticks(tick : tick*10 : n*h);
        yticklabels(n : -10 : 0);
        set(gca,'TickDir','out');
        axis tight
        xlim(cfg.epoch.toi{1});
    end
end

if strcmp(iEMG{imarker},'no')
    subplot(1,2,2);
    set(gca,'TickLength',[0 0]);
    yticklabels([]);
    xticklabels([]);
    set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
    text(0.2,0.5,sprintf('No EMG data \nassociated with %s', cfg.LFP.name{imarker}),'Interpreter','none','Fontsize',15);
end


%% print to file
if saveplot
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    
    if strcmp(iEMG{imarker},'no')|| iEMG{imarker} == false
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG{imarker}}]),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG{imarker}}]),'-r600');
    else
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG{imarker}},'_',data{imarker}.label{iEMG{imarker}}]),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,'timecourse',[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG{imarker}},'_',data{imarker}.label{iEMG{imarker}}]),'-r600');
    end
    
    close all
end

end

