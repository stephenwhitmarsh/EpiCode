function dtx_plot_comparison_eeg_emg(cfg,data,iEEG,iEMG,imarker,saveplot)
%Plot OL and EMG data from selected channels
%iEEG et iEMG are the indexes of the channel in data.


data = data{imarker};
iEEG = iEEG{imarker};
iEMG = iEMG{imarker};

% à faire avec xlim 5 et 2
abscisse_scale = 2;
fig = figure;
fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap


%% Overdraw OL
subplot(3,1,1);
hold;
fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
for itrial = 1 : size(data.trial,2)
    plot(data.time{itrial},data.trial{itrial}(iEEG,:),'color',[0.6 0.6 0.6]); %first on top
end

title(sprintf('%s : EEG average',cfg.LFP.name{imarker}),'Fontsize',18,'Interpreter','none');
ylabel(sprintf('EEG %s (µV)',data.label{iEEG}),'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
axis tight;
xlim([-abscisse_scale, abscisse_scale]);
%ylim([-150 250]);
set(gca,'ycolor','r');
set(gca,'Fontsize',15);


%% Average OL

cfgtemp                 = [];
cfgtemp.vartrllength    = 2;
data_rptavg        = ft_timelockanalysis(cfgtemp,data);

plot(data_rptavg.time,data_rptavg.avg(iEEG,:),'r','LineWidth', 2);


%% Overdraw and rectified EMG
subplot(3,1,2);
hold;


for itrial = 1 : size(data.trial,2)
    plot(data.time{itrial},data.trial{itrial}(iEMG,:),'color',[0.6 0.6 0.6]); %first on top
end

yyaxis left
ylabel(sprintf('%s (µV)',data.label{iEMG}),'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor',[0.6 0.6 0.6]);

yyaxis right


data_abs                = data;
data_abs.label          = data.label(iEMG);
for itrial = 1:length(data.trial)
    data_abs.trial{itrial}                = abs(data.trial{itrial}(iEMG,:));
end
cfgtemp                 = [];
data_abs_avg            = ft_timelockanalysis(cfgtemp,data_abs);

plot(data_abs_avg.time,data_abs_avg.avg,'b','LineWidth', 2);
%plot(data_abs_avg.time,envelope(data_abs_avg.avg,15,'rms'),'-b','LineWidth', 2);


title('Average of rectified EMG','Fontsize',18);
ylabel(sprintf('%s (rect avg)',data_abs.label{1}),'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);



%% Comparaison EEG-EMG
subplot (3,1,3)
hold;

plot(data_rptavg.time,data_rptavg.avg(iEEG,:),'r','LineWidth', 2);

yyaxis left
set(gca,'ycolor','r');
ylabel(sprintf('EEG %s (µV)',data.label{iEEG}), 'Fontsize',15);

ylim_eeg = get(gca,'ylim');
ylim_rapport = ylim_eeg(2)/ylim_eeg(1); %to set automatically EMG scale with same 0 as eeg

yyaxis right

%plot(data_abs_avg.time,envelope(data_abs_avg.avg,15,'rms'),'b','LineWidth', 2);
plot(data_abs_avg.time,data_abs_avg.avg,'b','LineWidth', 2);

axis tight
xlim([-abscisse_scale, abscisse_scale]);
ylim_emg = get(gca,'ylim');
ylim_emg(1)=ylim_emg(2)/ylim_rapport; %to set automatically EMG scale with same 0 as EEG
ylim(ylim_emg);

title('EEG-EMG comparison','Fontsize',18);
ylabel(sprintf('%s (rect avg)',data.label{iEMG}), 'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');

xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);

xlabel('Time (s)','Fontsize',18);


%% sava data
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'comparisoneegemg_',cfg.LFP.name{imarker},'_',data.label{iEEG},'_',data.label{iEMG},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'comparisoneegemg_',cfg.LFP.name{imarker},'_',data.label{iEEG},'_',data.label{iEMG},'.png']),'-r600');
    close all
end


end

