function dtx_plot_comparison_eeg_emg(cfg,data,iEEG,iEMG,imarker)
%Plot OL and EMG data from selected channels
%iEEG et iEMG are the indexes of the channel in data.
% if iEMG = false : ignore the steps of the function which need EMG data and do only EEG/LFP analysis


data = data{imarker};
iEEG = iEEG(imarker);
iEMG = iEMG(imarker);

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


%% Overdraw and TFR EMG
subplot(3,1,2);
hold;

%yyaxis left

fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
for itrial = 1 : size(data.trial,2)
    plot(data.time{itrial},data.trial{itrial}(iEMG,:),'color',[0.6 0.6 0.6]); %first on top
end

yyaxis left
ylabel('EMG (µV)','Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
%axis tight
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor',[0.6 0.6 0.6]);





% TFR EMG
%yyaxis left
cfgtemp                         = [];
cfgtemp.channel                 = iEMG;%'all'; %ichannel;
cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
cfgtemp.output                  = 'pow';
cfgtemp.taper                   = 'hanning';
cfgtemp.pad                     = 'nextpow2';
cfgtemp.keeptrials              = 'yes'; %return average
cfgtemp.foi                     = 80:0.1:125;
cfgtemp.t_ftimwin               = 0.4.*true(1,length(cfgtemp.foi));%40./cfgtemp.foi;
%cfgtemp.t_ftimwin               = 20./cfgtemp.foi;
%cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
%cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;

cfgtemp.toi                     = [-5:0.01:25];
TFR_macro                       = ft_freqanalysis(cfgtemp,data);


% TFR_macro_log = TFR_macro;
% TFR_macro_log.powspctrm = log(TFR_macro.powspctrm);
%
% cfgtemp = [];
% cfgtemp.colorbar        = 'no';
% cfgtemp.baseline        = [-4, -2];
% % cfgtemp.baselinetype    = 'relchange';
% % cfgtemp.zlim            = 'maxabs';
% cfgtemp.parameter       = 'powspctrm';
% cfgtemp.colormap        = parula(5000);
% cfgtemp.renderer        = 'painters';
% ft_singleplotTFR(cfgtemp, TFR_macro);
%
% colormap_axis = caxis;
% %caxis([0,colormap_axis(2)/2]);
% title('Mean frequency power of all EMG responses','Fontsize',11);
% xlabel('Time from SlowWave (s)','Fontsize',15);
% ylabel('Hz');
% set(gca,'FontWeight','bold' );
% set(gca,'TickDir','out');
% %axis tight
% xlim([-abscisse_scale, abscisse_scale]);


% Avg of EMG powspctrm

% power(i) = squeeze(TFR_macro.powspctrm(1,1,1,:)));




for i = 1:length(TFR_macro.time)
    TFR_macro.powspctrm(1,1,1,i) = average(TFR_macro.powspctrm(1,1,:,i));
end

yyaxis right

plot(TFR_macro.time,squeeze(TFR_macro.powspctrm(1,1,1,:)),'b','LineWidth', 2);

title('EMG : power average','Fontsize',18);
ylabel('EMG (pow avg)','Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
%axis tight
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);

%save

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

plot(TFR_macro.time,squeeze(TFR_macro.powspctrm(1,1,1,:)),'b','LineWidth', 2);
xlim([-abscisse_scale, abscisse_scale]);
%ylim([-4000 7000]); %attention adapté pour les diapos mais pas pour généraliser
ylim_emg = get(gca,'ylim');
ylim_emg(1)=ylim_emg(2)/ylim_rapport; %to set automatically EMG scale with same 0 as eeg
ylim(ylim_emg);

title('EEG-EMG comparison','Fontsize',18);
ylabel('EMG (pow avg)', 'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
%axis tight
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);

xlabel('Time (s)','Fontsize',18);


%% sava data
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,'comparison_eegemg',[cfg.prefix,'comparisoneegemg_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG},'_',data{imarker}.label{iEMG},'.pdf']),'-r600');
print(fig, '-dpng', fullfile(cfg.imagesavedir,'comparison_eegemg',[cfg.prefix,'comparisoneegemg_',cfg.LFP.name{imarker},'_',data{imarker}.label{iEEG},'_',data{imarker}.label{iEMG},'.png']),'-r600');
close all


end

