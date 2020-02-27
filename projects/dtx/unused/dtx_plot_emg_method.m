function dtx_plot_emg_method(cfg,data,ipart,imarker,saveplot)
%iEEG et 1 are the indexes of the channels in data.
% if 1 = 'no' : subplot indicating absence of EMG
%if 1 = false : ignoring everything concerning EMG

%% Setting parameters
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

cfgtemp = [];
cfgtemp.channel = cfg.LFP.emg{imarker};
data = ft_selectdata(cfgtemp,data{ipart}{imarker});

abscisse_scale = 2;%s
envelope_method = 'rms';
envelope_parameter = 30;
timewindowlength=0.05;%s, for TFR

nb_trials = size(data.trial,2);


%h automatic setting :
for itrial = 1 : nb_trials
    h_temp_max = max(data.trial{itrial}(1,:));
    h_temp_min = min(data.trial{itrial}(1,:));
    h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
end
h = mean(h_temp_amplitude);


%% Plot timecourse EMG
fig=figure;
subplot(1,3,1);
hold;

plot([0 0],[h/2 (nb_trials+1)*h], 'r', 'Linewidth', 2);

for itrial = 1 : nb_trials
    plot(data.time{itrial},data.trial{itrial}(1,:)+(nb_trials+1)*h- itrial*h,'k'); %first on top
end

%xlabel(sprintf('Time from \n%s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('Timecourse of \n%s (%d trials)', data.label{1}, size(data.trial,2)),'Interpreter','none','Fontsize',15);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(tick : tick*10 : nb_trials*h);
yticklabels(nb_trials : -10 : 0);
set(gca,'TickDir','out');
axis tight
xlim([-abscisse_scale abscisse_scale]);



%% Plot timecourse rectified EMG
h = h/2; %because no more negative values

subplot(1,3,2);
hold;

plot([0 0],[h/2 (nb_trials+1)*h], 'r', 'Linewidth', 2);

for itrial = 1 : nb_trials
    plot(data.time{itrial},abs(data.trial{itrial}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
end

xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
%ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('Timecourse of \nrectified %s (%d trials)', data.label{1}, size(data.trial,2)),'Interpreter','none','Fontsize',15);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(tick : tick*10 : nb_trials*h);
yticklabels(nb_trials : -10 : 0);
set(gca,'TickDir','out');
axis tight
xlim([-abscisse_scale abscisse_scale]);

%% Plot timecourse of envelope of rectified EMG
subplot(1,3,3);
hold;

plot([0 0],[h/2 (nb_trials+1)*h], 'r', 'Linewidth', 2);

env_upper = [];
for itrial = 1 : nb_trials
    [env_upper{itrial}, ~] = envelope(abs(data.trial{itrial}(1,:)),envelope_parameter,envelope_method);
    plot(data.time{itrial},abs(data.trial{itrial}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
    plot(data.time{itrial},env_upper{itrial}+ (nb_trials+1)*h - itrial*h,'c','LineWidth',2);
end

%xlabel(sprintf('Time from \n%s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
%ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('Timecourse of envelope \nof rectified %s (%d trials)', data.label{1}, size(data.trial,2)),'Interpreter','none','Fontsize',15);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
tick = h;
yticks(tick : tick*10 : nb_trials*h);
yticklabels(nb_trials : -10 : 0);
set(gca,'TickDir','out');
axis tight
xlim([-abscisse_scale abscisse_scale]);


%% save data
if saveplot
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('Create folder %s',cfg.imagesavedir);
    end
    
    if ~(exist (fullfile(cfg.imagesavedir, 'emg_method'))==7)
        mkdir(fullfile(cfg.imagesavedir, 'emg_method'));
        warning('Create folder %s',fullfile(cfg.imagesavedir, 'emg_method'));
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_emg_method_',data.label{1}]),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_emg_method_',data.label{1}]),'-r600');
    
    close all
    
end


% %% TFR
% fig2=figure;
% fig2.Renderer    = 'Painters'; % Else pdf is saved to bitmap
% subplot(4,1,1)
% 
% 
% 
% % TFR EMG
% %yyaxis left
% cfgtemp                         = [];
% cfgtemp.channel                 = 1;%'all'; %ichannel;
% cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
% cfgtemp.output                  = 'pow';
% cfgtemp.taper                   = 'hanning';
% cfgtemp.pad                     = 'nextpow2';
% cfgtemp.keeptrials              = 'yes'; %return average
% cfgtemp.foi                     = 0:1:125;
% cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*timewindowlength;%40./cfgtemp.foi;
% %cfgtemp.t_ftimwin               = 20./cfgtemp.foi;
% %cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
% %cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;
% 
% cfgtemp.toi                     = [-abscisse_scale:0.01:abscisse_scale];
% TFR_macro                       = ft_freqanalysis(cfgtemp,data);
% 
% cfgtemp = [];
% cfgtemp.colorbar        = 'no';
% %cfgtemp.baseline        = [-4, -2];
% % cfgtemp.baselinetype    = 'relchange';
% % cfgtemp.zlim            = 'maxabs';
% cfgtemp.parameter       = 'powspctrm';
% cfgtemp.colormap        = parula(5000);
% cfgtemp.renderer        = 'painters';
% ft_singleplotTFR(cfgtemp, TFR_macro);
% 
% title(sprintf('%s : TFR and average power',data.label{1}),'Interpreter','none','Fontsize',18);
% 
% 
% 
% for i = 1:length(TFR_macro.time)
%     TFR_macro.powspctrm(1,1,1,i) = average(TFR_macro.powspctrm(1,1,:,i));
% end
% 
% yyaxis right
% hold
% 
% plot(TFR_macro.time,squeeze(TFR_macro.powspctrm(1,1,1,:)),'r','LineWidth', 2);
% 
% 
% 
% %% average of hilbert of abs
% 
% %plot([0 0],[0 n*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
% subplot(4,1,2)
% title(sprintf('Average of envelopes of abs of %s',data.label{1}),'Interpreter','none','Fontsize',18);
% 
% hold
% 
% env_upper = [];
% for itrial = 1 : n
%     t = data.time{itrial};
%     [env_upper{itrial} ~] = envelope(abs(data.trial{itrial}(1,:)),envelope_parameter,envelope_method);
%     plot(t,abs(data.trial{itrial}(1,:)),'color',[0.6 0.6 0.6]); %first on top
%     plot(t,env_upper{itrial},'b');
% end
% 
% for ioffset = 1:length(env_upper{1})
%     for itrial = 1:length(env_upper)
%         envupper_by_offset(itrial) = env_upper{itrial}(ioffset);
%     end
%     envupper_avg(ioffset) = mean(envupper_by_offset);
% end
% 
% 
% plot(data.time{1},envupper_avg,'r','LineWidth',2);
% xlim([-abscisse_scale abscisse_scale]);
% 
% 
% %% hilbert of average of abs
% subplot(4,1,3)
% title(sprintf('Envelope of average of abs of %s',data.label{1}),'Interpreter','none','Fontsize',18);
% 
% hold
% 
% data_abs                = data;
% data_abs.label          = data.label(1);
% for itrial = 1:length(data.trial)
%     data_abs.trial{itrial}                = abs(data.trial{itrial}(1,:));
% end
% cfgtemp                 = [];
% data_abs_avg            = ft_timelockanalysis(cfgtemp,data_abs);
% 
% plot(data_abs_avg.time,data_abs_avg.avg,'k');
% plot(data_abs_avg.time,envelope(data_abs_avg.avg,envelope_parameter,envelope_method),'r','LineWidth', 2);
% xlim([-abscisse_scale abscisse_scale]);
% 
% 
% %% hilbert of average of rms
% subplot(4,1,4)
% title(sprintf('Envelope of average of rms of %s',data.label{1}),'Interpreter','none','Fontsize',18);
% hold
% 
% data_by_offset=[];
% for ioffset = 1:length(data.time{1})
%     for itrial = 1:length(data.trial)
%         data_by_offset(itrial) = data.trial{itrial}(1,ioffset);
%     end
%     data_rms(ioffset) = rms(data_by_offset);
% end
% 
% 
% plot(data.time{1},data_rms,'k');
% plot(data.time{1},envelope(data_rms,envelope_parameter,envelope_method),'r','LineWidth', 2);
% xlim([-abscisse_scale abscisse_scale]);
% xlabel(sprintf('Time from \n%s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
% 
% 
% %% overdraw of all methods
% fig3 = figure;
% hold
% title(sprintf('%s : comparison of quantification methods (normalized)',data.label{1}),'Interpreter','none','Fontsize',18);
% 
% %normaliser
% envupper_avg_norm = envupper_avg/max(envupper_avg);
% env_avgabs_norm = envelope(data_abs_avg.avg,envelope_parameter,envelope_method)/max(envelope(data_abs_avg.avg,envelope_parameter,envelope_method));
% env_avgrms_norm = envelope(data_rms,envelope_parameter,envelope_method)/max(envelope(data_rms,envelope_parameter,envelope_method));
% TFR_norm = squeeze(TFR_macro.powspctrm(1,1,1,:))/max(squeeze(TFR_macro.powspctrm(1,1,1,:)));
% %TFR_env_norm = envelope(squeeze(TFR_macro.powspctrm(1,1,1,:)),envelope_parameter,envelope_method)/max(envelope(squeeze(TFR_macro.powspctrm(1,1,1,:)),envelope_parameter,envelope_method));
% 
% t=data.time{1}; %all trials same size
% plot(t,envupper_avg_norm,t,env_avgabs_norm,t,env_avgrms_norm,'LineWidth',2);
% plot(TFR_macro.time,TFR_norm,'LineWidth',2);
% %plot(TFR_macro.time,TFR_env_norm);
% xlim([-abscisse_scale abscisse_scale]);
% legend('Avg of env of abs', 'Env of avg of abs', 'Env of avg of rms', 'TFR avg power');%,'Env of TFR avg power');
% 
% 
% xlabel(sprintf('Time from \n%s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
% 
% if saveplot
%     %% print to file
%     set(fig1,'PaperOrientation','landscape');
%     set(fig1,'PaperUnits','normalized');
%     set(fig1,'PaperPosition', [0 0 1 1]);
%     set(fig1,'Renderer','Painters');
%     print(fig1, '-dpdf', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_env_of_abs']),'-r600');
%     print(fig1, '-dpng', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_env_of_abs']),'-r600');
% 
%     set(fig2,'PaperOrientation','landscape');
%     set(fig2,'PaperUnits','normalized');
%     set(fig2,'PaperPosition', [0 0 1 1]);
%     set(fig2,'Renderer','Painters');
%     print(fig2, '-dpdf', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_all']),'-r600');
%     print(fig2, '-dpng', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_all']),'-r600');
% 
%     set(fig3,'PaperOrientation','landscape');
%     set(fig3,'PaperUnits','normalized');
%     set(fig3,'PaperPosition', [0 0 1 1]);
%     set(fig3,'Renderer','Painters');
%     print(fig3, '-dpdf', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_comparison']),'-r600');
%     print(fig3, '-dpng', fullfile(cfg.imagesavedir,'emg_method',[cfg.prefix,cfg.LFP.name{imarker},'_',data.label{1},'_emgmethod_comparison']),'-r600');
%  
%     close all
% 
% end

end

