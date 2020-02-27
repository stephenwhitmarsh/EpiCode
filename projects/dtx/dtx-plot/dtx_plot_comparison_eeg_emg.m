function [data_EEG, data_EMG] = dtx_plot_comparison_eeg_emg(cfg,data,ipart,imarker,saveplot)
abscisse_scale = 2;

data = data{ipart};

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


%select the channels of interest
cfgtemp = [];
cfgtemp.channel = cfg.align.channel{imarker};
data_EEG = ft_selectdata(cfgtemp,data{imarker});


data_EMG = [];
if isfield(cfg.LFP, 'emg')
    if ~strcmp(cfg.LFP.emg{imarker},'no')
        if ~(cfg.LFP.emg{imarker} == false)
            cfgtemp = [];
            cfgtemp.channel = cfg.LFP.emg{imarker};
            data_EMG = ft_selectdata(cfgtemp,data{imarker});
        else
            warning('cfg.LFP.emg{imarker} == false : no EMG data loaded, plot not made');
            return
        end
    else
        warning('cfg.LFP.emg{imarker} == ''no'' : no EMG data loaded, plot not made');
        return
    end
else
    warning('Field cfg.LFP.emg does not exist : no EMG data loaded, plot not made');
    return
end


fig = figure;

%% Overdraw OL
subplot(3,1,1);
hold;

%find good Y lim to avoid flattening by too much noise
for itrial = 1 : length(data_EEG.trial)
t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data_EEG.fsample; % offset for which t = 0;
y_max_temp(itrial) = max(data_EEG.trial{itrial}(1,round(-0.5*data_EEG.fsample)+t_0: round(0.5*data_EEG.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
y_min_temp(itrial) = min(data_EEG.trial{itrial}(1,round(-0.5*data_EEG.fsample)+t_0: round(0.5*data_EEG.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
end
y_max = max(y_max_temp);
y_min = min(y_min_temp);


fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
for itrial = 1 : size(data_EEG.trial,2)
    plot(data_EEG.time{itrial},data_EEG.trial{itrial}(1,:),'color',[0.6 0.6 0.6]); %first on top
end

title(sprintf('%s - EEG/EMG comparison : %d trials \n\n%s : EEG average',cfg.prefix(1:end-1),length(data{imarker}.trial),cfg.LFP.name{imarker}(1:end-4)),'Fontsize',18,'Interpreter','none');
ylabel(sprintf('EEG %s (µV)',data_EEG.label{1}),'Fontsize',15,'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
axis tight;
xlim([-abscisse_scale, abscisse_scale]);
%ylim([y_min y_max]);
set(gca,'ycolor','r');
set(gca,'Fontsize',15);


%% Average OL

cfgtemp                 = [];
cfgtemp.vartrllength    = 2;
data_EEG_rptavg        = ft_timelockanalysis(cfgtemp,data_EEG);

plot(data_EEG_rptavg.time,data_EEG_rptavg.avg(1,:),'r','LineWidth', 2);


%% Overdraw and rectified EMG
subplot(3,1,2);
hold;


title(sprintf('Average of envelopes of abs of %s',data_EMG.label{1}),'Interpreter','none','Fontsize',18);

env = [];

%plot trial by trial : rect and enveloppe

for itrial = 1 : size(data_EMG.trial,2)
    t = data{imarker}.time{itrial};
    rect_emg = abs(data_EMG.trial{itrial}(1,:));
    plot(t,rect_emg,'color',[0.6 0.6 0.6]); %first on top
end

for itrial = 1 : size(data_EMG.trial,2)
    t = data{imarker}.time{itrial};
    rect_emg = abs(data_EMG.trial{itrial}(1,:));
    [env{itrial}, ~] = envelope(rect_emg,20,'rms');
    plot(t,env{itrial},'color','c');
end

%plot avg of enveloppe
for ioffset = 1:length(env{1}) %all trials must have the same size
    for itrial = 1:length(env)
        env_by_offset(itrial) = env{itrial}(ioffset);
    end
    env_avg(ioffset) = mean(env_by_offset);
end
plot(t,env_avg,'b','LineWidth',2);


ylabel(sprintf('%s (µV)',data_EMG.label{1}),'Fontsize',15,'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');



% %  Old method : avg of abs
% data_EMG_abs                = data_EMG;
% data_EMG_abs.label          = data_EMG.label(1);
% for itrial = 1:length(data_EMG.trial)
%     data_EMG_abs.trial{itrial}                = abs(data_EMG.trial{itrial}(1,:));
% end
% cfgtemp                 = [];
% data_EMG_abs_avg            = ft_timelockanalysis(cfgtemp,data_EMG_abs);
% 
% plot(data_EMG_abs_avg.time,data_EMG_abs_avg.avg,'b','LineWidth', 2);
%plot(data_EMG_abs_avg.time,envelope(data_EMG_abs_avg.avg,15,'rms'),'-b','LineWidth', 2);



title('Average enveloppe of rectified EMG','Fontsize',18);
ylabel(sprintf('%s (µV)',data_EMG.label{1}),'Fontsize',15,'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);



%% Comparaison EEG-EMG
subplot (3,1,3)
hold;

plot(data_EEG_rptavg.time,data_EEG_rptavg.avg(1,:),'r','LineWidth', 2);

yyaxis left
set(gca,'ycolor','r');
ylabel(sprintf('EEG %s (µV)',data_EEG.label{1}), 'Fontsize',15);

axis tight
xlim([-abscisse_scale, abscisse_scale]);
ylim_eeg = get(gca,'ylim');
ylim_rapport = -ylim_eeg(1)/ylim_eeg(2); %to set automatically EMG scale with same 0 as eeg

yyaxis right

%plot(data_EMG_abs_avg.time,envelope(data_EMG_abs_avg.avg,15,'rms'),'b','LineWidth', 2);
plot(t,env_avg,'b','LineWidth', 2);

%set lower y of emg equal to y=0 eeg
ylim_emg = get(gca,'ylim');
ylim_emg_rapport(1)=ylim_emg(1)-ylim_emg(2)*ylim_rapport; %to set automatically EMG scale with same 0 as EEG
ylim_emg_rapport(2)=ylim_emg(2);
ylim(ylim_emg_rapport);


title('EEG-EMG comparison','Fontsize',18);
ylabel(sprintf('%s (avg env)',data_EMG.label{1}), 'Fontsize',15);
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');

set(gca,'ycolor','b');
set(gca,'Fontsize',15);

xlabel('Time (s)','Fontsize',18);


%% sava data
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprinf('Create forlder %s',cfg.imagesavedir);
    end
    
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_',data_EEG.label{1},'_',data_EMG.label{1},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_',data_EEG.label{1},'_',data_EMG.label{1},'.png']),'-r600');
    close all
end

%% rename for saving over several patients
data_EEG.label{1} = 'chan_SlowWave';
data_EEG.ID = cfg.prefix(1:end-1);
data_EMG.label{1} = 'EMG';
data_EMG.ID = cfg.prefix(1:end-1);

end

