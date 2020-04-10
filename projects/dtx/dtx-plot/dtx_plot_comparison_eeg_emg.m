function envEMG = dtx_plot_comparison_eeg_emg(cfg,data,ipart,imarker,saveplot,varargin)

abscisse_scale = 2;

data = data{ipart};

data_EEG = [];
data_EMG = [];

%get EEG and EMG channels
if isempty(varargin)
    chanEEG = cfg.align.channel{imarker};
    if isfield(cfg.LFP, 'emg')
        if ~strcmp(cfg.LFP.emg{imarker},'no')
            if ~(cfg.LFP.emg{imarker} == false)
                chanEMG = cfg.LFP.emg{imarker};
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
elseif length(varargin) == 2
    chanEEG = varargin{1};
    chanEMG = varargin{2};
else
    error('varargin must have 0 argument, or 2 arguments (chan name chanEEG and chanEMG).\nHere there are %d arguments',length(varargin));
end

%select the channels of interest
cfgtemp = [];
cfgtemp.channel = chanEEG;
data_EEG = ft_selectdata(cfgtemp,data{imarker});

cfgtemp = [];
cfgtemp.channel = chanEMG;
data_EMG = ft_selectdata(cfgtemp,data{imarker});

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

%plot trial by trial : rect and envelope

for itrial = 1 : size(data_EMG.trial,2)
    t = data{imarker}.time{itrial};
    rect_emg = abs(data_EMG.trial{itrial}(1,:));
    plot(t,rect_emg,'color',[0.6 0.6 0.6]); %first on top
end

for itrial = 1 : size(data_EMG.trial,2)
    t = data{imarker}.time{itrial};
    rect_emg = abs(data_EMG.trial{itrial}(1,:));
    [env{itrial}, ~] = envelope(rect_emg,cfg.EMG.envparam,cfg.EMG.envmethod);
    plot(t,env{itrial},'color','c');
end

%plot avg of envelope
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

%Temporaire Paul pour sortir la valeur
envEMG{ipart}{imarker}.time{1} = t;
envEMG{ipart}{imarker}.trial{1} = env_avg; 
envEMG{ipart}{imarker}.label{1} = 'EMG';

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



title('Average envelope of rectified EMG','Fontsize',18);
ylabel(sprintf('%s (µV)',data_EMG.label{1}),'Fontsize',15,'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlim([-abscisse_scale, abscisse_scale]);
set(gca,'ycolor','b');
set(gca,'Fontsize',15);



%% Comparaison EEG-EMG
subplot (3,1,3)
hold;

% plot(data_EEG_rptavg.time,data_EEG_rptavg.avg(1,:),'r','LineWidth', 2);
% %A REMETTRE

yyaxis left
set(gca,'ycolor','r');
ylabel(sprintf('EEG %s (µV)',data_EEG.label{1}), 'Fontsize',15);

axis tight
xlim([-abscisse_scale, abscisse_scale]);
%ylim_eeg = get(gca,'ylim');
%ylim_rapport = -ylim_eeg(1)/ylim_eeg(2); %to set automatically EMG scale with same 0 as eeg

yyaxis right

%plot(data_EMG_abs_avg.time,envelope(data_EMG_abs_avg.avg,15,'rms'),'b','LineWidth', 2);
plot(t,env_avg,'b','LineWidth', 2);

%set lower y of emg equal to y=0 eeg
% ylim_emg = get(gca,'ylim');
% ylim_emg_rapport(1)=ylim_emg(1)-ylim_emg(2)*ylim_rapport; %to set automatically EMG scale with same 0 as EEG
% ylim_emg_rapport(2)=ylim_emg(2);
% ylim(ylim_emg_rapport);


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
    
    
    set(fig,'PaperOrientation','portrait');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,['POURDIAPO',cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_',data_EEG.label{1},'_',data_EMG.label{1},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,['POURDIAPO',cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_',data_EEG.label{1},'_',data_EMG.label{1},'.png']),'-r600');
    close all
end


end

