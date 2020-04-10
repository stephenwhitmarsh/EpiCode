function dtx_plot_timecourse_eeg_emg(cfg,data,ipart,imarker,datatype,electrodeToPlot,toi,saveplot)
%data{ipart} 
%This script is for one ipart and one imarker
%datatype : 'eeg' or 'emg'

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


data = data{ipart};

isEEG = 0;
isEMG = 0;
if strcmp(datatype, 'eeg')
    isEEG = 1;
elseif strcmp(datatype, 'emg')
    isEMG = 1;
else
    error('Error in the function arguments : datatype must be ''eeg'' or ''emg''');
end




%select the channel of interest
if isEEG
    cfgtemp = [];
    cfgtemp.channel = electrodeToPlot;
    data = ft_selectdata(cfgtemp,data{imarker});
   
    
elseif isEMG
    if isfield(cfg.LFP, 'emg')
        if ~strcmp(cfg.LFP.emg{imarker},'no')
            if ~(cfg.LFP.emg{imarker} == false)
                cfgtemp = [];
                cfgtemp.channel = electrodeToPlot;
                data = ft_selectdata(cfgtemp,data{imarker});
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
    
end

nb_trials = size(data.trial,2);

%% EEG
fig1 = figure;
hold;

%h automatic setting :
for itrial = 1 : nb_trials
    if isEMG %max around t0
        t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data.fsample; % offset for which t = 0;
        h_temp(itrial) = max(data.trial{itrial}(1,round(-0.5*data.fsample)+t_0: round(0.5*data.fsample)+t_0)); %amplitude of peak vs baseline
    elseif isEEG %amplitude t0-baseline
        bl_period_inf = data.time{itrial}>cfg.align.toibaseline{imarker}(1);
        bl_period_sup = data.time{itrial}<cfg.align.toibaseline{imarker}(2);
        bl_period = logical(bl_period_inf .* bl_period_sup);
        h_temp(itrial) = max(data.trial{itrial}(1, data.time{itrial}>-1 & data.time{itrial}<1)) - nanmean(data.trial{itrial}(1, bl_period));
    end
end

h = mean(h_temp)*2;

% if isEEG && isMicromed
%     h = mean(h_temp)/2;
% elseif isEEG && isBrainvision
%     h = mean(h_temp)*2;
% elseif isEMG
%     h = mean(h_temp)*2;
% end

%plot([0 0],[0 nb_trials*h+h], 'r', 'Linewidth', 2); %MODIF TEMPORAIRE PAUL
%plot([0 0],[0 n*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
for itrial = 50:75%1 : nb_trials 
    plot(data.time{itrial},data.trial{itrial}(1,:)+ (nb_trials+1)*h - itrial*h,'k'); %first on top
end


xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
ylabel('Number of seizures', 'Fontsize',15);
title(sprintf('%s : aligned data from %s (%d trials)', cfg.LFP.name{imarker}, data.label{1}, size(data.trial,2)),'Interpreter','none','Fontsize',20);
set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);

tick_interval = round(nb_trials/10);
if tick_interval == 0
    tick_interval = 1;
end
yticks(h : h*tick_interval : nb_trials*h);
yticklabels(nb_trials : -tick_interval : 0);

set(gca,'TickDir','out');
axis tight

xlim(toi);


% %% EMG
% 
% if isfield(cfg.LFP, 'emg')
%     if ~strcmp(cfg.LFP.emg{imarker},'no')
%         if ~(cfg.LFP.emg{imarker} == false)
%             
%             fig2 = figure;
%             hold;
%             
%             %h automatic setting :
%             for itrial = 1 : nb_trials
%                 h_temp_max = max(data_EMG.trial{itrial}(1,:));
%                 h_temp_min = min(data_EMG.trial{itrial}(1,:));
%                 h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
%             end
%             h = mean(h_temp_amplitude);
%             
%             
%             
%             %plot([0 0],[0 n*h+h], 'color',[0.6 0.6 0.6], 'Linewidth', 2);
%             plot([0 0],[0 nb_trials*h+h], 'r', 'Linewidth', 2);
%             
%             for itrial = 1 : nb_trials
%                 t = data_EMG.time{itrial};
%                 plot(t,data_EMG.trial{itrial}(1,:)+ (nb_trials+1)*h - itrial*h,'k'); %first on top
%             end
%             
%             
%             
%             xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
%             ylabel('Number of seizures', 'Fontsize',15);
%             title(sprintf('Aligned data from %s', data_EMG{imarker}.label{iEMG{imarker}}),'Interpreter','none','Fontsize',18);
%             set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
%             tick = h;
%             yticks(tick : tick*10 : nb_trials*h);
%             yticklabels(nb_trials : -10 : 0);
%             set(gca,'TickDir','out');
%             axis tight
%             xlim(cfg.epoch.toi{1});
%         end
%     end
% end
% 
% if isfield(cfg.LFP, 'emg')
%     if ~strcmp(cfg.LFP.emg{imarker},'no')
%         subplot(1,2,2);
%         set(gca,'TickLength',[0 0]);
%         yticklabels([]);
%         xticklabels([]);
%         set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
%         text(0.2,0.5,sprintf('No EMG data \nassociated with %s', cfg.LFP.name{imarker}),'Interpreter','none','Fontsize',15);
%     end
% end


%% print to file
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('Create folder %s',cfg.imagesavedir);
    end
        
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition', [0 0 1 1]);
    set(fig1,'Renderer','Painters');
     
    print(fig1, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_timecourse_',data.label{1},'_scale[',num2str(toi),']']),'-r600');
    print(fig1, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_timecourse_',data.label{1},'_scale[',num2str(toi),']']),'-r600');
    

%     %EMG
%     if isfield(cfg.LFP, 'emg')
%         if ~strcmp(cfg.LFP.emg{imarker},'no')
%             if ~(cfg.LFP.emg{imarker} == false)
%                 set(fig2,'PaperOrientation','landscape');
%                 set(fig2,'PaperUnits','normalized');
%                 set(fig2,'PaperPosition', [0 0 1 1]);
%                 set(fig2,'Renderer','Painters');
%                 
%                 print(fig2, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_EMG_',data_EMG.label{1}]),'-r600');
%                 print(fig2, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'timecourse_',cfg.LFP.name{imarker},'_EMG_',data_EMG.label{1}]),'-r600');
%             end
%         end
%     end
    
    close all
    
end

end



