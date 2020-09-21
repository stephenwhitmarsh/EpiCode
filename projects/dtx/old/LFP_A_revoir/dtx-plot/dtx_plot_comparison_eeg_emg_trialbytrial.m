function dtx_plot_comparison_eeg_emg_trialbytrial(cfg,data,ipart,imarker,saveplot,varargin)
%one fig of comparison per trial. Made for plotting avg data patient by
%patient

abscisse_scale = 2;

data = data{ipart}{imarker};

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
data_EEG = ft_selectdata(cfgtemp,data);

cfgtemp = [];
cfgtemp.channel = chanEMG;
data_EMG = ft_selectdata(cfgtemp,data);

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
for itrial = 1:length(data.trial)
    
    subplot(2,round(length(data.trial)/2),itrial)
    hold;
    
    t = data_EEG.time{1};
    
    %plot eeg
    yyaxis left
    plot(t, data_EEG.trial{itrial},'r','LineWidth',2);
    
    set(gca,'ycolor','r');
    ylabel(sprintf('EEG %s avg', chanEEG), 'Fontsize',15);
    axis tight
    xlim([-abscisse_scale, abscisse_scale]);
%     ylim_eeg = get(gca,'ylim');
%     ylim_rapport = -ylim_eeg(1)/ylim_eeg(2); %to set automatically EMG scale with same 0 as eeg
%     
    %plot emg
    yyaxis right
    plot(t, data_EMG.trial{itrial},'b','LineWidth',2);
%     ylim_emg = get(gca,'ylim');
%     ylim_emg_rapport(1)=ylim_emg(1)-ylim_emg(2)*ylim_rapport; %to set automatically EMG scale with same 0 as EEG
%     ylim_emg_rapport(2)=ylim_emg(2);
%     ylim(ylim_emg_rapport);
%     ylabel('EMG avg envelope', 'Fontsize',15);

    set(gca,'FontWeight','bold' );
    set(gca,'TickDir','out');
    
    set(gca,'ycolor','b');
    set(gca,'Fontsize',15);
    
    xlabel('Time (s)','Fontsize',18);
    xlim([-abscisse_scale, abscisse_scale]);

    title(sprintf('Patient n°%d', itrial),'Fontsize',18);

end

if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprintf('Create forlder %s',cfg.imagesavedir);
    end
    
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_trial_by_trial_',data_EEG.label{1},'_',data_EMG.label{1},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_comparisoneegemg_trial_by_trial_',data_EEG.label{1},'_',data_EMG.label{1},'.png']),'-r600');
    close all
end


end

