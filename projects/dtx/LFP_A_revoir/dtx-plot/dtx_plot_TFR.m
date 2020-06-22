function TFR_data = dtx_plot_TFR(cfg,data,ipart,imarker,dograndavg,saveplot)
%TFR of all trials of one electrode, with time of interest defined in cfg
%For eeg TFR analysis : data{ipart}{imarker}
%For TFR grand avg : data{imarker}{irat} 
%ipart is needed only for eeg TFR.

if dograndavg
    isdata = ~cellfun('isempty',data{imarker});
    if any(isdata)
        [TFR_data] = ft_freqgrandaverage([], data{imarker}{isdata});
    end
    
else
    data = data{ipart}{imarker};
    
    %Select EEG channel
    cfgtemp = [];
    cfgtemp.channel = cfg.LFP.electrodetoplot{imarker};
    data = ft_selectdata(cfgtemp,data);
    
    cfgtemp                         = [];
    cfgtemp.channel                 = 'all';
    cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
    cfgtemp.output                  = 'pow';
    cfgtemp.taper                   = 'hanning';
    cfgtemp.pad                     = 'nextpow2';
    cfgtemp.keeptrials              = 'no'; 
    cfgtemp.foi                     = 0:0.2:50;%80:0.2:125;
    cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
    %cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;
    cfgtemp.toi                     = 'all';
    TFR_data                        = ft_freqanalysis(cfgtemp,data);
end

TFR_macro_log = TFR_data;
TFR_macro_log.powspctrm = log(TFR_data.powspctrm);

fig = figure;
subplot(2,1,1);

cfgtemp                 = [];
cfgtemp.baseline        = cfg.LFP.TFR.baseline;
cfgtemp.baselinetype    = cfg.LFP.TFR.baselinetype;
cfgtemp.colorbar        = 'yes';
cfgtemp.zlim            = 'maxabs';
cfgtemp.xlim            = cfg.epoch.toi{imarker};
cfgtemp.ylim            = 'maxmin';
cfgtemp.parameter       = 'powspctrm';
cfgtemp.colormap        = parula(5000);
cfgtemp.renderer        = 'painters';

ft_singleplotTFR(cfgtemp, TFR_data);

title(sprintf('%s%s : Frequency power over time',cfg.prefix,TFR_data.label{1}),'Interpreter','none','Fontsize',20);
%xlim([-1, 1]);
%ylim([80 125]);
xlabel(sprintf('Time from %s (s)',cfg.LFP.name{imarker}),'Interpreter','none');
ylabel('Frequency (Hz)');
set(gca, 'FontWeight','bold', 'Fontsize',15);
colormap_axis = caxis;
caxis([0,colormap_axis(2)]);

subplot(2,1,2);
cfgtemp                 = [];
cfgtemp.baseline        = cfg.LFP.TFR.baseline;
cfgtemp.baselinetype    = cfg.LFP.TFR.baselinetype;
cfgtemp.colorbar        = 'yes';
cfgtemp.zlim            = 'maxabs';
cfgtemp.xlim            = cfg.epoch.toi{imarker};
%cfgtemp.ylim            = [0 foi_max];
cfgtemp.parameter       = 'powspctrm';
cfgtemp.colormap        = parula(5000);
cfgtemp.renderer        = 'painters';

ft_singleplotTFR(cfgtemp, TFR_macro_log);

title(sprintf('%s%s : Frequency log power over time',cfg.prefix,TFR_data.label{1}),'Interpreter','none','Fontsize',20);
%xlim([-5, 5]);
xlabel(sprintf('Time from %s (s)',cfg.LFP.name{imarker}),'Interpreter','none');
ylabel('Frequency (Hz)');
set(gca, 'FontWeight','bold', 'Fontsize',15);

colormap_axis = caxis;
caxis([0,colormap_axis(2)]);
 


%% print to file
if saveplot
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end

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

    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    
    if dograndavg
        method = 'TFRgrandavg_';
    else
        method = 'TFR_';
    end
    
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,method, cfg.LFP.name{imarker}, '_', TFR_data.label{1}]),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,method, cfg.LFP.name{imarker}, '_', TFR_data.label{1}]),'-r600');
    close all
end

end

