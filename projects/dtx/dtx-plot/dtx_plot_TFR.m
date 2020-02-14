function dtx_plot_TFR(cfg,data,iEEG,imarker,saveplot)
%TFR of all trials of one electrode, with time of interest defined in cfg

data = data{imarker};
iEEG = iEEG{imarker};


cfgtemp                         = [];
cfgtemp.channel                 = iEEG;%'all'; 
cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
cfgtemp.output                  = 'pow';
cfgtemp.taper                   = 'hanning';
cfgtemp.pad                     = 'nextpow2';
cfgtemp.keeptrials              = 'yes'; %return average
cfgtemp.foi                     = 0:0.2:45;%80:0.2:125;
cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
%cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;

cfgtemp.toi                     = [-10:0.01:15];
TFR_macro                       = ft_freqanalysis(cfgtemp,data);

TFR_macro_log = TFR_macro;
TFR_macro_log.powspctrm = log(TFR_macro.powspctrm);

fig = figure;
subplot(2,1,1);

cfgtemp                 = [];
cfgtemp.baseline        = [-10, -5];
cfgtemp.baselinetype    = 'relchange';
cfgtemp.colorbar        = 'yes';
cfgtemp.zlim            = 'maxabs';
cfgtemp.xlim            = cfg.epoch.toi{imarker};
%cfgtemp.ylim            = [0 foi_max];
cfgtemp.parameter       = 'powspctrm';
cfgtemp.colormap        = parula(5000);
cfgtemp.renderer        = 'painters';

ft_singleplotTFR(cfgtemp, TFR_macro);

title(sprintf('%s%s : Frequency power over time',cfg.prefix,data.label{iEEG}),'Interpreter','none');
%xlim([-1, 1]);
%ylim([80 125]);
xlabel(sprintf('Time from %s (s)',cfg.LFP.name{imarker}));
ylabel('Frequency (Hz)');
% colormap_axis = caxis;
% caxis([0,colormap_axis(2)]);

subplot(2,1,2);
cfgtemp                 = [];
cfgtemp.baselinetype    = 'relchange';
cfgtemp.colorbar        = 'yes';
cfgtemp.zlim            = 'maxabs';
cfgtemp.xlim            = cfg.epoch.toi{imarker};
%cfgtemp.ylim            = [0 foi_max];
cfgtemp.parameter       = 'powspctrm';
cfgtemp.colormap        = parula(5000);
cfgtemp.renderer        = 'painters';

ft_singleplotTFR(cfgtemp, TFR_macro_log);

title(sprintf('%s%s : Frequency log power over time',cfg.prefix,data.label{iEEG}),'Interpreter','none');
%xlim([-5, 5]);
xlabel(sprintf('Time from %s (s)',cfg.LFP.name{imarker}));
ylabel('Frequency (Hz)');
% colormap_axis = caxis;
% caxis([0,colormap_axis(2)]);
 


%% print to file
if saveplot
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'TFR_', cfg.LFP.name{imarker}, '_', data.label{iEEG}]),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'TFR_', cfg.LFP.name{imarker}, '_', data.label{iEEG}]),'-r600');
    close all
end

end

