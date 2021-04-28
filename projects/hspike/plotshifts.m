fig = figure;
set(fig,'position', get(0,'ScreenSize'));
set(fig, 'PaperPositionMode', 'auto');

markername = "combined1";

ncol = max(SpikeTrials_timelocked{ipatient}{1}.(markername).trialinfo.idir);

for idir = unique(SpikeTrials_timelocked{ipatient}{1}.(markername).trialinfo.idir)'
    
    cfgtemp = [];
    cfgtemp.trials = LFP{ipatient}{1}.(markername).trialinfo.idir == idir;
    cfgtemp.avgoverrpt = 'yes';
    sel = ft_selectdata(cfgtemp, LFP{ipatient}{1}.(markername));
    subplot(2, ncol, idir);
    plot(sel.time{1}, sel.trial{1}'); 
    hold on
    ax = axis;
    xlim([-0.5, 0.5]);
    plot([0, 0], [ax(3), ax(4)], ':k');
    
    cfg_raster                  = [];
    % cfgtemp.latency       = [cfg.stats.toi.(markername)(1), cfg.stats.toi.(markername)(2)];
    cfg_raster.trialborders     = 'yes';
    cfg_raster.linewidth        = 1;
    cfg_raster.topplotsize      = 0.5;
    cfg_raster.trialborders     = 'no';
    cfg_raster.trials           = SpikeTrials_timelocked{ipatient}{1}.(markername).trialinfo.idir == idir;
    
    subplot(2, ncol, idir + ncol);
    
    ft_spike_plot_raster(cfg_raster, SpikeTrials_timelocked{ipatient}{1}.(markername));
    hold on
    plot([0, 0], [ax(3), ax(4)], ':k');
    
    ax = axis;
    xlim([-0.5, 0.5]);
    plot([0, 0], [ax(3), ax(4)], ':k');
    
end
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
fname = fullfile(config{ipatient}.imagesavedir, [config{ipatient}.prefix, '_shift_per_directory']);
%     exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff'], 'Resolution', '600');

%%%
fig = figure;
set(fig,'position', get(0,'ScreenSize'));
set(fig, 'PaperPositionMode', 'auto');

cfgtemp = [];
cfgtemp.avgoverrpt = 'yes';
sel = ft_selectdata(cfgtemp, LFP{ipatient}{1}.(markername));
subplot(2, 2, 1);
plot(sel.time{1}, sel.trial{1}');
hold on
ax = axis;
xlim([-0.5, 0.5]);
plot([0, 0], [ax(3), ax(4)], ':k');

cfg_raster                  = [];
cfg_raster.trialborders     = 'yes';
cfg_raster.linewidth        = 1;
cfg_raster.topplotsize      = 0.5;
cfg_raster.trialborders     = 'no';
subplot(2, 2, 3);

ft_spike_plot_raster(cfg_raster, SpikeTrials_timelocked{ipatient}{1}.(markername));
hold on
plot([0, 0], [ax(3), ax(4)], ':k');

ax = axis;
xlim([-0.5, 0.5]);
plot([0, 0], [ax(3), ax(4)], ':k');

fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0.1 0.1 1 1]);
fname = fullfile(config{ipatient}.imagesavedir, [config{ipatient}.prefix, 'shift_all_directories']);
%     exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff'], 'Resolution', '600');