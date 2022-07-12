function PET_plot_LFP_aligned(cfg, ipart, markername)

MuseStruct      = readMuseMarkers(cfg, true);
M.non_aligned   = MuseStruct;
M.aligned_xcorr = alignMuseMarkersXcorr(cfg, MuseStruct, false);

close all
fig1 = figure;
sgtitle(sprintf('%s : %s', cfg.prefix(1:end-1), markername), 'interpreter', 'none');
fig2 = figure;
sgtitle(sprintf('%s : %s', cfg.prefix(1:end-1), markername), 'interpreter', 'none');
iplot = 0;

for alignment = ["non_aligned", "aligned_xcorr"]
    
    iplot = iplot + 1;
    cfgtemp = cfg;
    cfgtemp.LFP.write = false; %do not write LFP to disk because it is read several times, it would be confusing which one is left on the disk
    cfgtemp.LFP.name = {markername};
    LFP_temp = readLFP(cfgtemp, M.(alignment), true);
    LFP_temp = remove_artefacted_trials(cfg, LFP_temp);
    
    cfgtemp         = [];
    cfgtemp.channel = cfg.alignpeak.channel.(markername);
    data            = ft_selectdata(cfgtemp, LFP_temp{ipart}.(markername));
    dataavg         = ft_timelockanalysis([], data);
    
    %raw LFP
    figure(1);
    subplot(2,1,iplot); hold on;
    for itrial = 1:size(data.trial, 2)
        p = plot(data.time{itrial}, data.trial{itrial}, 'k');
        p.Color(4) = 0.2;
    end
    plot(dataavg.time, dataavg.avg, 'color', 'r', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel(cfg.alignpeak.channel.(markername), 'interpreter', 'none');
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    title(sprintf('%s (%d trials)', alignment, size(LFP_temp{ipart}.(markername).trial, 2)), 'interpreter', 'none');
    xlim(cfg.epoch.toi.(markername));
    
    % average LFP
    figure(2);
    subplot(2,1,iplot); hold on;
    patch_std(dataavg.time, dataavg.avg, sqrt(dataavg.var), 'k');
    plot(dataavg.time, dataavg.avg, 'k', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel(cfg.alignpeak.channel.(markername), 'interpreter', 'none');
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    title(sprintf('%s (%d trials)', alignment, size(LFP_temp{ipart}.(markername).trial, 2)), 'interpreter', 'none');
    xlim(cfg.epoch.toi.(markername))
end

%save the 2 figures:
figname1 = fullfile(cfg.imagesavedir, 'LFP', sprintf('%sp%d_LFP_%s_raw', cfg.prefix, ipart, markername));
savefigure_own(fig1, figname1, 'png', 'pdf', 'close');
figname2 = fullfile(cfg.imagesavedir, 'LFP', sprintf('%sp%dLFP_%s_avg', cfg.prefix, ipart, markername));
savefigure_own(fig2, figname2, 'png', 'pdf', 'close');

