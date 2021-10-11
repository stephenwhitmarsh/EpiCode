function [cfg, LFP] = alignClusters(cfg, LFP)

cfg.orientation = ft_getopt(cfg, 'orientation', 'landscape');

t = [-0.15 0.10];

% t-zero LFPs
fig         = figure;
papersize   = 800;
set(fig, 'PaperPositionMode', 'auto');
if strcmp(cfg.orientation, 'landscape')
    set(fig, 'position', [20 -60 papersize*sqrt(2) papersize]);
else
    set(fig, 'position', [20 -60 papersize papersize*sqrt(2)]);
end
set(fig, 'Renderer', 'Painters');

for imarker = 1 : size(LFP, 2)
    
    if isempty(LFP{imarker})
        continue
    end
    
    cfg_temp = [];
    cfg_temp.latency = cfg.cluster.latency;
    LFP{imarker} = ft_selectdata(cfg_temp, LFP{imarker});
    
    chani       = find(contains(LFP{imarker}.label, cfg.align.zerochannel));
    subplot(2, size(LFP, 2), imarker); hold;
    
    plot(LFP{imarker}.time, LFP{imarker}.avg', 'k');
    plot(LFP{imarker}.time, LFP{imarker}.avg(chani, :)', 'r');
    timeindx    = LFP{imarker}.time > t(1) & LFP{imarker}.time < t(2);
    [~, LOCS]   = findpeaks(-LFP{imarker}.avg(chani, timeindx), 'SortStr', 'descend');
    if isempty(LOCS)
        continue
    end
    
    i0 = find(LFP{imarker}.time > t(1), 1, 'first') + LOCS(1) - 1;
    t0 = LFP{imarker}.time(i0);
    axis tight; ax = axis;
    plot([t0, t0], [1.2 * ax(3), 1.2 * ax(4)], ':r');
    title(sprintf('Shift = %.0fms', -t0*1000), 'FontSize', 6);
    LFP{imarker}.time = LFP{imarker}.time - t0;
    
    cfg_temp = [];
    cfg_temp.latency = cfg.align.latency.Hspike;
    LFP_sel = ft_selectdata(cfg_temp, LFP{imarker});
    
    subplot(2, size(LFP, 2), imarker + size(LFP, 2)); hold;
    plot(LFP_sel.time, LFP_sel.avg', 'k');
    plot(LFP_sel.time, LFP_sel.avg(chani, :)', 'r');

    dat(imarker, :) = LFP_sel.avg(chani, :);
    clear LFP_sel;
end
fname_fig = fullfile(cfg.imagesavedir, 'alignment_cluster', [cfg.prefix, 'clusters_tzeroed.png']);
isdir_or_mkdir(fileparts(fname_fig));
exportgraphics(fig, fname_fig);
close(fig)
