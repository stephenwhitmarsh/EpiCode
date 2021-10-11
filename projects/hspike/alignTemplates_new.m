function [cfg, LFP] = alignTemplates(cfg, LFP)

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
    chani       = find(contains(LFP{imarker}.label, cfg.align.zerochannel));
    subplot(1, size(LFP, 2), imarker); hold;
    
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
    
end
fname = fullfile(cfg.imagesavedir, 'alignment_clusters', [cfg.prefix, 'clusters_tzeroed']);
disp(['Exporting figure ', fname])
exportgraphics(fig, [fname, '.png'], 'resolution', 150);
close(fig)


% select templates; 
% FIX: THIS SHOULD BE BASED ON REALLIGNED DATA
dat     = permute(mean(dat, 1), [2, 3, 1]);
x       = corr(dat', mean(dat)');
labels  = [];

for i = 1 : size(LFP{ipart}, 2)
    labels{i} = ['template', num2str(i)];
end

for ipart = 1 : 3
    cfg.template.selected   = labels(x > 0.7);
    cfg.template.corr       = x;
end

% x = corr(dat');
%
% for itemp = 1 : size(x, 1)
%     t = x(itemp, :);
%     t(itemp) = [];
%     c(itemp) = mean(t);
%     s(itemp) = std(t);
% % end
%
% x = corr(dat', mean(dat)');
% t = mean(x) - std(x)


