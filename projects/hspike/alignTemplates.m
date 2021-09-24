function [cfg, MuseStruct] = alignTemplates(cfg, MuseStruct, LFP)

% t-zero LFPs
fig = figure;
for ipart = 1 : 3
    for imarker = 1 : size(LFP{ipart}, 2)
        if isempty(LFP{ipart}{imarker})
            continue
        end
        chani       = find(strcmp(LFP{ipart}{imarker}.label, cfg.align.zerochannel));
        subplot(3, 6, imarker + (ipart-1) * 6); hold;
        plot(LFP{ipart}{imarker}.time, LFP{ipart}{imarker}.avg', 'k');
        plot(LFP{ipart}{imarker}.time, LFP{ipart}{imarker}.avg(chani, :)', 'r');
        timeindx    = LFP{ipart}{imarker}.time > -0.1 & LFP{ipart}{imarker}.time < 0.1;
        [~, LOCS] = findpeaks(-LFP{ipart}{imarker}.avg(chani, timeindx), 'SortStr', 'descend');
        if isempty(LOCS)
            continue
        end
        i0 = find(LFP{ipart}{imarker}.time > -0.1, 1, 'first') + LOCS(1) - 1;
        t0 = LFP{ipart}{imarker}.time(i0);
        axis tight; ax = axis;
        plot([t0, t0], [1.2 * ax(3), 1.2 * ax(4)], ':r');
        title(sprintf('Shift = %.0fms', -t0*1000), 'FontSize', 6);
        for idir = 1 : size(MuseStruct{ipart}, 2)
            MuseStruct{ipart}{idir}.markers.(['template', num2str(imarker)]).t0_shift = t0;
            MuseStruct{ipart}{idir}.markers.(['template', num2str(imarker)]).synctime = MuseStruct{ipart}{idir}.markers.(['template', num2str(imarker)]).synctime + t0;
            MuseStruct{ipart}{idir}.markers.(['template', num2str(imarker)]).clock    = MuseStruct{ipart}{idir}.markers.(['template', num2str(imarker)]).clock + seconds(t0);
        end
    end
end
fname = fullfile(cfg.imagesavedir, 'alignment', [cfg.prefix, 'tzeroed']);
disp(['Exporting figure ', fname])
% exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff'], 'resolution', 150);
close(fig)

for ipart = 1 : 3
    for imarker = 1 : size(LFP{ipart}, 2)
        if isempty(LFP{ipart}{imarker})
            continue
        end
        chani       = find(strcmp(LFP{ipart}{imarker}.label, cfg.align.zerochannel));
        dat(ipart, imarker, :) = LFP{ipart}{imarker}.avg(chani, :);
        
    end
end

% select templates
dat                             = permute(mean(dat, 1), [2, 3, 1]);
x                               = corr(dat', mean(dat)');
labels = [];
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


