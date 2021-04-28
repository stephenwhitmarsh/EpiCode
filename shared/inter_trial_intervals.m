function [intervals] = inter_trial_intervals(cfg, MuseStruct, force)

%% Make overview of data as segmented by markers aligned to (first) peak
binwidths = [1, 1, 0.25, 1, 1, 0.25, 1, 0.25];
binlimits = [60, 60, 15, 60, 60, 15, 60, 15];
fname = fullfile(cfg.datasavedir, [cfg.prefix,'timings_intervals.mat']);

if exist(fname, 'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed intervals ***\n');
    fprintf('************************************\n\n');
    load(fname, 'intervals');
    return
end

fprintf('********************************\n');
fprintf('*** (re-) computing intervals ***\n');
fprintf('********************************\n\n');

intervals = [];
ncols = size(cfg.muse.name, 2);
nrows = 3;

for ipart = 1 : size(MuseStruct, 2)
    
    for markername = string(cfg.muse.name)
        
        Starttime     = [];
        Endtime       = [];
        for idir = 1 : size(MuseStruct{ipart}, 2)
            if isfield(MuseStruct{ipart}{idir}, 'markers')
                if isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(markername))
                    if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).events)
                        Starttime = [Starttime; MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime'];
                        Endtime   = [Endtime;   MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime'];
                    end
                end
            end
        end
        intervals.(markername)           = [Starttime, Endtime, Endtime-Starttime];
        intervals.(markername)(:,4)      = [diff(Starttime); nan];
        intervals.(markername)(intervals.((markername))(:,4) < 0, :) = nan;
    end
end

save(fname,'intervals','-v7.3');



fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0.1 0.1 0.9 0.9]);

iplot = 1;

for ipart = 1 : size(MuseStruct, 2)
    
    %     subplot(nrows, ncols, nrows * ncols);
    tiledlayout(3, 3);
    
    for markername = string(cfg.muse.name)
        
        subplot(nrows, ncols, iplot + (ipart-1) * ncols);
        histogram(intervals.(markername)(:,4), 'BinLimits', [0 cfg.interval.histlim.(markername)], 'BinWidth', cfg.interval.histbin.(markername), 'EdgeColor', 'black', 'facecolor','k');
        %             histogram(intervals.(markername)(:,4), 'BinLimits',[0 10],'BinWidth', 0.25,'EdgeColor','black','facecolor','k');
        %             histogram(intervals.(markername)(:,4), 'EdgeColor','black','facecolor','k');
        me = nanmedian(intervals.(markername)(:,4));
        mo = mode(intervals.(markername)(:,4));
        m  = nanmean(intervals.(markername)(:,4));
        sd = nanstd(intervals.(markername)(:,4));
        axis tight
        title(sprintf('%s, Median: %1.2f, Mode: %1.2f, Mean: %1.2f, SD: %1.2f', markername, me, mo, m, sd));
        title(sprintf('%s', markername));
        xlabel('Seconds');
        ylabel('Count');
        box off
        iplot = iplot + 1;
    end
end


fname = fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_intervals'] );
%             exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff'], 'Resolution', 150);
disp('done');

