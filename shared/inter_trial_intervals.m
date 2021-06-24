function [intervals] = inter_trial_intervals(cfg, MuseStruct, force)

%% Make overview of data as segmented by markers aligned to (first) peak
binwidths = [1, 1, 0.25, 1, 1, 0.25, 1, 0.25];
binlimits = [60, 60, 15, 60, 60, 15, 60, 15];

fname = fullfile(cfg.datasavedir, [cfg.prefix,'timings_intervals.mat']);

if nargin == 1
    if exist(fname, 'file')
        fprintf('Reading %s\n', fname);
        load(fname, 'intervals');
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        return
    end
end

if exist(fname, 'file') && force == false
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
                    if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).events) ...
                            && MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).events > 0
                        Starttime = [Starttime; MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime'];
                        Endtime   = [Endtime;   MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime'];
                    end
                end
            end
        end
        temp.(markername){ipart}           = [Starttime, Endtime, Endtime-Starttime];
        temp.(markername){ipart}(:,4)      = [diff(Starttime); nan];
    end
end

for markername = string(cfg.muse.name)
    intervals.table.(markername) = table;
    intervals.table.(markername).Starttime  = cat(1, temp.(markername){:}(:,1));
    intervals.table.(markername).Endtime    = cat(1, temp.(markername){:}(:,2));
    intervals.table.(markername).interval   = cat(1, temp.(markername){:}(:,4));
    toclear = intervals.table.(markername).interval < 0;
    intervals.table.(markername)(toclear, :) = [];
    intervals.median.(markername) = nanmedian(intervals.table.(markername).interval);
    intervals.mean.(markername) = nanmean(intervals.table.(markername).interval);
    intervals.min.(markername)  = min(intervals.table.(markername).interval);
    intervals.max.(markername)  = max(intervals.table.(markername).interval);
    intervals.mode.(markername) = mode(intervals.table.(markername).interval);
    intervals.std.(markername)  = nanstd(intervals.table.(markername).interval);
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
        histogram(intervals.table.(markername).interval, 'BinLimits', [0 cfg.interval.histlim.(markername)], 'BinWidth', cfg.interval.histbin.(markername), 'EdgeColor', 'black', 'facecolor','k');
        axis tight
        title(sprintf('%s\nMedian: %1.2fs\nMode: %1.2fs\nMean: %1.2fs\nSD: %1.2fs', markername, intervals.median.(markername), intervals.mode.(markername), intervals.mean.(markername), intervals.std.(markername)));
        xlabel('Seconds');
        ylabel('Count');
        box off
        iplot = iplot + 1;
    end

    fname = fullfile(cfg.imagesavedir, [cfg.prefix, 'p', num2str(ipart), '_intervals'] );
    exportgraphics(fig, [fname, '.tiff'], 'Resolution', 150);
    
end

disp('done');

