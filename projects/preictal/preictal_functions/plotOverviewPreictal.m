function plotOverviewPreictal(config, MuseStruct, SpikeTrials_windowed, SpikeStats_windowed)

% threshold artefacted based on total time artefacted per window
for ipart = 1 : size(config.directorylist, 2)
    
    for iunit = 1 : size(SpikeStats_windowed{ipart}.window, 2)
        SpikeStats_windowed{ipart}.window{iunit}.trialinfo.artefact = ...
            SpikeStats_windowed{ipart}.window{iunit}.trialinfo.BAD_sec > ...
            config.plot.artthresh; % in seconds
    end
    
    % create time-axis relative to seizures
    t_seizure = MuseStruct{1}{2}.markers.CriseStart.clock;
    
    time = seconds(SpikeTrials_windowed{ipart}.window.trialinfo.starttime ...
        + (SpikeTrials_windowed{ipart}.window.trialinfo.endtime ...
        - SpikeTrials_windowed{ipart}.window.trialinfo.starttime) / 2 ...
        - t_seizure);
    
    % artefact timings
    i = 1;
    clear bad_t1 bad_t2
    for idir = 1 : size(MuseStruct{ipart}, 2)
        for ibad = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.clock, 2)
            bad_t1(i) = seconds(MuseStruct{ipart}{idir}.markers.BAD__START__.clock(ibad)    - t_seizure);
            bad_t2(i) = seconds(MuseStruct{ipart}{idir}.markers.BAD__END__.clock(ibad)      - t_seizure);
            i = i + 1;
        end
    end
    
    % colormap
    cm = lines(size(SpikeTrials_windowed{ipart}.window.label, 2));
    
    for iunit = 1 : size(SpikeStats_windowed{ipart}.window, 2)
        
        artefact = SpikeStats_windowed{ipart}.window{iunit}.trialinfo.artefact;
        
        % plot selection
        fig = figure;
        
        subplot(4,1,1);
        title(sprintf('FR unit %d "%s"', iunit, deblank(SpikeTrials_windowed{ipart}.window.cluster_group{iunit})));
        hold on;
        plot(time, SpikeStats_windowed{ipart}.window{iunit}.trialfreq, 'k');
        plot(time(artefact), SpikeStats_windowed{ipart}.window{iunit}.trialfreq(artefact), 'xr');
        
        xlim([time(1), time(end)]);
        ylim([0, 10]);
        y = ylim;
        for ibad = 1 : size(bad_t1, 2)
            fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'k', 'edgecolor', 'none', 'facealpha', 0.2);
        end
        
        subplot(4,1,2);
        title('CV2');
        hold on;
        plot([time(1), time(end)], [1 1],':k');
        plot(time, SpikeStats_windowed{ipart}.window{iunit}.CV2_trial, 'k');
        plot(time(artefact), SpikeStats_windowed{ipart}.window{iunit}.CV2_trial(artefact), 'xr');
        xlim([time(1), time(end)]);  ylim([.2, 1.8]);
        y = ylim;
        for ibad = 1 : size(bad_t1, 2)
            fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'k', 'edgecolor', 'none', 'facealpha', 0.2);
        end
        
        subplot(4,1,3);
        title('Nr. of Bursts');
        hold on;
        plot(time, SpikeStats_windowed{ipart}.window{iunit}.burst_trialsum, 'k');
        plot(time(artefact), SpikeStats_windowed{ipart}.window{iunit}.burst_trialsum(artefact), 'xr');
        xlim([time(1), time(end)]);
        y = ylim;
        for ibad = 1 : size(bad_t1, 2)
            fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'k', 'edgecolor', 'none', 'facealpha', 0.2);
        end
        
        subplot(4,1,4);
        title(sprintf('Spike distance (from %s)', SpikeStats_windowed{ipart}.window{iunit}.label),'interpreter', 'none');
        hold on;
        for ichan = 1 : size(SpikeStats_windowed{ipart}.window{iunit}.dist, 2)
            ci = strcmp(SpikeTrials_windowed{ipart}.window.label, SpikeStats_windowed{ipart}.window{ichan}.label);
            plot(time, SpikeStats_windowed{ipart}.window{iunit}.dist(:, ichan), 'color', cm(ci, :));
        end
        plot(time(artefact), SpikeStats_windowed{ipart}.window{iunit}.dist(artefact, :), 'xr');
        
        xlim([time(1), time(end)]);
        ylim([0, 1]);
        y = ylim;
        for ibad = 1 : size(bad_t1, 2)
            fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'k', 'edgecolor', 'none', 'facealpha', 0.2);
        end
        legend(SpikeStats_windowed{ipart}.window{iunit}.dist_label, 'interpreter', 'none');
        
        % save to file
        fname = fullfile(config.imagesavedir, sprintf('%soverview_part%d_unit%d', config.prefix, ipart, iunit));
        exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
        exportgraphics(fig, strcat(fname, '.pdf'));
        
    end % iunit
    
end % ipart