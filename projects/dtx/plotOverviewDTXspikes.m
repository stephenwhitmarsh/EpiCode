function plotOverviewDTXspikes(cfg,SpikeTrials_timelocked, SpikeTrials_windowed, SpikeStats_windowed, SpikeDensity_timelocked, LFP, TFR, SpikeWaveforms, SpikeWaveforms_stats, MuseStruct_concat, seizure_infos)

if isempty(LFP)
    istimelocked = false;
    event_name  = "";
    baseline_name   = "Control";
    nplots_1 = 8;
    nplots_2 = 7;
else
    istimelocked = true;
    event_name      = "SlowWave";
    baseline_name   = "Interictal";
    nplots_1 = 9;
    nplots_2 = 7;
end

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])
set(0,'defaultfigurecolor',[1 1 1]);

%% go through each part and each units
for ipart = 1:size(SpikeTrials_windowed,2)
   
    if istimelocked
        %Prepare TFR
        %nanmean because fieldtrip does only mean
        if strcmp(TFR{ipart}.(baseline_name).dimord, 'rpt_chan_freq_time')
            TFR{ipart}.(baseline_name).powspctrm = nanmean(TFR{ipart}.(baseline_name).powspctrm,1);
            TFR{ipart}.(baseline_name).powspctrm = permute(TFR{ipart}.(baseline_name).powspctrm, [2 3 4 1]);
            TFR{ipart}.(baseline_name).dimord    = erase(TFR{ipart}.(baseline_name).dimord, 'rpt_');
        elseif ~strcmp(TFR{ipart}.(baseline_name).dimord, 'chan_freq_time')
            error('dimord is expected to be chan_freq_time for tfr data')
        end
    end
        
    for itemp = 1:size(SpikeTrials_windowed{ipart}.window.label, 2)
        
        fprintf('Plot overview for template %d from %d\n', itemp, size(SpikeTrials_windowed{ipart}.window.label, 2));
        fig = figure;
        
        if istimelocked
            %% Nb of events per hour
            subplot(nplots_1, nplots_2, 1:6); hold on;
            title(sprintf('"%s" rate', char(event_name)));
            %endtime so the plot is more precise about when the first seizures occur
            bar(seizure_infos.statsovertime.endtime, seizure_infos.statsovertime.nb_seizures, 1, 'EdgeColor', 'none', 'FaceColor', [0, 0, 0]);
            ylabel('event / hr.');
            axis tight
            yticklabels(yticks .* 2);
            x = xlim;
            xlim([MuseStruct_concat{ipart}.markers.Analysis_Start.clock(end), x(2)])
            y = ylim;
            x = cfg.injectiontime;
            plot([x x], y, '-r');
            xl = xlim;
            
            %% distrib freq of events
            subplot(nplots_1, nplots_2, 7); hold on;
            title(sprintf('"%s" per hr.', char(event_name)));
            
            m = -inf;
            t           = seizure_infos.statsovertime.endtime > cfg.injectiontime;
            y           = nanmean(seizure_infos.statsovertime.nb_seizures(t));
            stdev       = nanstd(seizure_infos.statsovertime.nb_seizures(t));
            hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
            m           = max([m, y, y+stdev]);
            he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
            xlim([-0.7, 2.7]); ylim([0, m * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        else
            xl = [MuseStruct_concat{ipart}.markers.Analysis_Start.clock(end), MuseStruct_concat{ipart}.markers.Analysis_End.clock(end)];
        end
        
        %% Firing rate along all data
        subplot(nplots_1, nplots_2, 8:13); hold on;
        title('Firingrate');
        idx = SpikeStats_windowed{ipart}.window{itemp}.trialfreq < 100;
        plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime(idx), SpikeStats_windowed{ipart}.window{itemp}.trialfreq(idx), 'k');
        ylabel('Hz');
        disp(xl);
        xlim(xl);
        set(gca,'TickLength',[0 0], 'Xticklabels', []);
        y = ylim;
        x = cfg.injectiontime;
        plot([x x], y, '-r');

        %% Average firing rate
        subplot(nplots_1, nplots_2, 14); hold on;
        title('Firing rate (Hz)');
        
        m = -inf;
        y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.trialfreq(idx));
        stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.trialfreq(idx));
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, m * 1.1]);
        set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% Amplitude over time
        subplot(nplots_1, nplots_2, 15:20); hold on;
        title('Amplitude');
        plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.amplitude, 'Color', [0, 0, 0]);
        ylabel('uV');
        axis tight
        xlim(xl);
        set(gca,'TickLength',[0 0], 'Xticklabels', []);
        y = ylim;
        x = cfg.injectiontime;
        plot([x x], y, '-r');
   
        %% Averaged amplitude
        subplot(nplots_1, nplots_2, 21); hold on;
        title('Amplitude (uV)');
        
        ymax = -inf;
        ymin = +inf;
        y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.amplitude);
        stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.amplitude);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        ymax        = max([ymax, y+stdev]);
        ymin        = min([ymin, y-stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([ymin * 0.9, ymax * 1.1]);
        set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');

        %% CV2 over time
        subplot(nplots_1, nplots_2, 22:27); hold on;
        title('CV2');
        plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.CV2_trial', 'Color', [0, 0, 0]);
        plot(xl,[1 1],'k--');
        ylabel('CV2');
        axis tight
        xlim(xl);
        set(gca,'TickLength',[0 0], 'Xticklabels', []);
        ylim([0 2]);
        y = ylim;
        x = cfg.injectiontime;
        plot([x x], y, '-r');
        
        %% Averaged CV2
        subplot(nplots_1, nplots_2, 28); hold on;
        title('CV2');
        
        y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.CV2_trial);
        stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.CV2_trial);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, 2]);
        set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');

        %% Nr of bursts per minute
        subplot(nplots_1, nplots_2, 29:34); hold on;
        title('Number of bursts per minute');
        plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum', 'Color', [0, 0, 0]);
        ylabel('burst / min.');
        axis tight
        xlim(xl);
        set(gca,'TickLength',[0 0], 'Xticklabels', []);
        y = ylim;
        x = cfg.injectiontime;
        plot([x x], y, '-r');
        
        %% Averaged Bursts per minute
        subplot(nplots_1, nplots_2, 35); hold on;
        title('Bursts per minute');
        
        m = -inf;
        y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum);
        stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); 
        try ylim([0, m * 1.1]); end
        set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
       
        %% spike waveform : plot average of all waveforms, and plot 1000 randomly selected trials
        if istimelocked
            subplot(nplots_1, nplots_2,  [42 49]); hold on;
        else
            subplot(nplots_1, nplots_2,  [45 52]); hold on;
        end
        title('Action Potential Waveform');
        if size(SpikeWaveforms{ipart}.window{itemp}.trial, 2) > 1000
            spikes_idx_sel = randperm(size(SpikeWaveforms{ipart}.window{itemp}.trial,2),1000);
        else
            spikes_idx_sel = 1:length(SpikeWaveforms{ipart}.window{itemp}.trial);
        end
        for itrial = spikes_idx_sel
            lh = plot(SpikeWaveforms{ipart}.window{itemp}.time{itrial}, SpikeWaveforms{ipart}.window{itemp}.trial{itrial}, 'color', [0.3 0.3 0.3]);
            lh.Color = [0,0,0,0.1];
        end
        if isfield(SpikeWaveforms_stats{ipart}.window.waveformavg{itemp}, 'time')
            plot(SpikeWaveforms_stats{ipart}.window.waveformavg{itemp}.time, SpikeWaveforms_stats{ipart}.window.waveformavg{itemp}.avg, 'color', [1 1 1 ], 'linewidth', 1);
            %scatter(SpikeWaveforms_stats{ipart}.window.halfwidth.x(itemp,:), SpikeWaveforms_stats{ipart}.window.halfwidth.y(itemp,:), 'x', 'MarkerEdgeColor', 'y');
            %scatter(SpikeWaveforms_stats{ipart}.window.peaktrough.x(itemp,:), SpikeWaveforms_stats{ipart}.window.peaktrough.y(itemp,:), 'x', 'MarkerEdgeColor', 'y');
            %scatter(SpikeWaveforms_stats{ipart}.window.troughpeak.x(itemp,:), SpikeWaveforms_stats{ipart}.window.troughpeak.y(itemp,:), 'x', 'MarkerEdgeColor', 'y');
            axis tight
            ylim([min(SpikeWaveforms_stats{ipart}.window.waveformavg{itemp}.avg)*2 max(SpikeWaveforms_stats{ipart}.window.waveformavg{itemp}.avg)]*2);
        end
        xticklabels(xticks.*1000); 
        %xlabel('Time (ms)'); 
        ylabel('uV');
        set(gca, 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
        
        %plot template used by spykng circus
        %         y = SpikeTrials_windowed{ipart}.window.template{itemp}(:,SpikeTrials_windowed{ipart}.window.template_maxchan(itemp)+1,:);%+1 because electrodes nr are zero-based
        %         x = ( (0 : size(SpikeTrials_windowed{ipart}.window.template{itemp},3) - 1) / SpikeTrials_windowed{ipart}.window.hdr.Fs )';
        %         x = x - x(end)/2;
        %         plot (squeeze(x)',squeeze(y)','r','LineWidth', 2);
        %         ylim([min(y)*2, max(y)]*2);

        %% ISI histogram
        if istimelocked
            subplot(nplots_1, nplots_2,  [56 63]); hold on;
        else
            subplot(nplots_1, nplots_2,  [47 54]); hold on;
        end
        title(sprintf('ISI (RPV = %.2f%%)',SpikeStats_windowed{ipart}.window{itemp}.RPV * 100));
        bar(SpikeStats_windowed{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats_windowed{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0 0 0], 'edgecolor', 'none');
        idx = SpikeStats_windowed{ipart}.window{itemp}.isi_avg_time <= cfg.spike.RPV;
        bar(SpikeStats_windowed{ipart}.window{itemp}.isi_avg_time(idx) * 1000, SpikeStats_windowed{ipart}.window{itemp}.isi_avg(idx), 1, 'facecolor', [1 0 0], 'edgecolor', 'none');
        set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
        axis tight
        xlim([0, 0.025] * 1000);
        xlabel('Time (ms)'); ylabel('Spikecount');
        
        %% Titles for events and baseline
        if istimelocked
            subplot(nplots_1, nplots_2, 36:38);
            axis off
            str = sprintf('%s : %d spikes', event_name, size(SpikeTrials_timelocked{ipart}.(event_name).time{itemp},2));
            text(0.5, 0.4, str, 'HorizontalAlignment', 'center', 'FontSize', sgtitle_obj.FontSize, 'FontWeight', sgtitle_obj.FontWeight);
            
            subplot(nplots_1, nplots_2, 39:41);
            axis off
            str = sprintf('%s : %d spikes', baseline_name, size(SpikeTrials_timelocked{ipart}.(baseline_name).time{itemp},2));
            text(0.5, 0.4, str, 'HorizontalAlignment', 'center', 'FontSize', sgtitle_obj.FontSize, 'FontWeight', sgtitle_obj.FontWeight);
        else
            subplot(nplots_1, nplots_2, 36:42);
            axis off
            str = sprintf('%s : %d spikes', baseline_name, size(SpikeTrials_windowed{ipart}.window.time{itemp},2));
            text(0.5, 0.4, str, 'HorizontalAlignment', 'center', 'FontSize', sgtitle_obj.FontSize, 'FontWeight', sgtitle_obj.FontWeight);
        end
          
        if istimelocked
            %% Event LFP & TFR
            h = subplot(nplots_1, nplots_2,  [43, 50, 57]); hold on;
            
            cfgtemp             = [];
            cfgtemp.latency     = cfg.epoch.toi.(event_name);
            LFP_avg             = ft_timelockanalysis(cfgtemp, LFP{ipart}.(event_name));
            
            chansel = find(contains(TFR{ipart}.(event_name).label,channel));
            cfgtemp                 = [];
            cfgtemp.channel         = chansel;
            cfgtemp.baseline        = cfg.LFP.baselinewindow.(event_name);
            cfgtemp.baselinetype    = 'relchange';
            cfgtemp.ylim            = [10, 100];
            cfgtemp.zlim            = 'maxabs';
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.figure          = h;
            cfgtemp.interactive     = 'no';
            ft_singleplotTFR(cfgtemp, TFR{ipart}.(event_name));
            xlabel('Time (s)'); ylabel('Frequency');
            c = colorbar;
            set(c, 'Location', 'southoutside', 'color', [0 0 0]);
            c.Title.String = 'Relative change in power';
            pos = get(c, 'Position');
            pos(2) = 0.03;
            set(c, 'pos', pos);
            colormap(jet(1000));
            title('TFR & LFP');
            
            % plot single average LFP
            yyaxis right
            set(gca,'YColor','k');
            hold on
            n = 1;
            ytick = [];
            label = [];
            l = [];
            chansel = contains(LFP_avg.label,channel);
            x = LFP_avg.time;
            y = LFP_avg.avg(chansel, :);
            ystd = sqrt(LFP_avg.var(chansel, :));
            patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], [0 0 0], 'edgecolor', 'none', 'facealpha', 0.2);
            plot(x, y, 'color', [0 0 0],'linewidth', 0.5);
            h = ylabel('uV');
            h.Rotation = -90;
            ymax = max(abs([y + ystd, y - ystd]));
            ylim([-ymax, ymax]);
            xlim([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)]);
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
            
            %% Event raster
            subplot(nplots_1, nplots_2,  [44, 45, 51 52]); hold on;
            cfgtemp                 = [];
            cfgtemp.spikechannel    = itemp;
            cfgtemp.latency         = [cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)];
            cfgtemp.trialborders    = 'yes';
            cfgtemp.linewidth       = 1;
            ft_spike_plot_raster(cfgtemp, SpikeTrials_timelocked{ipart}.(event_name));
            set(gca, 'Yticklabels', [], 'Ytick', [], 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
            title('Rasterplot');
            ylabel('');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
            xlabel([]);
            
            %% Event spike density and stats
            subplot(nplots_1, nplots_2,  [58 59]); hold on;
            title('Spike Density (Hz)');
            bar(SpikeDensity_timelocked{ipart}.psth.(event_name).time, squeeze(nanmean(SpikeDensity_timelocked{ipart}.psth.(event_name).trial(:, itemp, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
            %ylabel('Firing rate (Hz)');
            xlabel('Time (s)');
            
            % plot positive clusters
            lag = size(SpikeDensity_timelocked{ipart}.psth.(event_name).time, 2) - size(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.time, 2);
            
            if isfield(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}, 'posclusters')
                for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.posclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.posclusterslabelmat == ipos);
                        bar(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.time(sel), SpikeDensity_timelocked{ipart}.psth.(event_name).avg(itemp, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
                    end
                end
            end
            
            % plot negative clusters
            if isfield(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}, 'negclusters')
                for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.negclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.negclusterslabelmat == ineg);
                        bar(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.time(sel), SpikeDensity_timelocked{ipart}.psth.(event_name).avg(itemp, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
                    end
                end
            end
            
            % plot smoothed
            plot(SpikeDensity_timelocked{ipart}.sdf_bar.(event_name).time, squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_bar.(event_name).trial(:, itemp, :), 1))', 'Color', [0 0 0]);
            
            % plot baseline patch
            x = [cfg.stats.bl.(event_name)(1) cfg.stats.bl.(event_name)(2) cfg.stats.bl.(event_name)(2) cfg.stats.bl.(event_name)(1)];
            ax = axis;
            y = [ax(3) ax(3) ax(4) ax(4)];
            patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
            set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
            
            % plot baseline indicator line
            y = mean(SpikeDensity_timelocked{ipart}.stat.(event_name){itemp}.baseline(:, itemp));
            plot([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)], [y, y], ':k');
            clear y
            
            % add LFP
            yyaxis right
            x = LFP_avg.time;
            y = LFP_avg.avg(chansel, :);
            plot(x,y,'-g');
            set(gca,'YColor','k', 'TickDir', 'out');
            ylabel('uV');
            xlim([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)]);
            
            %% Baseline TFR
            h = subplot(nplots_1, nplots_2,  [46, 53, 60]); hold on;
            
            chansel = find(contains(TFR{ipart}.(baseline_name).label,channel));
            cfgtemp                 = [];
            cfgtemp.channel         = chansel;
            cfgtemp.baseline        = [-100 -80];%cfg.LFP.baselinewindow.(baseline_name);
            cfgtemp.baselinetype    = 'relchange';
            cfgtemp.ylim            = [10, 100];
            cfgtemp.zlim            = 'maxabs';
            if TFR{ipart}.(baseline_name).time(1) < -100
                cfgtemp.xlim            = [-100, TFR{ipart}.(baseline_name).time(end)];
            else
                cfgtemp.xlim            = [TFR{ipart}.(baseline_name).time(1), TFR{ipart}.(baseline_name).time(end)];
            end
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.figure          = h;
            cfgtemp.interactive     = 'no';
            ft_singleplotTFR(cfgtemp, TFR{ipart}.(baseline_name));
            xlabel('Time (s)'); ylabel('Frequency');
            c = colorbar;
            set(c, 'Location', 'southoutside', 'color', [0 0 0]);
            caxis([-1 1]);
            c.Title.String = 'Power relative change';
            pos = get(c, 'Position');
            pos(2) = 0.03;
            set(c, 'pos', pos);
            colormap(jet(1000));
            title('TFR & LFP');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
            
            %% baseline : all densities lines or raster
            subplot(nplots_1, nplots_2,  [47, 48, 54, 55]); hold on;
            
%             cfgtemp                 = [];
%             cfgtemp.spikechannel    = itemp;
%             cfgtemp.latency         = 'all';
%             cfgtemp.trialborders    = 'yes';
%             cfgtemp.linewidth       = 1;
%             ft_spike_plot_raster(cfgtemp, SpikeTrials_timelocked{ipart}.(baseline_name));
%             set(gca, 'Yticklabels', [], 'Ytick', [], 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
%             title('Rasterplot');
%             ylabel('');
%             set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
%             xlabel([]);
            
            for itrial = 1:size(SpikeDensity_timelocked{ipart}.sdf_bar.(baseline_name).trial, 1)
                p = plot(SpikeDensity_timelocked{ipart}.sdf_bar.(baseline_name).time, squeeze(SpikeDensity_timelocked{ipart}.sdf_bar.(baseline_name).trial(itrial,itemp,:)), 'k');
                p.Color(4) = 0.2;
            end
            title('Spike Density per trial');
            axis tight
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
            ylabel('Firing rate (Hz)');
            
            %% baseline densities
            subplot(nplots_1, nplots_2,  [61 62]); hold on;
            title('Spike Density averaged');
            bar(SpikeDensity_timelocked{ipart}.psth.(baseline_name).time, squeeze(nanmean(SpikeDensity_timelocked{ipart}.psth.(baseline_name).trial(:, itemp, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
            ylabel('Firing rate (Hz)');
            
%             % plot positive clusters
%             lag = size(SpikeDensity_timelocked{ipart}.psth.(baseline_name).time, 2) - size(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.time, 2);
%             
%             if isfield(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}, 'posclusters')
%                 for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.posclusters, 2)
%                     if SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.posclusters(ipos).prob < cfg.stats.alpha
%                         sel = find(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.posclusterslabelmat == ipos);
%                         bar(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.time(sel), SpikeDensity_timelocked{ipart}.psth.(baseline_name).avg(itemp, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
%                     end
%                 end
%             end
%             
%             % plot negative clusters
%             if isfield(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}, 'negclusters')
%                 for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.negclusters, 2)
%                     if SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.negclusters(ineg).prob < cfg.stats.alpha
%                         sel = find(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.negclusterslabelmat == ineg);
%                         bar(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.time(sel), SpikeDensity_timelocked{ipart}.psth.(baseline_name).avg(itemp, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
%                     end
%                 end
%             end
            
            % plot smoothed
            plot(SpikeDensity_timelocked{ipart}.sdf_bar.(baseline_name).time, squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_bar.(baseline_name).trial(:, itemp, :), 1))', 'Color', [0 0 0]);
            
%             % plot baseline patch
%             x = [cfg.stats.bl.(baseline_name)(1) cfg.stats.bl.(baseline_name)(2) cfg.stats.bl.(baseline_name)(2) cfg.stats.bl.(baseline_name)(1)];
%             ax = axis;
%             y = [ax(3) ax(3) ax(4) ax(4)];
%             patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
%             set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
            
%             % plot baseline indicator line
%             y = nanmean(SpikeDensity_timelocked{ipart}.stat.(baseline_name){itemp}.baseline(:, itemp));
%             x = xlim;
%             plot(x, [y, y], ':k');
%             clear y
            
        end
        %% save figure
        fname = fullfile(cfg.imagesavedir,'spike_overview',sprintf('%sp%d-%s_overview_%s_%s', cfg.prefix, ipart,SpikeTrials_windowed{ipart}.window.label{itemp}, char(event_name), char(baseline_name)));
        if ~isfolder(fileparts(fname))
            fprintf('Creating directory %s\n', fileparts(fname));
            mkdir(fileparts(fname));
        end
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        fig.PaperType = 'A2';
        
        print(fig, '-dpdf', [fname, '.pdf'],'-r600');
        print(fig, '-dpng', [fname, '.png'],'-r600');
        
    end %itemp
end %ipart