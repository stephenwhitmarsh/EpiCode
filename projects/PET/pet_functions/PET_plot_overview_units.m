function PET_plot_overview_units(cfg,SpikeTrials, SpikeStats, SpikeDensity, LFP, TFR, FFT, SpikeWaveforms, SpikeWaveforms_stats, MuseStruct_concat)

% artefactlist if optional, and is a structure of 2 datetime arrays of the
% same size :
% artefacts_list.starttime
% artefacts_list.endtime

event_name     = 'IED';
baseline_name  = 'interIED';
do_plot_artefacts = false;
do_plot_events    = false;

nplots_1    = 12;
nplots_2    = 7;
title_size  = 7;
legend_size = 7;

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])
set(0,'defaultfigurecolor',[1 1 1]);


%% go through each part and each units
for ipart = 1:size(SpikeTrials,2)
    
    %% get artefact times if required
    i = 1;
    bad_t1 = datetime.empty;
    bad_t2 = datetime.empty;
    if ~isfield(MuseStruct_concat{ipart}.markers, 'BAD__START__')
        %do nothing
    elseif ~isfield(MuseStruct_concat{ipart}.markers.BAD__START__, 'clock')
        % do nothing
    else
        for ibad = 1 : size(MuseStruct_concat{ipart}.markers.BAD__START__.clock, 2)
            bad_t1(i) = MuseStruct_concat{ipart}.markers.BAD__START__.clock(ibad);
            bad_t2(i) = MuseStruct_concat{ipart}.markers.BAD__END__.clock(ibad);
            i = i + 1;
        end
    end

    
    %% get events time if required
    if do_plot_events
        i = 1;
        if ~isfield(MuseStruct_concat{ipart}.markers, cfg.muse.startmarker.(event_name))
            event_t1 = [];
            event_t2 = [];
        elseif ~isfield(MuseStruct_concat{ipart}.markers.(cfg.muse.startmarker.(event_name)), 'clock')
            event_t1 = [];
            event_t2 = [];
        else
            for ievent = 1 : size(MuseStruct_concat{ipart}.markers.(cfg.muse.startmarker.(event_name)).clock, 2)
                event_t1(i) = MuseStruct_concat{ipart}.markers.(cfg.muse.startmarker.(event_name)).clock(ievent) + ...
                    seconds(cfg.epoch.toi.(event_name)(1));
                event_t2(i) = MuseStruct_concat{ipart}.markers.(cfg.muse.endmarker.(event_name)).clock(ievent) + ...
                    seconds(cfg.epoch.toi.(event_name)(1));
                i = i + 1;
            end
        end
    else
        event_t1 = [];
        event_t2 = [];
    end
    
    %% go trough each unit
    
    for itemp = 1:size(SpikeTrials{ipart}.(baseline_name).label, 2)
        
        fprintf('***** Plot overview for template %d from %d *****\n', itemp, size(SpikeTrials{ipart}.(baseline_name).label, 2));
        fig = figure;
        unit_name     = SpikeTrials{ipart}.(baseline_name).label{itemp};
        maxchan       = SpikeTrials{ipart}.(baseline_name).template_maxchan(itemp) + 1; %+1 because starts at zero
        channel       = cfg.circus.channel{maxchan};
        cluster_group = strrep(SpikeTrials{ipart}.(baseline_name).cluster_group{itemp}, ' ', '');
        sgtitle_obj   = sgtitle(sprintf('%s : %s (%s)', cfg.prefix(1:end-1), unit_name, cluster_group), 'Interpreter', 'none', 'FontWeight', 'bold', 'Fontsize', 18);
        
        xl = [MuseStruct_concat{ipart}.starttime MuseStruct_concat{ipart}.endtime];
        
        %% title for baseline
        subplot(nplots_1, nplots_2, 1:2);
        axis off
        if do_plot_artefacts
            str = sprintf('%s : %d APs in %d trials (>%ds between IED, and no artefacts; red = artefacted periods)', baseline_name, size(SpikeTrials{ipart}.(baseline_name).time{itemp},2), ...
                size(SpikeTrials{ipart}.(baseline_name).trialinfo,1), cfg.min_trial_length.(baseline_name));
        else
            str = sprintf('%s : %d APs in %d trials (>%ds between IED, and no artefacts)', baseline_name, size(SpikeTrials{ipart}.(baseline_name).time{itemp},2), ...
                size(SpikeTrials{ipart}.(baseline_name).trialinfo,1), cfg.min_trial_length.(baseline_name));
        end
        text(0, 0.2, str, 'HorizontalAlignment', 'left', 'FontSize', 12, 'fontweight', 'bold');%, 'FontWeight', sgtitle_obj.FontWeight);
        
        
        %% FFT
        subplot(nplots_1, nplots_2, 8:13); hold on;
        % power 1-15 Hz
        
        % select and average data
        cfgtemp             = [];
        cfgtemp.avgoverfreq = 'yes';
        cfgtemp.channel     = ['*', channel, '*'];
        data                = ft_selectdata(cfgtemp, FFT{ipart}.(baseline_name));
        assert(length(data.label) == 1, 'error in channel name'); %with the defined parameters we should only have one channel.
        time                = FFT{ipart}.(baseline_name).trialinfo.starttime + (FFT{ipart}.(baseline_name).trialinfo.endtime - FFT{ipart}.(baseline_name).trialinfo.starttime) / 2;
        bar(time, data.powspctrm, 1, 'facecolor', 'k', 'edgecolor', 'k');
        t = ylabel('uV²');
        xlim(xl);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'box', 'off');
        plot_events(true, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        if do_plot_artefacts
            title(sprintf('LFP power (%d-%d Hz)', cfg.FFT.foi.(baseline_name)(1), cfg.FFT.foi.(baseline_name)(end)), 'verticalalignment', 'top', 'fontsize', title_size);
        else
            title(sprintf('LFP power (%d-%d Hz) and artefacted periods (in red)', cfg.FFT.foi.(baseline_name)(1), cfg.FFT.foi.(baseline_name)(end)), 'verticalalignment', 'top', 'fontsize', title_size);
        end
        
        %% mean power
        subplot(nplots_1, nplots_2, 14); hold on;
        title(sprintf('Mean power (%d-%d Hz)', cfg.FFT.foi.(baseline_name)(1), cfg.FFT.foi.(baseline_name)(end)), 'fontsize', title_size);
        
        m = -inf;
        y           = mean(data.powspctrm);
        stdev       = std(data.powspctrm);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, m * 1.1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none',...
            'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        
        %% Firing rate along all data
        subplot(nplots_1, nplots_2, 15:20); hold on;
        title('Firing rate', 'verticalalignment', 'top', 'fontsize', title_size);
        idx = SpikeStats{ipart}.(baseline_name){itemp}.trialfreq < 100;
        plot(SpikeStats{ipart}.(baseline_name){itemp}.trialinfo.starttime(idx), SpikeStats{ipart}.(baseline_name){itemp}.trialfreq(idx), '-k');
        ylabel('Hz');
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', []);
        plot_events(do_plot_artefacts, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        xlim(xl);
        
        %% Average firing rate
        subplot(nplots_1, nplots_2, 21); hold on;
        title('Firing rate (Hz)', 'fontsize', title_size);
        
        m = -inf;
        y           = nanmean(SpikeStats{ipart}.(baseline_name){itemp}.trialfreq(idx));
        stdev       = nanstd(SpikeStats{ipart}.(baseline_name){itemp}.trialfreq(idx));
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, m * 1.1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% Amplitude over time
        subplot(nplots_1, nplots_2, 22:27); hold on;
        title('Amplitude', 'verticalalignment', 'top', 'fontsize', title_size);
        plot(SpikeStats{ipart}.(baseline_name){itemp}.trialinfo.starttime, SpikeStats{ipart}.(baseline_name){itemp}.amplitude, '-k');
        ylabel('uV');
        axis tight
        xlim(xl);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', []);
        y = ylim;
        plot_events(do_plot_artefacts, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        
        %% Averaged amplitude
        subplot(nplots_1, nplots_2, 28); hold on;
        title('Amplitude (uV)', 'fontsize', title_size);
        
        ymax = -inf;
        ymin = +inf;
        y           = nanmean(SpikeStats{ipart}.(baseline_name){itemp}.amplitude);
        stdev       = nanstd(SpikeStats{ipart}.(baseline_name){itemp}.amplitude);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        ymax        = max([ymax, y+stdev]);
        ymin        = min([ymin, y-stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([ymin * 0.9, ymax * 1.1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% CV2 over time
        subplot(nplots_1, nplots_2, 29:34); hold on;
        title('CV2', 'verticalalignment', 'top', 'fontsize', title_size);
        plot(SpikeStats{ipart}.(baseline_name){itemp}.trialinfo.starttime, SpikeStats{ipart}.(baseline_name){itemp}.CV2_trial', '-k');
        plot(xl,[1 1],'k--');
        ylabel('CV2');
        axis tight
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', []);
        xlim(xl);
        ylim([0 2]);
        plot_events(do_plot_artefacts, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        
        %% Averaged CV2
        subplot(nplots_1, nplots_2, 35); hold on;
        title('CV2', 'fontsize', title_size);
        
        y           = nanmean(SpikeStats{ipart}.(baseline_name){itemp}.CV2_trial);
        stdev       = nanstd(SpikeStats{ipart}.(baseline_name){itemp}.CV2_trial);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, 2]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% Nr of bursts per minute
        subplot(nplots_1, nplots_2, 36:41); hold on;
        title('Number of bursts per minute', 'verticalalignment', 'top', 'fontsize', title_size);
        plot(SpikeStats{ipart}.(baseline_name){itemp}.trialinfo.starttime, SpikeStats{ipart}.(baseline_name){itemp}.burst_trialsum', '-k');
        ylabel('burst / min.');
        axis tight
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', []);
        xlim(xl);
        plot_events(do_plot_artefacts, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        
        %% Averaged Bursts per minute
        subplot(nplots_1, nplots_2, 42); hold on;
        title('Bursts per minute', 'fontsize', title_size);
        
        m = -inf;
        y           = nanmean(SpikeStats{ipart}.(baseline_name){itemp}.burst_trialsum);
        stdev       = nanstd(SpikeStats{ipart}.(baseline_name){itemp}.burst_trialsum);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]);
        try ylim([0, m * 1.1]); end
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% spike distance over time
        subplot(nplots_1, nplots_2, 43:48); hold on;
        title('Spike distance to each other unit', 'verticalalignment', 'top', 'fontsize', title_size);
        
        cm  = lines(size(SpikeStats{ipart}.interIED, 2));
        for i = 1 : size(SpikeStats{ipart}.interIED, 2)
            l{i} = SpikeStats{ipart}.interIED{i}.label;
        end
        
        x = SpikeStats{ipart}.(baseline_name){itemp}.trialinfo.starttime;
        % plot data
        y = SpikeStats{ipart}.(baseline_name){itemp}.dist;
        
        leg = [];
        for itarget = 1:size(SpikeStats{ipart}.interIED{itemp}.dist_label, 2)
            ci = strcmp(SpikeStats{ipart}.interIED{itemp}.dist_label{itarget}, SpikeTrials{ipart}.interIED.label);
            plot(x, y(itarget,:), 'color', cm(ci, :));
            leg = [leg; string(SpikeStats{ipart}.interIED{itemp}.dist_label{itarget})];
        end
        
        xlim(xl);
        %         y = ylim;
        %         ylim([0 y(2)]);
        ylim([0 1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0]);
        xlabel('Time');
        plot_events(do_plot_artefacts, bad_t1, bad_t2, 'r');
        plot_events(do_plot_events, event_t1, event_t2, 'b');
        
        %% spike distance : mean per unit
        subplot(nplots_1, nplots_2, 49); hold on;
        title('Spike distance', 'fontsize', title_size);
        
        y = SpikeStats{ipart}.(baseline_name){itemp}.dist;
        
        leg = [];
        for itarget = 1:size(SpikeStats{ipart}.interIED{itemp}.dist_label, 2)
            ci = strcmp(SpikeStats{ipart}.interIED{itemp}.dist_label{itarget}, SpikeTrials{ipart}.interIED.label);
            errorbar(2/size(SpikeStats{ipart}.interIED{itemp}.dist_label, 2)*itarget, mean(y(itarget,:)), std(y(itarget,:)), '.', 'color', cm(ci, :));
            leg = [leg; string(SpikeStats{ipart}.interIED{itemp}.dist_label{itarget})];
        end
        
        xlim([-0.3, 2.3]);
        %         y = ylim;
        %         ylim([0 y(2)]);
        ylim([0 1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        %% Title for events
        subplot(nplots_1, nplots_2, 50:51);
        axis off
        str = sprintf('%s : %d APs in %d trials (without artefacts)', event_name, size(SpikeTrials{ipart}.(event_name).time{itemp},2), ...
            size(SpikeTrials{ipart}.(event_name).trialinfo,1));
        text(0, 0.2, str, 'HorizontalAlignment', 'left', 'FontSize', 12, 'fontweight', 'bold');%, 'FontWeight', sgtitle_obj.FontWeight);
        
        %% Nb of events per minute
        s = subplot(nplots_1, nplots_2, 57:60); hold on;
        %[left bottom width height]
        %         s.Position(2) = s.Position(2)+0.05;
        title(sprintf('"%s" rate', char(event_name)), 'verticalalignment', 'top', 'fontsize', title_size);
        bins = 0 : 60 : seconds(MuseStruct_concat{ipart}.endtime - MuseStruct_concat{ipart}.starttime);
        [n, e] = histcounts(MuseStruct_concat{ipart}.markers.(cfg.muse.startmarker.(event_name)).synctime, bins);
        x = seconds(e(2:end)) + MuseStruct_concat{ipart}.starttime;
        bar(x, n, 1, 'EdgeColor', 'none', 'FaceColor', [0, 0, 0]);
        ylabel('event / min');
        axis tight
        set(gca, 'fontsize', legend_size);
        
        %% distrib freq of events
        s = subplot(nplots_1, nplots_2, 61); hold on;
        %[left bottom width height]
        %         s.Position(2) = s.Position(2)+0.05;
        t = title(sprintf('"%s" per min', char(event_name)), 'fontsize', title_size);
        
        m = -inf;
        y           = mean(n);
        stdev       = std(n);
        hb          = bar(1, y, 1); set(hb, 'FaceColor', [0.6 0.6 0.6]);
        m           = max([m, y, y+stdev]);
        he          = errorbar(1, y, stdev, 'clipping','off','color', 'k');
        xlim([-0.7, 2.7]); ylim([0, m * 1.1]);
        set(gca, 'fontsize', legend_size,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
        
        
        %% Event LFP & TFR
        h = subplot(nplots_1, nplots_2,  [64,65,71,72,78,79]); hold on;
        
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
        c.Title.FontSize = legend_size;
        
        y = ylim;
        x = xlim;
        t_x = x(1) + (x(2) - x(1))/2;
        %         t = title('TFR & LFP', 'Position', [t_x, y(2), 1],  'verticalalignment', 'top');
        title([]);
        text(t_x, y(2), 1, 'TFR & LFP', 'color', 'k',  'verticalalignment', 'top', 'horizontalalignment', 'center',...
            'fontweight', 'bold', 'fontsize', title_size, 'FontName', t.FontName);
        % plot single average LFP
        yyaxis right
        set(gca,'YColor','k');
        hold on
        n = 1;
        ytick = [];
        label = [];
        l = [];
        chansel = contains(LFP_avg.label,channel);
        assert(sum(chansel) == 1); %with the defined parameters we should only have one channel.
        x = LFP_avg.time;
        y = LFP_avg.avg(chansel, :);
        ystd = sqrt(LFP_avg.var(chansel, :));
        patch_std(x, y, ystd, [0 0 0])
        %         patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], [0 0 0], 'edgecolor', 'none', 'facealpha', 0.2);
        plot(x, y, '-', 'color', 'k','linewidth', 1);
        h = ylabel('uV');
        h.Rotation = -90;
        ymax = max(abs([y + ystd, y - ystd]));
        ylim([-ymax, ymax]);
        xlim([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)]);
        set(gca, 'fontsize', legend_size, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
        
        %         t = title('TFR & LFP', 'verticalalignment', 'bottom');
        % yyaxis left
        %         t = title('TFR & LFP');%, 'verticalalignment', 'bottom');
        %         t.Position(3) = 100;
        
        %% Event raster
        subplot(nplots_1, nplots_2,  [66,67,68,73,74,75]); hold on;
        cfgtemp                 = [];
        cfgtemp.spikechannel    = itemp;
        cfgtemp.latency         = [cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)];
        cfgtemp.trialborders    = 'yes';
        cfgtemp.linewidth       = 1;
        ft_spike_plot_raster(cfgtemp, SpikeTrials{ipart}.(event_name));
        set(gca, 'Yticklabels', [], 'Ytick', [], 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
        title('Rasterplot', 'verticalalignment', 'top', 'fontsize', title_size);
        ylabel('');
        set(gca, 'fontsize', legend_size, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
        xlabel([]);
        xticks([]);
        
        %% Event psth and stats
        subplot(nplots_1, nplots_2,  [80,81,82]); hold on;
        title('PSTH', 'verticalalignment', 'top', 'fontsize', title_size);
        bar(SpikeDensity{ipart}.psth.(event_name).time, squeeze(nanmean(SpikeDensity{ipart}.psth.(event_name).trial(:, itemp, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
        set(gca, 'fontsize', legend_size, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
        %ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        
        % plot positive clusters
        lag = size(SpikeDensity{ipart}.psth.(event_name).time, 2) - size(SpikeDensity{ipart}.stat.(event_name){itemp}.time, 2);
        
        if isfield(SpikeDensity{ipart}.stat.(event_name){itemp}, 'posclusters')
            for ipos = 1 : size(SpikeDensity{ipart}.stat.(event_name){itemp}.posclusters, 2)
                if SpikeDensity{ipart}.stat.(event_name){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                    sel = find(SpikeDensity{ipart}.stat.(event_name){itemp}.posclusterslabelmat == ipos);
                    if length(sel) == 1
                        if sel == 1
                            continue
                        end
                        w = SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel) - SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel-1);
                    else
                        w = 1;
                    end
                    bar(SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel), SpikeDensity{ipart}.psth.(event_name).avg(itemp, sel+lag), w, 'facecolor', 'r', 'edgecolor', 'none');
                end
            end
        end
        
        % plot negative clusters
        if isfield(SpikeDensity{ipart}.stat.(event_name){itemp}, 'negclusters')
            for ineg = 1 : size(SpikeDensity{ipart}.stat.(event_name){itemp}.negclusters, 2)
                if SpikeDensity{ipart}.stat.(event_name){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                    sel = find(SpikeDensity{ipart}.stat.(event_name){itemp}.negclusterslabelmat == ineg);
                    if length(sel) == 1
                        if sel == 1
                            continue
                        end
                        w = SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel) - SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel-1);
                    else
                        w = 1;
                    end
                    bar(SpikeDensity{ipart}.stat.(event_name){itemp}.time(sel), SpikeDensity{ipart}.psth.(event_name).avg(itemp, sel+lag), w, 'facecolor', 'b', 'edgecolor', 'none');
                end
            end
        end
        axis tight
        
        % plot baseline patch
        x = [cfg.stats.bl.(event_name)(1) cfg.stats.bl.(event_name)(2) cfg.stats.bl.(event_name)(2) cfg.stats.bl.(event_name)(1)];
        ax = axis;
        y = [ax(3) ax(3) ax(4) ax(4)];
        patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
        set(gca, 'fontsize', legend_size,'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
        
        % plot baseline indicator line
        y = mean(SpikeDensity{ipart}.stat.(event_name){itemp}.baseline(:, itemp));
        plot([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)], [y, y], ':k');
        clear y
        
        % add LFP
        yyaxis right
        x = LFP_avg.time;
        y = LFP_avg.avg(chansel, :);
        plot(x,y,'-g');
        set(gca, 'fontsize', legend_size,'YColor','k', 'TickDir', 'out');
        %         ylabel('uV');
        yticks([]);
        xlim([cfg.epoch.toi.(event_name)(1), cfg.epoch.toi.(event_name)(2)]);
        x_length = cfg.epoch.toi.(event_name)(2) - cfg.epoch.toi.(event_name)(1);
        xtemp = x >= (x(end) - x_length/4);
        ytemp = mean(y(xtemp));
        text(x(end) - x_length/4, ytemp + abs(ytemp)/10, 'LFP', 'color', 'g', 'fontsize', 7, 'verticalalignment', 'bottom');
        
        
        %% spike waveform : plot average of all waveforms, and plot 1000 randomly selected trials
        subplot(nplots_1, nplots_2,  [63,70]); hold on;
        flip = -SpikeWaveforms_stats{ipart}.peak_direction(itemp);
        title('Action Potential Waveform', 'fontsize', title_size);
        if size(SpikeWaveforms{ipart}{itemp}.trial, 2) > 1000
            spikes_idx_sel = randperm(size(SpikeWaveforms{ipart}{itemp}.trial,2),1000);
        else
            spikes_idx_sel = 1:length(SpikeWaveforms{ipart}{itemp}.trial);
        end
        for itrial = spikes_idx_sel
            lh = plot(SpikeWaveforms{ipart}{itemp}.time{itrial}, SpikeWaveforms{ipart}{itemp}.trial{itrial}*flip, 'color', [0.3 0.3 0.3]);
            lh.Color = [0,0,0,0.1];
        end
        if isfield(SpikeWaveforms_stats{ipart}.waveformavg{itemp}, 'time')
            plot(SpikeWaveforms_stats{ipart}.waveformavg{itemp}.time, SpikeWaveforms_stats{ipart}.waveformavg{itemp}.avg*flip, 'color', [1 1 1 ], 'linewidth', 1);
            axis tight
            %ylim([min(SpikeWaveforms_stats{ipart}.waveformavg{itemp}.avg)*2, max(SpikeWaveforms_stats{ipart}.waveformavg{itemp}.avg)*2]);
        end
        xticklabels(xticks.*1000);
        xlabel('Time (ms)');
        ylabel('uV');
        %adapt y scale in case of aberrant raw spikes:
        yinf = min(SpikeWaveforms_stats{ipart}.waveformavg{itemp}.avg*flip)*2;
        ysup = max(SpikeWaveforms_stats{ipart}.waveformavg{itemp}.avg*flip)*10;
        y = ylim;
        if y(1) < yinf * 3
            ylim([yinf, y(2)]);
        end
        y = ylim;
        if y(2) > ysup * 3
            ylim([y(1), ysup]);
        end
        set(gca, 'fontsize', legend_size, 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
        
        
        %% ISI histogram
        subplot(nplots_1, nplots_2,  [77,84]); hold on;
        bar(SpikeStats{ipart}.(baseline_name){itemp}.isi_avg_time * 1000, SpikeStats{ipart}.(baseline_name){itemp}.isi_avg, 1, 'facecolor', [0 0 0], 'edgecolor', 'none');
        idx = SpikeStats{ipart}.(baseline_name){itemp}.isi_avg_time <= cfg.spike.RPV;
        bar(SpikeStats{ipart}.(baseline_name){itemp}.isi_avg_time(idx) * 1000, SpikeStats{ipart}.(baseline_name){itemp}.isi_avg(idx), 1, 'facecolor', [1 0 0], 'edgecolor', 'none');
        set(gca, 'fontsize', legend_size, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
        axis tight
        xlim([0, 0.025] * 1000);
        xlabel(sprintf('ISI (ms)\nRPV = %.2f%%', SpikeStats{ipart}.(baseline_name){itemp}.RPV * 100));
        ylabel('Spikecount');
        
        %% save figure
        temp      = split(unit_name, '_');
        temp_name = strcat(temp{1:2});
        channame  = temp{3};
        fname = fullfile(cfg.imagesavedir,'unit_overview',sprintf('%sp%d-%s_%s_overview_%s_%s', cfg.prefix, ipart,channame, temp_name, char(event_name), char(baseline_name)));
        if ~isfolder(fileparts(fname))
            fprintf('Creating directory %s\n', fileparts(fname));
            mkdir(fileparts(fname));
        end
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        fig.PaperType = 'A4';
        
        print(fig, '-dpdf', [fname, '.pdf'],'-r600');
        print(fig, '-dpng', [fname, '.png'],'-r300');
        
        close all
    end %itemp
end %ipart

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_events(do_plot, event_start, event_end, color)
if do_plot
    yplot = ylim;
    xplot = xlim;
    for event = 1 : size(event_start, 2)
        fill([event_start(event), event_end(event), event_end(event), event_start(event)], [yplot(1) yplot(1) yplot(2) yplot(2)], [-1 -1 -1 -1], 'facecolor', color, 'edgecolor', 'none', 'facealpha', 0.4, 'edgealpha', 0.4);
        %fill([event_start(event), event_end(event), event_end(event), event_start(event)], [yplot(1) yplot(1) yplot(2) yplot(2)], [-1 -1 -1 -1], 'facecolor', color, 'edgecolor', color, 'facealpha', 0.4, 'edgealpha', 0.1);
        %patch('XData', [event_start(event), event_end(event), event_end(event), event_start(event)], 'YData', [yplot(1) yplot(1) yplot(2) yplot(2)], 'ZData', [-1 -1 -1 -1], 'facecolor', color, 'edgecolor', color, 'facealpha', 0.4, 'edgealpha', 0.4);
    end
    ylim(yplot);
    xlim(xplot);
end
end
