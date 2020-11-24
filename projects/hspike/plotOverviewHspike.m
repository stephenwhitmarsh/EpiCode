function plotOverviewHspike(cfg, marker, hypnogram, hypmusestat, SpikeTrials_timelocked, SpikeTrials_windowed, SpikeStats_windowed, SpikeDensity_timelocked, LFP, TFR, SpikeWaveforms)

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])
set(0,'defaultfigurecolor',[1 1 1]);

cfg.visible = ft_getopt(cfg, 'visible', 'on');

% colormap
cm          = cool(5);

% hypnogram labels to use
hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

% order of plots and colors 
labelorder  = ["REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

for ipart = 1 : size(cfg.directorylist,2)
    
    if isempty(SpikeTrials_windowed{ipart}.window)
        continue
    end
    
    SpikeTrials_windowed{ipart}.window.trialinfo.hyplabel(SpikeTrials_windowed{ipart}.window.trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
    
    for markername = unique(marker.name)'
        
        fprintf('Building overview for part %d marker %s\n', ipart, markername);
        
        LFP{ipart}.(markername).trialinfo.hyplabel(LFP{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
        SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
  
        % baseline correction
        cfgtemp             = [];
        cfgtemp.baseline    = 'yes';
        cfgtemp.baseline    = cfg.LFP.baselinewindow.(markername);
        LFP_avg.all         = ft_timelockbaseline(cfgtemp, LFP{ipart}.(markername));
        
        % select timeperiod and calculate variance
        cfgtemp             = [];
        cfgtemp.latency     = cfg.epoch.toi.(markername);
        LFP_avg.all         = ft_timelockanalysis(cfgtemp, LFP_avg.all);       
        [~, chansel]        = max(var(LFP_avg.all.avg, 0, 2));
        
        % prepare LFPs
        maxrange = 0;
        for hyplabel = hyplabels
            
            % baseline correction
            cfgtemp             = [];
            cfgtemp.baseline    = 'yes';
            cfgtemp.baseline    = cfg.LFP.baselinewindow.(markername);
            LFP_avg.(hyplabel)  = ft_timelockbaseline(cfgtemp,  LFP{ipart}.(markername));
            
            % select timeperiod and calculate variance
            cfgtemp             = [];
            cfgtemp.latency     = cfg.epoch.toi.(markername);
            cfgtemp.trials      = LFP{ipart}.(markername).trialinfo.hyplabel == hyplabel;
            LFP_avg.(hyplabel)  = ft_timelockanalysis(cfgtemp, LFP_avg.(hyplabel));
            
            % get scaling parameter for all channels
            for ichan = 1 : size(LFP_avg.(hyplabel).label, 1)
                y           = LFP_avg.(hyplabel).avg(ichan,:);
                ystd        = sqrt(LFP_avg.(hyplabel).var(ichan,:));
                maxrange    = max(abs([maxrange, y+ystd, y-ystd]));
            end
        end
        
        for itemp = 1 : size(SpikeTrials_timelocked{ipart}.(markername).label, 2)
            
            fig = figure('visible', cfg.visible);
            set(gcf,'position', get(0,'ScreenSize'));
            set(fig, 'PaperOrientation', 'landscape');
            set(fig, 'PaperUnits', 'normalized');
            set(fig, 'PaperPosition', [0 0 1 1]);

            %% LFP over all channels
            subplot(9,7, [43, 50, 57]); hold on;
            title('LFP macroelectrodes');
            n = 1;
            ytick = [];
            label = [];
            l = [];
            for ichan = 1 : size(LFP_avg.(hyplabel).label, 1)
                ytick = [ytick, n*maxrange];
                icolor = 1;
                for hyplabel = hyplabels(end:-1:1)
                    colloc = hyplabel == labelorder;
                    x = LFP_avg.(hyplabel).time;
                    y = LFP_avg.(hyplabel).avg(ichan, :);
                    ystd = sqrt(LFP_avg.(hyplabel).var(ichan,:));
                    patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)] + n*maxrange, cm(colloc,:), 'edgecolor', 'none', 'facealpha', 0.5);
                end
                for hyplabel = hyplabels(end:-1:1)
                    colloc = hyplabel == labelorder;
                    x = LFP_avg.(hyplabel).time;
                    y = LFP_avg.(hyplabel).avg(ichan, :);
                    plot(x,y + n*maxrange, 'color', cm(colloc,:),'linewidth', 0.5);
                    l{colloc} = sprintf('%s: n=%d rpt', hyplabel, sum(LFP_avg.(hyplabel).cfg.trials));
                end
                n = n + 1;
                label{ichan} = LFP_avg.(hyplabel).label{ichan}(2:end);
               
                if ichan == chansel
                    label{ichan} = [label{ichan}, '*'];
                end
            end
            yticks(ytick);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off',  'XGrid', 'on', 'yticklabels', label, 'TickDir', 'out')
            xlabel('Time (s)');
            axis tight;
            xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            hl = legend(p, l, 'location', 'southoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(2) = 0.02;
            set(hl, 'position', pos_legend);

            %% TFR & LFP
            
            try
                h = subplot(9,7, [44, 51, 58]); hold off;
                title(sprintf('TFR & LFP electrode %s', LFP_avg.(hyplabel).label{chansel}(2:end)),'Interpreter', 'none');
                
                cfgtemp                 = [];
                cfgtemp.channel         = 1;
                cfgtemp.baseline        = [-0.5 0.2];
                cfgtemp.baselinetype    = 'relchange';
                cfgtemp.ylim            = [10, 100];
                cfgtemp.zlim            = 'maxabs';
                cfgtemp.title           = '';
                cfgtemp.colorbar        = 'no';
                cfgtemp.figure          = h;
                cfgtemp.interactive     = 'no';
                ft_singleplotTFR(cfgtemp, TFR{ipart}.(markername));
                xlabel('Time (s)'); ylabel('Frequency');
                c = colorbar;
                set(c, 'Location', 'southoutside', 'color', [0 0 0]);
                c.Title.String = 'Relative change in power';
                pos = get(c, 'Position');
                pos(2) = 0.03;
                set(c, 'pos', pos);
                colormap(jet(1000));
                title(sprintf('TFR & LFP electrode %s', LFP_avg.(hyplabel).label{chansel}(2:end)),'Interpreter', 'none');
                
                % plot single average LFP
                yyaxis right
                set(gca,'YColor','k');
                hold on
                n = 1;
                ytick = [];
                label = [];
                l = [];
                [~, chansel] = max(var(LFP_avg.all.avg, 0, 2));
                x = LFP_avg.all.time;
                y = LFP_avg.all.avg(chansel, :);
                ystd = sqrt(LFP_avg.all.var(chansel, :));
                %             patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], [0 0 0], 'edgecolor', 'none', 'facealpha', 0.2);
                plot(x, y, 'color', [0 0 0],'linewidth', 0.5);
                h = ylabel('uV');
                %             h.Position(1) = 1.5;
                h.Rotation = -90;
                
                ymax = max(abs([y + ystd, y - ystd]));
                ylim([-ymax, ymax]);
                xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
                
                set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
            catch
            end
            
            %% Spike Density and stats
            subplot(9, 7, 53); hold on;
            title('Spike Density');
            bar(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).time, squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).trial(:, itemp, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
            ylabel('Firing rate (Hz)');

            % plot positive clusters
            lag = size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).time, 2) - size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.time, 2);
            
            if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'posclusters')
                for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);      
                        bar(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
                        
                        % plot percentage
                        %                             si = sel+lag;
                        %                             [Y, I] = max(sdf_bar_bar{ipart}.avg(sel+lag));
                        %                             x = sdf_bar_bar{ipart}.time(si(I));
                        %                             y = sdf_bar_bar{ipart}.avg(si(I));
                        %                             d = (sdf_bar_bar{ipart}.avg(si(I)) / stats_binned{ipart}.(markername).clusterstat{itemp}.bl.avg) * 100 - 100;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.maxcluster.perc{ipos} = d;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.maxcluster.x{ipos} = x;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.maxcluster.y{ipos} = y;
                        %                             text(x+0.05, y+0.25, sprintf('+%.1f%%\n', d), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    end
                end
            end
            
            % plot negative clusters
            if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'negclusters')
                for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                        bar(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
                        
                        %                             % plot percentage
                        %                             si = sel+lag;
                        %                             [Y, I] = min(sdf_bar{ipart}.avg(sel+lag));
                        %                             x = sdf_bar{ipart}.time(si(I));
                        %                             y = sdf_bar{ipart}.avg(si(I));
                        %                             d = 100 - (sdf_bar{ipart}.avg(si(I)) / stats_binned{ipart}.(markername).clusterstat{itemp}.bl.avg) * 100;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.mincluster.perc{ineg} = d;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.mincluster.x{ineg} = x;
                        %                             stats_binned{ipart}.(markername).clusterstat{itemp}.mincluster.y{ineg} = y;
                        %                             text(x+0.05, y-0.25, sprintf('-%.1f%%\n', d), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                    end
                end
            end
            
            % plot baseline patch
            x = [cfg.stats.bl.(markername)(1) cfg.stats.bl.(markername)(2) cfg.stats.bl.(markername)(2) cfg.stats.bl.(markername)(1)];
            ax = axis;
            y = [ax(3) ax(3) ax(4) ax(4)];
            patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
            set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out');

            % plot baseline indicator line
            y = mean(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.baseline(:, itemp));
            plot([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)], [y, y], ':k');
            clear y
            
            % add LFP 
            yyaxis right
            x = LFP_avg.all.time;
            y = LFP_avg.all.avg(chansel, :);
            plot(x,y,'-g');            
            set(gca,'YColor','k', 'TickDir', 'out');
            ylabel('uV');
            xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);

            %% Spike Density for each stage overlapping
            subplot(9, 7, 60); hold on;
            title('Spike Density x Sleep Stage');
            for hyplabel = hyplabels
                cindx               = hyplabel == labelorder;
                trials              = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                x                   = SpikeDensity_timelocked{ipart}.sdf_bar.(markername).time;
                y                   = squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).trial(trials, itemp, :), 1))';
                temp                = squeeze(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).trial(trials, itemp, :));
                ystd                = nanstd(temp);
                bar(x, y, 1, 'facecolor', cm(cindx,:), 'edgecolor', 'none', 'facealpha', 0.4);
                l{cindx} = sprintf('%s: n=%d rpt', hyplabel, sum(trials));
                
            end
            xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            xlabel('Time (s)'); ylabel('Firing rate (Hz)');
            set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out'); %, 'Xticklabels', []);
            
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            hl = legend(p, l, 'location', 'southoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(2) = 0.025;
            set(hl, 'position', pos_legend);
            
            %% Raster
            try
            subplot(9,7, [45, 52, 59]); % hold on;
            cfgtemp                 = [];
            cfgtemp.spikechannel    = itemp;
            cfgtemp.latency         = [cfg.stats.toi.(markername)(1), cfg.stats.toi.(markername)(2)];
            cfgtemp.trialborders    = 'yes';
            cfgtemp.linewidth       = 1;
            ft_spike_plot_raster(cfgtemp, SpikeTrials_timelocked{ipart}.(markername));
            set(gca, 'Yticklabels', [], 'Ytick', [], 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
            
            title('Rasterplot');

            lines = findobj(gcf,'Type','Line');
            for i = 1:numel(lines)
                try
                    y = get(lines(i), 'y');
                    height = y(2) - y(1);
                    y(1) = y(1) - height * 2;
                    y(2) = y(2) + height * 2;
                    set(lines(i), 'y', y);
                catch
                end
            end
            
            state = "xxx";
            startend = "end";
            clear x y
            y = [0 0];
            for itrial = 1 : size(SpikeTrials_timelocked{ipart}.(markername).trialinfo, 1)
                
                % detect change
                if ~strcmp(SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(itrial), state)
                    % if already passed "end", we are now at a new start
                    if strcmp(startend, "end")
                        startend = "start";
                        
                        % set current state as new state
                        state = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(itrial);
                        
                        % set starting point
                        y(1) = itrial;
                        
                        % set color corresponding to label
                        ci = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(itrial) == labelorder;
                        
                        % plot
                        x = [cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2), cfg.epoch.toi.(markername)(2), cfg.epoch.toi.(markername)(1)];
                        patch(x, [y(1) y(1) y(2) y(2)], cm(ci, :) ,'edgecolor', 'none', 'facealpha',0.4);
%                         fprintf('Start %s, y = %d', state, itrial);
                        
                        % if already passed "start", we are now at the end
                    else
                        startend = "end";
                        
                        % set current state as new state
                        state = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel(itrial);
                        
                        y(2) = itrial;
                        patch(x, [y(1) y(1) y(2) y(2)], cm(ci, :) ,'edgecolor', 'none', 'facealpha',0.4);
                        
%                         fprintf('; End %s, y = %d\n', state, itrial);
                        
                    end
                end
            end
            set(gca,'children',flipud(get(gca,'children')))
            ylabel('');
            xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');
            catch
            end
                
            %% ISI windowed
            subplot(9,7, 47); hold on;
            title('ISI');            
            bar(SpikeStats_windowed{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats_windowed{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0 0 0], 'edgecolor', 'none');
            set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');  
            axis tight
            xlim([0, 0.025] * 1000);
            xlabel('Time (ms)'); ylabel('Spikecount');
                       
            %% autocorr windowed
%             subplot(9,7, 54); hold on;
%             title('Autocorrelation');     
%             bar(SpikeStats_windowed{ipart}.window{itemp}.xcorr_time * 1000, squeeze(SpikeStats_windowed{ipart}.window{itemp}.xcorr(itemp, itemp, :)), 1, 'facecolor', [0 0 0], 'edgecolor', 'none');
%             axis tight;
%             ax = axis;
%             patch([-2 2 2 -2], [ax(3), ax(3), ax(4)*1.1, ax(4)*1.1], [1 0 0], 'edgecolor', 'none');
%             set(gca,'children',flipud(get(gca,'children')))
%             
%             xlim([-0.025, 0.025] * 1000);
%             set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');  
%             ylabel('Correlation');
     
            %% autocorr per sleepstage       
%             subplot(9,7, 61); hold on;
%             title('Autocorrelation x Sleep Stage');        
%             for hyplabel = hyplabels(end:-1:1)
%                 colloc = hyplabel == labelorder;
%                 x = SpikeStats_windowed{ipart}.window{itemp}.xcorr_time * 1000;
%                 y = squeeze(SpikeStats_windowed{ipart}.window{itemp}.(hyplabel).xcorr(itemp, itemp, :));
%                 bar(x, y, 1, 'facecolor', cm(colloc,:), 'edgecolor', 'none', 'facealpha', 0.5);
%             end
%             axis tight  
%             xlim([-0.025, 0.025]  * 1000);
%             set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out');  
%             xlabel('Time (ms)'); ylabel('Correlation');
            
            %% FR x IED rate

            % need to order groups to make sure the colormap is correct
            clear G
            for itrial = 1 : size(SpikeStats_windowed{ipart}.window{itemp}.trialinfo, 1)
                G(itrial) = find(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.hyplabel(itrial) == labelorder);
            end     
            [G, order]  = sort(G);
            spikefreq   = SpikeStats_windowed{ipart}.window{itemp}.trialfreq(order)';
            IEDfreq     = SpikeTrials_windowed{ipart}.window.trialinfo.(char(markername))(order);
            amplitude   = SpikeStats_windowed{ipart}.window{itemp}.amplitude(order);
            
            % plotting
            subplot(9,7, [55, 62]);
            title('Firing rate x IED rate');
            s = scatterhistogram(spikefreq, IEDfreq, 'GroupData', G, 'binwidth', [0.25; 1], 'HistogramDisplayStyle', 'bar', ...
                'Legendvisible', 'off', 'LineStyle', {'-','-','-','-','-'}, 'Color', cm, 'MarkerStyle', 'o', ...
                'MarkerSize', 10, 'MarkerAlpha', 0.5, 'ScatterPlotProportion', 0.5, ...
                'XLimits', [0; ceil(max(spikefreq))], ...
                'YLimits', [0; ceil(max(IEDfreq))]);
            xlabel('Firing rate (Hz)'); ylabel('IED rate (per minute)');
            pos1 = s.OuterPosition;
            
            newpos1 = pos1;
            set(s,'OuterPosition', newpos1);
            
            %% FR x Amplitude
            
            subplot(9,7, [56, 63]); 
            s = scatterhistogram(spikefreq, amplitude, 'GroupData', G, 'binwidth', [0.25; 0.25], 'HistogramDisplayStyle', 'bar', ...
                'Legendvisible', 'off', 'LineStyle', {'-','-','-','-','-'}, 'Color', cm, 'MarkerStyle', 'o', ...
                'MarkerSize', 10, 'MarkerAlpha', 0.5, 'ScatterPlotProportion', 0.7, ...
                'XLimits', [0; ceil(max(spikefreq))], ...
                'YLimits', [floor(min(amplitude)); ceil(max(amplitude))]);
            xlabel('Firing rate(Hz)'); ylabel('Amplitude (uV)');
            
            pos = get(s, 'Position');
            pos(1) = 0.875;
            set(s, 'Position', pos);
            
%             pos2 = s.OuterPosition;
%             newpos2 = pos2;
%             newpos2(3) = 0.3;            
%             newpos2(4) = 0.3;
%             set(s,'OuterPosition', newpos2);
            
            %% Hypnogram
            subplot(9,7,1:6); hold on;
            title(sprintf('Hypnogram Patient %s, night #%d', cfg.prefix(1:end-1), ipart));
            hypnogram_sel = hypnogram(hypnogram.part == ipart, :);
            plot_hyp_lines(hypnogram_sel);
            axis tight;
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            xl = xlim;
            
            %% Time per sleep stage
            subplot(9,7,7); hold on;
            title('Time per sleep stage (hrs.)');
            m = -inf;
            for hyplabel = hyplabels
                barloc      = find(hyplabel ==labelorder);
                y           = hypmusestat{ipart}.(char(markername)).duration.(hyplabel);
                hb          = bar(barloc, y, 1);
                l{barloc}   = sprintf('%s=%0.2fhrs', hyplabel, y);
                m           = max(m, y);
                set(hb, 'FaceColor', cm(barloc,:));
            end
            
            xlim([0.5, 5.5]); ylim([0, m * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
            
            %% IED rate over time
            subplot(9,7,8:13); hold on;
            title(sprintf('IED rate "%s"', char(markername)));
            marker_sel      = marker(marker.name == char(markername) & marker.ipart == ipart, :);
            [n, edges, ~]   = histcounts(marker_sel.clock, 'BinWidth', seconds(60));
            bar(edges(1:end-1), n, 1, 'EdgeColor', 'none', 'FaceColor', [0, 0, 0]);
            ylabel('IED / min.');
            axis tight
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickDir', 'out');
            xlim(xl);
            
            %% IED rate normalized to wake
            subplot(9,7,14); hold on;
            title('IED rate vs. wake (per minute)');
            m = -inf;
            for hyplabel = hyplabels
                barloc      = find(hyplabel ==labelorder);
                y           = hypmusestat{ipart}.(char(markername)).IEDrateNorm.(hyplabel);
                hb          = bar(barloc, y, 1);
                l{barloc}   = sprintf('%s=%0.2f', hyplabel, y);
                m           = max(m, y);
                set(hb, 'FaceColor', cm(barloc,:));
            end
            xlim([0.5, 5.5]); ylim([0, m * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');
            
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
            
            %% Firing rate over time
            subplot(9,7,15:20); hold on;
            title(sprintf('Firingrate Unit #%d', itemp));
            bar(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.trialfreq, 1, 'EdgeColor', 'none', 'FaceColor', [0, 0, 0]);
            plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.trialfreq_corrected, 'color', [1, 1, 1]);
            ylabel('Firing rate original (bar) and corrected (line)');
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            xlim(xl);
            set(gca,'TickLength',[0 0], 'Xticklabels', []);
 
            %% Average firing rates per stage
            subplot(9,7,21); hold on;
            title('Firing rate (Hz)');
            
            m = -inf;
            for hyplabel = hyplabels
                barloc      = find(hyplabel == labelorder);
                trials      = SpikeStats_windowed{ipart}.window{itemp}.trialinfo.hyplabel == hyplabel;
                y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.trialfreq(trials));
                stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.trialfreq(trials));
                hb          = bar(barloc, y, 1); set(hb, 'FaceColor', cm(barloc,:));
                l{barloc}   = sprintf('%s=%0.1f(%0.1f)p/s', hyplabel, y, stdev);
                m           = max([m, y, y+stdev]);
                he          = errorbar(barloc, y, stdev, 'clipping','off','color', 'k');
            end
            
            xlim([0.5, 5.5]); ylim([0, m * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');

            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
             
            %% Amplitude over time
            subplot(9,7,22:27); hold on;
            title(sprintf('Amplitude Unit #%d', itemp));
            plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.amplitude, 'Color', [0, 0, 0]);
            ylabel('Amplitude (uV)');
            axis tight
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            xlim(xl);
            set(gca,'TickLength',[0 0], 'Xticklabels', []);

            %% Average amplitude per stage
            subplot(9,7,28); hold on;
            title('Amplitude (uV)');
            
            ymax = -inf;
            ymin = +inf;
            for hyplabel = hyplabels
                barloc      = find(hyplabel == labelorder);
                trials      = SpikeStats_windowed{ipart}.window{itemp}.trialinfo.hyplabel == hyplabel;
                y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.amplitude(trials));
                stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.amplitude(trials));
                hb          = bar(barloc, y, 1); set(hb, 'FaceColor', cm(barloc,:));
                l{barloc}   = sprintf('%s=%0.1f(%0.1f)p/s', hyplabel, y, stdev);
                ymax        = max([ymax, y+stdev]);
                ymin        = min([ymin, y-stdev]);
                he          = errorbar(barloc, y, stdev, 'clipping','off','color', 'k');
            end
            
            xlim([0.5, 5.5]); ylim([ymin * 0.9, ymax * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
            
            %% CV2 over time
            subplot(9,7,29:34); hold on;
            title(sprintf('CV2 Unit #%d', itemp));
            plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.CV2_trial', 'Color', [0, 0, 0]);
            plot(xl,[1 1],'k--');
            ylabel('CV2');
            axis tight
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            plot(xl,[1 1],'k--');
            xlim(xl);
            set(gca,'TickLength',[0 0], 'Xticklabels', []);

            %% Average CV2
            subplot(9,7,35); hold on;
            title('CV2');
            
            ymax = -inf;
            ymin = +inf;
            for hyplabel = hyplabels
                barloc      = find(hyplabel == labelorder);
                trials      = SpikeStats_windowed{ipart}.window{itemp}.trialinfo.hyplabel == hyplabel;
                y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.CV2_trial(trials));
                stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.CV2_trial(trials));
                hb          = bar(barloc, y, 1); set(hb, 'FaceColor', cm(barloc,:));
                l{barloc}   = sprintf('%s=%0.1f(%0.1f)', hyplabel, y, stdev);
                ymax        = max([ymax, y+stdev]);
                ymin        = min([ymin, y-stdev]);
                he          = errorbar(barloc, y, stdev, 'clipping','off','color', 'k');
            end
            xlim([0.5, 5.5]); ylim([ymin * 0.9, ymax * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
            
            %% nr. of bursts over time
            subplot(9,7,36:41); hold on;
            title(sprintf('Number of bursts Unit #%d', itemp));
            plot(SpikeStats_windowed{ipart}.window{itemp}.trialinfo.starttime, SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum', 'Color', [0, 0, 0]);
            ylabel('Nr. of bursts per Minute');
            axis tight
            plot_hyp_colors(hypnogram_sel, cm, ylim);
            set(gca,'children',flipud(get(gca,'children')))
            xlim(xl);
            set(gca,'TickLength',[0 0], 'Xticklabels', []);

            %% Average nr. of bursts
            subplot(9,7,42); hold on;
            title('Bursts per Minute');
            
            m = 1;
            for hyplabel = hyplabels
                barloc      = find(hyplabel == labelorder);
                trials      = SpikeStats_windowed{ipart}.window{itemp}.trialinfo.hyplabel == hyplabel;
                y           = nanmean(SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum(trials));
                stdev       = nanstd(SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum(trials));
                hb          = bar(barloc, y, 1); set(hb, 'FaceColor', cm(barloc,:));
                l{barloc}   = sprintf('%s=%0.1f(%0.1f)', hyplabel, y, stdev);
                m           = max([m, y, y+stdev]);
                he          = errorbar(barloc, y, stdev, 'clipping','off','color', 'k');
            end
            xlim([0.5, 5.5]); ylim([0, m * 1.1]);
            set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off');
            clear p
            for ii = 1:size(cm,1)
                p(ii) = patch(NaN, NaN, cm(ii,:));
            end
            
            hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
            pos_legend = get(hl, 'position');
            pos_legend(1) = 0.92;
            set(hl, 'position', pos_legend);
            
            %% Spike waveform
            subplot(9, 7, 46); hold on;
            title(sprintf('Action Potential Waveform (%s)', SpikeTrials_windowed{ipart}.window.cluster_group{itemp}(1:3)));
            y = cat(1, SpikeWaveforms{ipart}.window{itemp}.trial{:});
            avg_waveshape = mean(y, 1);
            for itrial = 1 : size(SpikeWaveforms{ipart}.window{itemp}.trial, 2)
                lh = plot(SpikeWaveforms{ipart}.window{itemp}.time{itrial}*1000, SpikeWaveforms{ipart}.window{itemp}.trial{itrial}, 'color', [0.3 0.3 0.3]);
                lh.Color = [0,0,0,0.1];
            end
            plot(SpikeWaveforms{ipart}.window{itemp}.time{itrial}*1000, avg_waveshape, 'color', [1 1 1 ], 'linewidth', 1);
            axis tight
            ylim( [-max(abs(avg_waveshape)) max(abs(avg_waveshape))]*1.5);
            xlabel('Time (ms)'); ylabel('uV');
            set(gca, 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');

            %% print to file
            
            set(findall(fig, '-property', 'fontsize'), 'fontsize', 8);
            fname = fullfile(cfg.imagesavedir, strcat(cfg.prefix, sprintf('overview_night%d_spike%d_%s_%s', ipart, itemp, markername, SpikeTrials_windowed{ipart}.window.cluster_group{itemp}(1:3))));
            export_fig(fname, '-png'); % need to install https://www.ghostscript.com/download/gsdnld.html           
            close all
            
        end
    end
end

function plot_hyp_colors(h, cm, y)
% plot hypnogram in colors
for im = 1 : height(h)
    switch h.hyplabel{im}
        case 'NO_SCORE'
            ci = 1;        
        case 'REM'
            ci = 1;
        case 'AWAKE'
            ci = 2;
        case 'PHASE_1'
            ci = 3;
        case 'PHASE_2'
            ci = 4;
        case 'PHASE_3'
            ci = 5;
    end
    fill([h.starttime(im), h.endtime(im), h.endtime(im), h.starttime(im)],[y(1) y(1) y(2) y(2)], cm(ci, :) ,'edgecolor', 'none', 'facealpha',0.4);
end

function plot_hyp_lines(h)

% plot hypnogram lines
X = [];
Y = [];
for im = 1 : height(h)
    if ~isempty(X)
        % if there's a gap, 'fill' with 0
        if h.starttime(im) ~= X(end)
            X = [X, X(end) h.starttime(im)];
            Y = [Y, 0,  0];
        end
    end
    X = [X, h.starttime(im), h.endtime(im)];
    
    % height in hypnogram is based on order of cfg.hyp.contains
    switch cell2mat(h.hyplabel(im))
        case 'NO_SCORE'
            y = 5;        
        case 'AWAKE'
            y = 5;
        case 'REM'
            y = 4;
        case 'PHASE_1'
            y = 3;
        case 'PHASE_2'
            y = 2;
        case 'PHASE_3'
            y = 1;
    end
    Y = [Y, y, y];
end

for i = 1 : length(X)-1
    if Y(i) ~= 0 && Y(i+1) ~= 0
        if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
            plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
        else
            plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
        end
    end
end
set(gca,'Layer','top');
set(gca,'Ytick', 1 : 5, 'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'});
