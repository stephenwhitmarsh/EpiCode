function plot_patterns_multilevel(cfg)

% set LFP config to read requested patterns
LFP                     = readLFP(cfg);
TFR                     = TFRtrials(cfg);
SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg);
SpikeDensity_timelocked = spikeTrialDensity(cfg);

cfg_TFR                 = [];
cfg_TFR.channel         = 1;
cfg_TFR.colorbar        = 'no';
cfg_TFR.zlim            = 'maxabs';
cfg_TFR.title           = ' ';
cfg_TFR.baselinetype    = 'relchange';
cfg_TFR.interactive     = 'no';

ncols = max(3, size(cfg.LFP.name, 2));
nrows = 6;

for ipart = 1 : size(LFP, 2)
    
    fig = figure('visible', true);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'position', get(0,'ScreenSize'));
    set(fig, 'PaperOrientation', 'landscape');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    
    imarker = 1;
    
    % first print LFP and TFR
    for markername = string(fields(LFP{ipart}))'
        
        % LFP
        LFPavg = ft_timelockanalysis([], LFP{ipart}.(markername));
        subaxis(nrows, ncols, imarker, 'SpacingVert', 0.02); hold;
        
        n = 1; ytick = []; label = [];
        maxrange = max(max(abs(LFPavg.avg))) /2;
        for ichan = 1 : size(LFPavg.label, 1)
            ytick = [ytick, n*maxrange];
            x       = LFPavg.time;
            y       = LFPavg.avg(ichan, :);
            plot(x, y + n*maxrange, 'k');
            label{ichan} = LFPavg.label{ichan};
            n = n + 1;
        end
        yticks(ytick);
        set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabels', label, 'TickDir', 'out')
        axis tight;
        xlim(cfg.epoch.toi.(markername));
        
        % TFR
        h = subaxis(nrows, ncols, ncols + imarker, 'SpacingVert', 0.02); hold;
        
        cfg_TFR.figure      = h;
        cfg_TFR.baseline    = cfg.TFR.bl.(markername);
        cfg_TFR.ylim        = [1, 200];
        ft_singleplotTFR(cfg_TFR, TFR{ipart}.(markername));
        set(gca, 'XGrid', 'on', 'box', 'off',  'xticklabel', [], 'TickDir', 'out');
        
        ylabel('Frequency');
        c                   = colorbar;
        set(c, 'Location', 'southoutside', 'color', [0 0 0]);
        c.Title.String      = 'Relative change in power';
        pos                 = get(c, 'Position');
        pos(2)              = 0.03;
        set(c, 'pos', pos);
        xlim(cfg.epoch.toi.(markername));
        
        imarker = imarker + 1;
    end
    
    % now units
    
    if ~isempty(SpikeDensity_timelocked{ipart})
        nrunits = size(SpikeDensity_timelocked{ipart}.stat.(markername), 2);
    else
        nrunits = 0;
    end
    
    for iunit = 0 : nrunits
   
        imarker = 1;
        isubplot = 1;
        for markername = string(fields(LFP{ipart}))'
            
            % Firing rate
            s{isubplot} = subaxis(nrows, ncols, ncols * 2 + imarker, 'SpacingVert', 0.02); hold on;
            isubplot    = isubplot + 1;
            
            for iplot = 1 : size(SpikeDensity_timelocked{ipart}.psth.(markername).avg, 1)
                
                if iunit == 0 || iplot == iunit
                    alpha = 1;
                else
                    alpha = 0.2;
                end
                
                if SpikeDensity_timelocked{ipart}.psth.(markername).corr_pval(iplot) < 0.025
                    if SpikeDensity_timelocked{ipart}.psth.(markername).corr_rho(iplot) < 0
                        color  = [1, 0, 0];
                    else
                        color  = [0, 0, 1];
                    end
                else
                    color  = [0, 0, 0];
                end
                
                lh = plot(SpikeDensity_timelocked{ipart}.psth.(markername).time, SpikeDensity_timelocked{1}.psth.(markername).avg(iplot, :), 'linewidth', 1);
                lh.Color = [color, alpha];
            end
            
            ylabel('Firing rate (Hz)');
            xlim(cfg.epoch.toi.(markername));
            set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'TickDir', 'out');
            
            % plot raster
            s{isubplot} = subaxis(nrows, ncols, [ncols * 3 + imarker, ncols * 4 + imarker, ncols * 5 + imarker], 'SpacingVert', 0.02);
            isubplot    = isubplot + 1;
            
            cfg_raster              = [];
            cfg_raster.trialborders = 'no';
            if iunit == 0
                cfg_raster.spikechannel = find(SpikeDensity_timelocked{ipart}.psth.(markername).corr_pval(iplot) < 0.025);
            else
                cfg_raster.spikechannel = iunit;
            end
            
            %                 if size(SpikeTrials_timelocked{ipart}.(markername).trialinfo, 1) > 100
            %                     cfg_raster.trials = randi(size(SpikeTrials_timelocked{ipart}.(markername).trialinfo, 1), 1, 100);
            %                 else
            %                     cfg_raster.trials = 'all';
            %                 end
            ft_spike_plot_raster(cfg_raster, SpikeTrials_timelocked{ipart}.(markername));
            xlim(cfg.epoch.toi.(markername));
            
            set(gca, 'XGrid', 'on', 'box', 'off', 'TickDir', 'out', 'xlabel', []);
            %             lines_obj = findobj(gcf,'Type','Line');
            %             for i = 1:numel(lines_obj)
            %                 try
            %                     y = get(lines_obj(i), 'y');
            %                     height = y(2) - y(1);
            %                     y(1) = y(1) - height * 2;
            %                     y(2) = y(2) + height * 2;
            %                     set(lines_obj(i), 'y', y);
            %                 catch
            %                 end
            %             end
            
            
            imarker = imarker + 1;
            
        end % markername
        
        fname = fullfile(cfg.imagesavedir, 'rasterplots', [cfg.prefix, 'p', num2str(ipart), 'u', num2str(iunit), '_rasterplot'] );
        exportgraphics(fig, [fname, '.pdf']);
        exportgraphics(fig, [fname, '.tiff'], 'Resolution', 600);
        exportgraphics(fig, [fname, '.jpg'],  'Resolution', 600);
        
        % clear subplots
        for i = 1 : size(s, 2)
            cla(s{i});
        end
        
    end % iunit
    
end % ipart




disp('done');
