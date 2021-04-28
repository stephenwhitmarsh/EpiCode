function plot_patterns_multilevel_examples(cfg)

% cfg.plot.unit       = ft_getopt(cfg.plot, 'unit',  0); % which means all
% cfg.plot.trial      = ft_getopt(cfg.plot, 'trial', 1:5); % 0 means all
% cfg.plot.dir        = ft_getopt(cfg.plot, 'dir',   0); % which means all


% set LFP config to read requested patterns
LFP = readLFP(cfg);
for ipart = 1 : size(LFP, 2)
    for markername = string(fields(LFP{ipart}))'
        LFPavg{ipart}.(markername) = ft_timelockanalysis([], LFP{ipart}.(markername));
    end
end

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

for ipart = 1 : size(LFPavg, 2)
    
    fig = figure('visible', true);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'position', get(0,'ScreenSize'));
    set(fig, 'PaperOrientation', 'portrait');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer', 'Painters');
    
    imarker = 1;
    
    % first print LFP and TFR
    for markername = string(cfg.LFP.name)
        
        % LFP average
        subaxis(nrows, ncols, imarker, 'SpacingVert', 0.01, 'SpacingHoriz', 0.01); hold;
        n = 1; ytick = []; label = [];
        maxrange = max(max(abs(LFPavg{ipart}.(markername).avg))) /2;
        for ichan = 1 : size(LFPavg{ipart}.(markername).label, 1)
            ytick = [ytick, n*maxrange];
            x       = LFPavg{ipart}.(markername).time;
            y       = LFPavg{ipart}.(markername).avg(ichan, :);
            %                 if ichan == SpikeDensity_timelocked{ipart}.psth.rho_chan
            %                     plot(x, y + n*maxrange, 'r');
            %                 else
            plot(x, y + n*maxrange, 'k');
            %                 end
            label{ichan} = LFPavg{ipart}.(markername).label{ichan};
            n = n + 1;
        end
        yticks(ytick);
        if imarker == 1
            set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabels', label, 'TickDir', 'out')
        else
            set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabel', [], 'TickDir', 'out')
        end
        
        axis tight;
        xlim(cfg.epoch.toi.(markername));
        
        % TFR average
        h = subaxis(nrows, ncols, ncols + imarker, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
        cfg_TFR.figure      = h;
        cfg_TFR.baseline    = cfg.TFR.bl.(markername);
        cfg_TFR.ylim        = [1, 200];
        ft_singleplotTFR(cfg_TFR, TFR{ipart}.(markername));
        set(gca, 'XGrid', 'on', 'box', 'off',  'xticklabel', [], 'TickDir', 'out');
        
        if imarker == 1
            ylabel('Frequency');
        else
            set(gca, 'yticklabel', [], 'TickDir', 'out')
        end
        c = colorbar; set(c, 'Location', 'southoutside', 'color', [0 0 0]);
        c.Title.String = 'Relative change in power';
        pos = get(c, 'Position'); pos(2) = 0.03; set(c, 'pos', pos);
        xlim(cfg.epoch.toi.(markername));
        imarker = imarker + 1;
    end
    
    imarker_lfp = 1;

    for itrial = 1 : size(cfg.plot.trial.(markername))
        
        for markername_trial = string(cfg.LFP.name)

            % LFP example
            cfgtemp = [];
            cfgtemp.trials = find(LFP{ipart}.(markername_trial).trialinfo.idir == cfg.plot.dir.(markername_trial)(itrial) ...
                & LFP{ipart}.(markername_trial).trialinfo.trialnr == cfg.plot.trial.(markername_trial)(itrial));
            
            if isempty(cfgtemp.trials)
                disp('*** No such trials! ***')
                continue
            end         
            
            LFP_example = ft_selectdata(cfgtemp, LFP{ipart}.(markername_trial));
            
            s_lfp{imarker_lfp} = subaxis(nrows, ncols, ncols * 2 + imarker_lfp, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
            
            n = 1; ytick = []; label = [];
            maxrange = max(max(abs(LFP_example.trial{1}))) /2;
            for ichan = 1 : size(LFP_example.label, 1)
                ytick = [ytick, n*maxrange];
                x       = LFP_example.time{1};
                y       = LFP_example.trial{1}(ichan, :);
                plot(x, y + n*maxrange, 'k');
                label{ichan} = LFP_example.label{ichan};
                n = n + 1;
            end
            yticks(ytick);
            if imarker_lfp == 1
                set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabels', label, 'TickDir', 'out')
            else
                set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabel', [], 'TickDir', 'out')
            end
            axis tight;
            xlim(cfg.epoch.toi.(markername_trial));
            
            imarker_lfp = imarker_lfp + 1;
        end
        
        for markername_trial = string(cfg.LFP.name)

            if isempty(SpikeDensity_timelocked{ipart})
                continue
            end
            
            if cfg.plot.unit.(markername_trial) == 0
                units = 1 : size(SpikeTrials_timelocked{ipart}.(markername_trial).label, 2);
            end
            
            for iunit = units
                
                imarker_unit = 1;
                isubplot = 1;
                
                for markername_unit = string(fields(LFP{ipart}))'
                    
                    % Example raw MUA data
                    
                    trialindx = find(SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.idir == cfg.plot.dir.(markername_unit)(itrial) ...
                        & SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.trialnr_dir == cfg.plot.trial.(markername_unit)(itrial));
                    
                    if isempty(trialindx)
                        continue
                    end
                    
                    s{isubplot} = subaxis(nrows, ncols, ncols * 3 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold on;
                    isubplot    = isubplot + 1;
                    
                    ichan = 1;
                    for chan = string(cfg.circus.channel)
                        temp                        = dir(fullfile(cfg.rawdir, SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.directory(trialindx, :), strcat('*', chan, ".ncs")));
                        cfgtemp                     = [];
                        cfgtemp.dataset             = fullfile(cfg.rawdir, SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.directory(trialindx, :), temp.name);
                        cfgtemp.trl(1)              = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.begsample(trialindx) - SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.fileoffset(trialindx);
                        cfgtemp.trl(2)              = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.endsample(trialindx) - SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.fileoffset(trialindx);
                        cfgtemp.trl(3)              = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.offset(trialindx);
                        cfgtemp.hpfilter            = 'yes';
                        cfgtemp.hpfreq              = 500;
                        cfgtemp.hpfiltord           = 3;
                        dat_chan{ichan}             = ft_preprocessing(cfgtemp);
                        dat_chan{ichan}.label{1}    = char(chan);
                        ichan = ichan + 1;
                    end
                    dat_MUA = ft_appenddata([], dat_chan{:}); clear dat_chan
                    
                    n = 1; ytick = []; label = [];
                    maxrange = max(max(abs(dat_MUA.trial{1}))) * 2;
                    for ichan = 1 : size(dat_MUA.label, 1)
                        ytick = [ytick, n*maxrange];
                        x       = dat_MUA.time{1};
                        y       = dat_MUA.trial{1}(ichan, :);
                        plot(x, y + n*maxrange, 'color', [0.5, 0.5, 0.5]);
                        label{ichan} = dat_MUA.label{ichan};
                        n = n + 1;
                    end
                    
                    cm = turbo(length(units));
                    ci = 1;
                    for iiunit = units
                        spikeidx        = SpikeTrials_timelocked{ipart}.(markername_unit).trial{iiunit} == SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.trialnr(trialindx);
                        spiketime       = SpikeTrials_timelocked{ipart}.(markername_unit).time{iiunit}(spikeidx);
                        [~, chanindx]   = max(rms(permute(SpikeTrials_timelocked{ipart}.(markername_unit).template{iiunit}(1,:,:), [2, 3, 1]), 2));
 
                        for ispike = 1 : size(spiketime, 2)
                            t1 = spiketime(ispike) - 0.001;
                            t2 = spiketime(ispike) + 0.001;
                            sel = dat_MUA.time{1} >= t1 & dat_MUA.time{1} <= t2;
                            if any(sel)
                                
                                plot(dat_MUA.time{1}(sel), dat_MUA.trial{1}(chanindx, sel) + chanindx*maxrange, 'color', cm(ci, :));
                                y = max(dat_MUA.trial{1}(chanindx, sel)) * 1.2;
                                if iiunit == iunit
                                    plot(spiketime(ispike), y + chanindx*maxrange, 'v', 'markerfacecolor', [1 0 0], 'markeredgecolor', [1 0 0], 'markersize', 2);
                                end
                            end
                        end
                        ci = ci + 1;
                        
                    end
                    yticks(ytick);
                    if imarker_unit == 1
                        set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabels', label, 'TickDir', 'out')
                    else
                        set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'yticklabel', [], 'TickDir', 'out')
                    end
                    axis tight;
                    xlim(cfg.epoch.toi.(markername_unit));
                    
                    clear dat_MUA
                    
                    % plot raster
                    s{isubplot}             = subaxis(nrows, ncols, ncols * 4 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01);
                    isubplot                = isubplot + 1;
                    cfg_raster              = [];
                    cfg_raster.trialborders = 'no';
                    if iunit == 0
                        cfg_raster.spikechannel = 'all';
                    else
                        cfg_raster.spikechannel = iunit;
                    end
                    ft_spike_plot_raster(cfg_raster, SpikeTrials_timelocked{ipart}.(markername_unit));
                    xlim(cfg.epoch.toi.(markername_unit));
                    %                 set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'TickDir', 'out', 'xlabel', []);
                    if imarker_unit == 1
                        set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'TickDir', 'out', 'xlabel', []);
                    else
                        set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'xlabel', [], 'yticklabel', [], 'ylabel', [], 'TickDir', 'out')
                    end
                    
                    % Firing rate
                    s{isubplot} = subaxis(nrows, ncols, ncols * 5 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold on;
                    isubplot    = isubplot + 1;
                    
                    for iplot = 1 : size(SpikeDensity_timelocked{ipart}.psth.(markername_unit).avg, 1)
                        
                        if iunit == 0 || iplot == iunit
                            alpha = 1;
                        else
                            alpha = 0.2;
                        end
                        
                        if SpikeDensity_timelocked{ipart}.psth.(markername_unit).corr_pval(iplot) < 0.025
                            if SpikeDensity_timelocked{ipart}.psth.(markername_unit).corr_rho(iplot) < 0
                                color  = [1, 0, 0];
                            else
                                color  = [0, 0, 1];
                            end
                        else
                            color  = [0, 0, 0];
                        end
                        
                        lh = plot(SpikeDensity_timelocked{ipart}.psth.(markername_unit).time, SpikeDensity_timelocked{1}.psth.(markername_unit).avg(iplot, :), 'linewidth', 1);
                        lh.Color = [color, alpha];
                    end
                    
                    ylabel('Firing rate (Hz)');
                    xlim(cfg.epoch.toi.(markername_unit));
                    if imarker_unit == 1
                        set(gca, 'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
                    else
                        set(gca,'XGrid', 'on', 'box', 'off', 'yticklabel', [], 'ylabel', [], 'TickDir', 'out')
                    end
                    xlabel('Time (s)');
                    imarker_unit     = imarker_unit + 1;
                    
                end % markername
                
                fname = fullfile(cfg.imagesavedir, 'rasterplots', [cfg.prefix, 'p', num2str(ipart), '_dir', num2str(cfg.plot.dir.(markername_unit)(itrial)), '_unit', num2str(iunit), '_(', strtrim(SpikeTrials_timelocked{ipart}.(markername_unit).cluster_group{iunit}) ,')', '_trial', num2str(cfg.plot.trial.(markername_unit)(itrial)), '_rasterplot'] );
%                 exportgraphics(fig, [fname, '.pdf']);
                %             exportgraphics(fig, [fname, '.tiff'], 'Resolution', 600);
                exportgraphics(fig, [fname, '.jpg'],  'Resolution', 600);
                
                % clear unit subplots
                for i = 1 : size(s, 2)
                    cla(s{i});
                end
                
            end % iunit

        end %imarker
        
        % clear unit subplots
        for i = 1 : size(s_lfp, 2)
            cla(s_lfp{i});
        end
        
    end %itrial
    
end % ipart

disp('done');
