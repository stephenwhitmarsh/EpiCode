function plot_patterns_multilevel_examples2(cfg)

cfg.plot.reref      = ft_getopt(cfg.plot, 'reref', 'no');
cfg.plot.refmethod  = ft_getopt(cfg.plot, 'refmethod', 'bipolar');
cfg.plot.postfix    = ft_getopt(cfg.plot, 'postfix', []);
cfg.spike.name      = cfg.plot.name;
cfg.LFP.name        = cfg.plot.name;
cfg.TFR.name        = cfg.plot.name;

% set LFP config to read requested patterns
LFP = readLFP(cfg);
for ipart = 1 : size(LFP, 2)
    for markername = string(fields(LFP{ipart}))'
        try
            LFPavg{ipart}.(markername) = ft_timelockanalysis([], LFP{ipart}.(markername));
        catch
        end
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

% ncols = max(3, size(cfg.LFP.name, 2));
ncols   = cfg.plot.ncols;
nrows   = 6;

for ipart = 1 : size(LFPavg, 2)
    
    fig = figure('visible', true);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'position', get(0,'ScreenSize'));
    % set(0, 'DefaultFigurePosition', [200 300  1000 500]);
    set(fig, 'PaperOrientation', 'portrait');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer', 'Painters');
    
    imarker = 1;
    
    % first print average LFP and TFR
    for markername = string(cfg.plot.name)
        
        if ~isfield(LFPavg{ipart}, markername)
            continue
        end
        
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
        title(markername);
        
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
    
    
    % loop over trials
    for itrial = 1 : length(cfg.plot.trial.(markername){ipart})
        
        % plot example trial for each marker
        imarker_lfp = 1;
        
        for markername_trial = string(cfg.plot.name)   
            
            if isempty(SpikeTrials_timelocked{ipart}.(markername_trial))
                continue
            end
            
            indx        = cfg.plot.trial.(markername_trial){ipart}(itrial);
            
            if indx > height(SpikeTrials_timelocked{ipart}.(markername_trial).trialinfo)
                fprintf('You asked for a higher trialcount for %s than possible (%d from max %d)\n', markername_trial, indx, height(SpikeTrials_timelocked{ipart}.(markername_trial).trialinfo));
                continue
            end
            
            directory   = SpikeTrials_timelocked{ipart}.(markername_trial).trialinfo.directory(indx, :);
            
            temp        = dir(fullfile(cfg.rawdir, directory, ['*', cfg.LFP.channel{1}, '.ncs']));
            dataset     = fullfile(cfg.rawdir, directory, temp.name);
            hdr         = ft_read_header(dataset);

            spikeFs     = SpikeDensity_timelocked{ipart}.psth.(markername_trial).cfg.previous.hdr.Fs;
            ss          = round(SpikeTrials_timelocked{ipart}.(markername_trial).trialinfo.t0(indx) / spikeFs * hdr.Fs);
            es          = round(SpikeTrials_timelocked{ipart}.(markername_trial).trialinfo.t1(indx) / spikeFs * hdr.Fs);
            
            Startsample = ss + cfg.epoch.toi.(markername_trial)(1) * hdr.Fs - cfg.epoch.pad.(markername_trial) * hdr.Fs;
            Endsample   = es + cfg.epoch.toi.(markername_trial)(2) * hdr.Fs + cfg.epoch.pad.(markername_trial) * hdr.Fs;
            Offset      = (cfg.epoch.toi.(markername_trial)(1) - cfg.epoch.pad.(markername_trial)) * hdr.Fs;
   
            for ifile = 1 : size(cfg.LFP.channel, 2)
                 
                temp = dir(fullfile(cfg.rawdir, directory, ['*', cfg.LFP.channel{ifile}, '.ncs']));
                if isempty(temp)
                    fprintf('Could not find %s\n', cfg.LFP.channel{ifile});
                    continue
                end
                 
                cfgtemp             = [];
                cfgtemp.dataset     = fullfile(cfg.rawdir, directory, temp.name);
                cfgtemp.trl         = round([Startsample; Endsample; Offset]');
                filedat{ifile}      = ft_preprocessing(cfgtemp);
            end
            
            cfgtemp                = [];
            cfgtemp.keepsampleinfo = 'no';
            LFP_example            = ft_appenddata(cfgtemp, filedat{:});
            clear filedat*
            
            % bipolar rereferencing if requested
            if strcmp(cfg.plot.reref, 'yes')
                if strcmp(cfg.plot.refmethod, 'bipolar')
                    labels_nonum    = regexprep(LFP_example.label, '[_0-9]', '_');
                    [~,~,indx]      = unique(labels_nonum);
                    clear group
                    for i = 1 : max(indx)
                        cfgtemp             = [];
                        cfgtemp.reref       = 'yes';
                        cfgtemp.refmethod   = 'bipolar';
                        cfgtemp.channel     = LFP_example.label(indx==i);
                        group{i}            = ft_preprocessing(cfgtemp, LFP_example);
                    end
                    LFP_example = ft_appenddata([], group{:});
                    clear group
                end
            end
            
            s_lfp{imarker_lfp} = subaxis(nrows, ncols, ncols * 2 + imarker_lfp, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01);
            hold on;
            
            n = 1; ytick = []; label = [];
            maxrange = max(max(abs(LFP_example.trial{1}))) / 2;
            for ichan = 1 : size(LFP_example.label, 1)
                ytick = [ytick, n * maxrange];
                x       = LFP_example.time{1};
                y       = LFP_example.trial{1}(ichan, :);
                plot(x, y + n * maxrange, 'k');
                label{ichan} = LFP_example.label{ichan}(end-6:end);
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
            title(sprintf("Trial %d (%s)", cfg.plot.trial.(markername_trial){ipart}(itrial), directory), 'interpreter', 'none');
            imarker_lfp = imarker_lfp + 1;
        end % markername
        
        % continue with units if present
        if isempty(SpikeDensity_timelocked{ipart})
            continue
        end
        
        % -1 means all units individually; 0 means all at once
        if cfg.plot.unit{ipart} == -1
            units = 1 : size(SpikeTrials_timelocked{ipart}.(markername_trial).label, 2);
        else
            units = cfg.plot.unit{ipart}(itrial);
        end
        
        % loop over units
        for iunit = units
            
            imarker_unit = 1;
            isubplot     = 1;
            
            % loop over markers again
            for markername_unit = string(cfg.plot.name)
                
                trialindx = cfg.plot.trial.(markername_unit){ipart}(itrial);

                % plot traces with units
                s{isubplot}         = subaxis(nrows, ncols, ncols * 3 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold on;
                isubplot            = isubplot + 1;
                ichan               = 1;
                
                if ~isfield(SpikeTrials_timelocked{ipart}, markername_unit)
                    continue
                end
                if isempty(SpikeTrials_timelocked{ipart}.(markername_unit))
                    continue
                end                
                if trialindx > height(SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo)
                    fprintf('You asked for a higher trialcount for %s than possible (%d from max %d)\n', markername_unit, trialindx, height(SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo));
                    continue
                end
                dirname             = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.directory(trialindx, :);
                
                cfgtemp             = [];
                cfgtemp.trl(1)      = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.begsample(trialindx) - SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.fileoffset(trialindx);
                cfgtemp.trl(2)      = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.endsample(trialindx) - SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.fileoffset(trialindx);
                cfgtemp.trl(3)      = SpikeTrials_timelocked{ipart}.(markername_unit).trialinfo.offset(trialindx);
                cfgtemp.hpfilter    = 'yes';
                cfgtemp.hpfreq      = 300;
                cfgtemp.hpfiltord   = 3;
                for chan = string(cfg.circus.channel)
                    temp                        = dir(fullfile(cfg.rawdir, dirname, strcat('*', chan, ".ncs")));
                    cfgtemp.dataset             = fullfile(cfg.rawdir, dirname, temp.name);
                    dat_chan{ichan}             = ft_preprocessing(cfgtemp);
                    dat_chan{ichan}.label{1}    = char(chan);
                    ichan = ichan + 1;
                end
                dat_MUA = ft_appenddata([], dat_chan{:}); clear dat_chan
                
                n = 1; ytick = []; label = [];
                maxrange = max(max(abs(dat_MUA.trial{1})));
                for ichan = 1 : size(dat_MUA.label, 1)
                    ytick = [ytick, n*maxrange];
                    x       = dat_MUA.time{1};
                    y       = dat_MUA.trial{1}(ichan, :);
                    plot(x, y + n*maxrange, 'color', [0.5, 0.5, 0.5]);
                    label{ichan} = dat_MUA.label{ichan};
                    n = n + 1;
                end

                if units == 0
                    units2color = 1 : length(SpikeTrials_timelocked{ipart}.(markername_unit).label);
                else
                    units2color = units;
                end
                
                cm = turbo(length(SpikeTrials_timelocked{ipart}.(markername_unit).label));
                ci = 1;
                
                if strcmp(markername_unit, "SEIZURE")
                    disp('wait');
                end
                                
                for iiunit = units2color
                    
                    spikeidx        = SpikeTrials_timelocked{ipart}.(markername_unit).trial{iiunit} == trialindx;
                    spiketime       = SpikeTrials_timelocked{ipart}.(markername_unit).time{iiunit}(spikeidx);
                    [~, chanindx]   = max(rms(permute(SpikeTrials_timelocked{ipart}.(markername_unit).template{iiunit}(1,:,:), [2, 3, 1]), 2));
                    
                    for ispike = 1 : size(spiketime, 2)
                        t1 = spiketime(ispike) - 0.001;
                        t2 = spiketime(ispike) + 0.001;
                        sel = dat_MUA.time{1} >= t1 & dat_MUA.time{1} <= t2;
                        if any(sel)
                            
                            plot(dat_MUA.time{1}(sel), dat_MUA.trial{1}(chanindx, sel) + chanindx * maxrange, 'color', cm(ci, :));
                            y = max(dat_MUA.trial{1}(chanindx, sel)) * 1.2;
                            if units > 0 & iiunit == iunit
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
                if iunit > 0 
                    title(sprintf('Unit %d', iunit));
                else
                    title('All units');
                end
                clear dat_MUA
                
                % plot raster
                s{isubplot}             = subaxis(nrows, ncols, ncols * 4 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01);
                isubplot                = isubplot + 1;
                cfg_raster              = [];
                cfg_raster.trialborders = 'no';
                cfg_raster.cmapneurons  = cm;
                if iunit == 0
                    cfg_raster.spikechannel = 'all';
                else
                    cfg_raster.spikechannel = iunit;
                end
                ft_spike_plot_raster(cfg_raster, SpikeTrials_timelocked{ipart}.(markername_unit));
                xlim(cfg.epoch.toi.(markername_unit));
                if imarker_unit == 1
                    set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'TickDir', 'out', 'xlabel', []);
                else
                    %                     set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'xlabel', [], 'yticklabel', [], 'ylabel', [], 'TickDir', 'out')
                    set(gca, 'XGrid', 'on', 'box', 'off', 'xticklabel', [], 'xlabel', [], 'yticklabel', [], 'TickDir', 'out')
                end
                
                % Firing rate
                s{isubplot} = subaxis(nrows, ncols, ncols * 5 + imarker_unit, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold on;
                isubplot    = isubplot + 1;
                
                if isfield(SpikeDensity_timelocked{ipart}.psth.(markername_unit), 'corr_pval')
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
                        lh = plot(SpikeDensity_timelocked{ipart}.psth.(markername_unit).time, SpikeDensity_timelocked{ipart}.psth.(markername_unit).avg(iplot, :), 'linewidth', 1);
                        lh.Color = [color, alpha];
                    end
                end
                
                ylabel('Firing rate (Hz)');
                xlim(cfg.epoch.toi.(markername_unit));
                if imarker_unit == 1
                    set(gca, 'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
                else
                    set(gca,'XGrid', 'on', 'box', 'off', 'yticklabel', [], 'ylabel', [], 'TickDir', 'out')
                end
                xlabel('Time (s)');
                imarker_unit = imarker_unit + 1;
                
            end % markername
            
            if iunit == 0
                fname = fullfile(cfg.imagesavedir, 'rasterplots', strcat(cfg.prefix, 'p', num2str(ipart), '_indx', num2str(itrial), '_allunits', cfg.plot.postfix));
            else
                fname = fullfile(cfg.imagesavedir, 'rasterplots', strcat(cfg.prefix, 'p', num2str(ipart), '_indx', num2str(itrial), '_unit', num2str(iunit), '_(', strtrim(SpikeTrials_timelocked{ipart}.(markername_unit).cluster_group{iunit}),')', cfg.plot.postfix));
            end
            exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 600);
            % exportgraphics(fig, [fname, '.pdf']);
            
            % clear unit subplots
            for i = 1 : size(s, 2)
                cla(s{i});
            end
            
        end % iunit
        
        % clear LFP subplots
        for i = 1 : size(s_lfp, 2)
            cla(s_lfp{i});
        end
        
    end %itrial
end % ipart

disp('done');
