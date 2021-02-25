function plotTimeCourses_seizures(cfg, LFP, TFR)

for ipart = 1 : size(LFP, 2)
    for markername = string(fieldnames(LFP{ipart}))'
        
        cfgtemp = [];
        cfgtemp.hpfilter = 'yes';
        cfgtemp.hpfreq = 1;
        LFP{ipart}.(markername) = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
    end
end

for ipart = 1 : size(LFP, 2)
    
    
    for markername = string(fieldnames(LFP{ipart}))'
        
        % select channels (and ignore missing)
        chanlist = unique(LFP{ipart}.(markername).trialinfo.chan1);
        
        for ichan = 1 : size(chanlist, 1)
            ind = contains(LFP{ipart}.(markername).trialinfo.chan1, chanlist{ichan});
%             contlist(ichan) = unique(LFP{ipart}.(markername).trialinfo.control(ind)); % so that it will break if there are more than 1 control channels
            contlist(ichan) = string(LFP{ipart}.(markername).trialinfo.control{find(ind,1,'first')}); % so that it will break if there are more than 1 control channels
            typelist(ichan) = string(LFP{ipart}.(markername).trialinfo.type{find(ind,1,'first')}); % so that it will break if there are more than 1 control channels

        end

        % get max range of timecourses in dataset
        cfgtemp         = [];
        cfgtemp.latency = [-0.5 10];
        sel             = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
        
        maxrange = 0;
        maxtrials = 0;
        for ichan = 1: size(chanlist, 1)
            trials      = find(strcmp(sel.trialinfo.chan1, chanlist{ichan}));
            maxtrials   = max([maxtrials, length(trials)]);
            channel     = strcmp(string(sel.label), strcat('_', chanlist{ichan}));
            for itrial = trials'
                y           = max(abs(sel.trial{itrial}(channel,:)));
                maxrange    = max(abs([maxrange, y]));
            end
        end
        
        maxrange = maxrange / 2; % extra scaling if needed
        ymax = maxrange * (maxtrials-1) + maxrange/2;
        ymin = -maxrange/2;
        
        % loop through channels
        for ichan = 1 : size(chanlist, 1)
            
            % select channel for timecourses
            cfgtemp         = [];
            cfgtemp.trials  = strcmp(LFP{ipart}.(markername).trialinfo.chan1, chanlist{ichan});
            cfgtemp.channel = find(strcmp(string(LFP{ipart}.(markername).label), strcat('_', chanlist{ichan})));
            LFP_chan_sel    = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
            
            cfgtemp.channel = find(strcmp(string(LFP{ipart}.(markername).label), strcat('_', contlist{ichan})));
            LFP_cont_sel    = ft_selectdata(cfgtemp, LFP{ipart}.(markername));

            % select channel and seizures for TFR
            cfgtemp                 = [];
            cfgtemp.trials          = strcmp(TFR{ipart}.(markername).trialinfo.chan1, chanlist{ichan});
            cfgtemp.channel         = find(strcmp(string(TFR{ipart}.(markername).label), strcat('_', chanlist{ichan})));
            cfgtemp.latency         = [-0.5 10];
            TFR_chan_sel            = ft_selectdata(cfgtemp, TFR{ipart}.(markername));

            cfgtemp.channel         = find(strcmp(string(TFR{ipart}.(markername).label), strcat('_', contlist{ichan})));
            cfgtemp.latency         = [-0.5 10];            
            TFR_cont_sel            = ft_selectdata(cfgtemp, TFR{ipart}.(markername));
            
            % baseline correction on indicidual seizures
            cfgtemp                 = [];
            cfgtemp.baseline        = [-0.5 0];
            cfgtemp.baselinetype    = 'relchange';
            TFR_chan_sel            = ft_freqbaseline(cfgtemp, TFR_chan_sel);
            TFR_cont_sel            = ft_freqbaseline(cfgtemp, TFR_cont_sel);
            
            % average over seizures
            cfgtemp                 = [];
            cfgtemp.avgoverrpt      = 'yes';
            TFR_chan_sel            = ft_selectdata(cfgtemp, TFR_chan_sel);
            TFR_cont_sel            = ft_selectdata(cfgtemp, TFR_cont_sel);
            
            % get maximum power to create equal colorbar between channels
            zmax                    = max(max(max(max(max(abs(TFR_chan_sel.powspctrm))))), max(max(max(max(abs(TFR_cont_sel.powspctrm))))));
            
            % get figure handle
            fig = figure;
            
            % plot timecourses channel
            subplot(3,2,[1, 3]); hold;
            
            n = 0;
            label = [];
            ytick = [];
      
            for itrial = 1 : size(LFP_chan_sel.trial, 2)
                ytick = [ytick, n*maxrange];
                x = LFP_chan_sel.time{itrial};
                y = LFP_chan_sel.trial{itrial};
                plot(x, y + n*maxrange, 'k');
                n = n + 1;
                label{itrial} = num2str(itrial);
            end % itrial
            
            yticks(ytick);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off',  'XGrid', 'on', 'yticklabels', label, 'TickDir', 'out')
            xlabel('Time (s)');
            ylabel('Seizure');
            axis tight;
            set(gca, 'YDir','reverse')
            title(strcat('Channel:'," ", chanlist{ichan}, '(', typelist{ichan}, ')'),'interpreter', 'none');            
            %             xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            xlim([-0.5, 10]);
            ylim([-maxrange, ymax]);
            
            top     = max(LFP_chan_sel.trial{end}) + (n-1)*maxrange;
            bottom  = min(LFP_chan_sel.trial{end});
            middle  = (top-bottom)/2;
            ymid    = (ymax-ymin) /2;
            ax      = axis;
            shift   = ymid - middle;            
            axis([ax(1), ax(2), ax(3)-shift, ax(4)-shift]);
            ax      = axis;
            plot([0, 0], [ax(3), ax(4)],'k:');
 
            % plot timecourses
            subplot(3, 2, [2, 4]); hold;
            
            n = 0;
            label = [];
            ytick = [];
            
            for itrial = 1 : size(LFP_cont_sel.trial, 2)
                ytick = [ytick, n*maxrange];
                x = LFP_cont_sel.time{itrial};
                y = LFP_cont_sel.trial{itrial};
                plot(x, y + n*maxrange, 'k');
                n = n + 1;
                label{itrial} = num2str(itrial);
            end % itrial
            
            yticks(ytick);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off',  'XGrid', 'on', 'yticklabels', label, 'TickDir', 'out')
            title(strcat('Channel:'," ", contlist{ichan}),'interpreter', 'none');

            xlabel('Time (s)');
            ylabel('Seizure');
            axis tight;
            set(gca, 'YDir','reverse')
            %             xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            xlim([-0.5, 10]);
            ylim([-maxrange, ymax]);
            
            top     = max(LFP_chan_sel.trial{end}) + (n-1)*maxrange;
            bottom  = min(LFP_chan_sel.trial{end});
            middle  = (top-bottom)/2;
            ymid    = (ymax-ymin) /2;
            ax      = axis;
            shift   = ymid - middle;            
            axis([ax(1), ax(2), ax(3)-shift, ax(4)-shift]);
            ax      = axis;
            plot([0, 0], [ax(3), ax(4)],'k:');
            
            % plot TFR channel
            h                       = subplot(3, 2, 5);
            cfgtemp                 = [];
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.interactive     = 'no';
            cfgtemp.colormap        = jet(1000);            
            cfgtemp.figure          = h;
            cfgtemp.title           = '';
            cfgtemp.zlim            = [-zmax, zmax];
            ft_singleplotTFR(cfgtemp, TFR_chan_sel);
            
            title('');
            xlabel('Time (s)'); ylabel('Frequency');
            c = colorbar;
            set(c, 'Location', 'southoutside', 'color', [0 0 0]);
            c.Title.String = 'Relative change in power';
            pos = get(c, 'Position');
            pos(2) = 0.03;
            set(c, 'pos', pos);
            xlim([-0.5, 10]);
            set(findall(fig, '-property', 'fontsize'), 'fontsize', 7);
    
            h                       = subplot(3, 2, 6);
            cfgtemp                 = [];
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.interactive     = 'no';
            cfgtemp.figure          = h;
            cfgtemp.colormap        = jet(1000);
            cfgtemp.title           = '';
            cfgtemp.zlim            = [-zmax, zmax];            
            ft_singleplotTFR(cfgtemp, TFR_cont_sel);
            
            title('');
            xlabel('Time (s)'); ylabel('Frequency');
            c = colorbar;
            set(c, 'Location', 'southoutside', 'color', [0 0 0]);
            c.Title.String = 'Relative change in power';
            pos = get(c, 'Position');
            pos(2) = 0.03;
            set(c, 'pos', pos);
            xlim([-0.5, 10]);
            set(findall(fig, '-property', 'fontsize'), 'fontsize', 7);            
            
            % print to file
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'channel_' ,chanlist{ichan},'_bipolar')), '-r300');
            print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'channel_' ,chanlist{ichan},'_bipolar')), '-r300');
            
            close all
        end % ichan
    end % markername
end % ipart
