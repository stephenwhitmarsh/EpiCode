function plotTimeCourses_new(cfg, LFP, TFR)

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
   
        % get scaling parameter for all channels
        maxrange = 0;
        for ichan = 1 : size(LFP{ipart}.(markername).label, 1)
            for itrial = 1 : size(LFP{ipart}.(markername).trial, 1)
                y           = LFP{ipart}.(markername).trial{itrial}(ichan,:);
                maxrange    = max(abs([maxrange, y]));
            end
        end
%         maxrange = maxrange / 4;
        
        for ichan = 1 : size(LFP{ipart}.(markername).label, 1)
            
            fig = figure; 
            subplot(2,1,1); hold;
            
            n = 0;
            label = [];
            ytick = [];

            for itrial = 1 : size(LFP{ipart}.(markername).trial, 2)
                title(['Seizures, electrode: ', LFP{ipart}.(markername).label{ichan}],'interpreter', 'none');
                ytick = [ytick, n*maxrange];
                x = LFP{ipart}.(markername).time{itrial};
                y = LFP{ipart}.(markername).trial{itrial}(ichan, :);
                plot(x, y + n*maxrange, 'k');
                
                n = n + 1;
                label{itrial} = num2str(itrial);
            end % itrial

            yticks(ytick);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off',  'XGrid', 'on', 'yticklabels', label, 'TickDir', 'out')
            xlabel('Time (s)');
            ylabel('Seizure');
            axis tight;
            ax = axis;
            plot([0, 0], [ax(3), ax(4)],'k:');            
            set(gca, 'YDir','reverse')
%             xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            xlim([-0.5, 10]);
            
            h = subplot(2, 1, 2);
            
            cfgtemp  = [];
            cfgtemp.channel         = ichan;
            cfgtemp.baseline        = [-0.5 0.2];
            cfgtemp.baselinetype    = 'relchange';
            cfgtemp.zlim            = 'maxabs';
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.interactive     = 'no';
            cfgtemp.figure          = h;
            cfgtemp.title           = '';
            ft_singleplotTFR(cfgtemp, TFR{ipart}.(markername));
            title('');
            xlabel('Time (s)'); ylabel('Frequency');
            c = colorbar;
            set(c, 'Location', 'southoutside', 'color', [0 0 0]);
            c.Title.String = 'Relative change in power';
            pos = get(c, 'Position');
            pos(2) = 0.03;
            set(c, 'pos', pos);
            colormap(jet(1000));            
            xlim([-0.5, 10]);
            
            set(findall(fig, '-property', 'fontsize'), 'fontsize', 7);

            % print to file
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'part-', num2str(ipart), '-pattern-', markername, '-channel-', LFP{ipart}.(markername).label{ichan}, '-timecourses.pdf')),'-r300');
            print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'part-', num2str(ipart), '-pattern-', markername, '-channel-', LFP{ipart}.(markername).label{ichan}, '-timecourses.png')),'-r300');
            
            close all
        end % ichan
    end % markername
end % ipart
