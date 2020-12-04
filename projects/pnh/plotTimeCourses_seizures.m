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
        
        chanlist = unique( LFP{ipart}.(markername).trialinfo.chan1( ~ismissing(LFP{ipart}.(markername).trialinfo.chan1) ))';
        chanlist(strcmp(chanlist, "diffuse")) = [];
        
        % get max range in dataset
        maxrange = 0;
        maxtrials = 0;
        for channame = chanlist
            trials  = find(strcmp(LFP{ipart}.(markername).trialinfo.chan1, channame));
            maxtrials = max([maxtrials, length(trials)]);
            channel = strcmp(string(LFP{ipart}.(markername).label), strcat('_',channame));
            for itrial = trials'
                y = max(abs(LFP{ipart}.(markername).trial{itrial}(channel,:)));
                maxrange = max(abs([maxrange, y]));
            end
        end
        maxrange = maxrange / 2;
        ymax = maxrange * (maxtrials-1) + maxrange/2;
        ymin = -maxrange/2;
        for channame = chanlist
            
            cfgtemp         = [];
            cfgtemp.trials  = strcmp(LFP{ipart}.(markername).trialinfo.chan1, channame);
            cfgtemp.channel = find(strcmp(string(LFP{ipart}.(markername).label), strcat('_',channame)));
            LFP_sel         = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
            TFR_sel         = ft_selectdata(cfgtemp, TFR{ipart}.(markername));
            
            % plot
            
            fig = figure;
            subplot(2,1,1); hold;
            
            n = 0;
            label = [];
            ytick = [];
            
            for itrial = 1 : size(LFP_sel.trial, 2)
                title(strcat('Channel:'," ", channame),'interpreter', 'none');
                ytick = [ytick, n*maxrange];
                x = LFP_sel.time{itrial};
                y = LFP_sel.trial{itrial};
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
            %             xlim([cfg.epoch.toi.(markername)(1), cfg.epoch.toi.(markername)(2)]);
            xlim([-0.5, 10]);
            ylim([-maxrange/2, ymax]);
            
            top     = max(LFP_sel.trial{end}) + (n-1)*maxrange;
            bottom  = min(LFP_sel.trial{end});
            middle  = (top-bottom)/2;
            ymid    = (ymax-ymin) /2;
            
            ax = axis;
%             plot([ax(1), ax(2)], [middle, middle], 'r');
%             plot([ax(1), ax(2)], [top, top] , 'r');
%             plot([ax(1), ax(2)], [bottom, bottom], 'r');
            
%             plot([ax(1), ax(2)], [ymax, ymax], 'b');
%             plot([ax(1), ax(2)], [ymin, ymin], 'b');
%             plot([ax(1), ax(2)], [ymid, ymid], 'b');
  
            shift = ymid - middle;            
            axis([ax(1), ax(2), ax(3)-shift, ax(4)-shift]);

            ax = axis;
            plot([0, 0], [ax(3), ax(4)],'k:');

            
            h = subplot(2, 1, 2);
            
            cfgtemp  = [];
            cfgtemp.baseline        = [-0.5 0.2];
            cfgtemp.baselinetype    = 'relchange';
            cfgtemp.zlim            = 'maxabs';
            cfgtemp.title           = '';
            cfgtemp.colorbar        = 'no';
            cfgtemp.interactive     = 'no';
            cfgtemp.figure          = h;
            cfgtemp.title           = '';
            ft_singleplotTFR(cfgtemp, TFR_sel);
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
            print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'channel_' ,channame)),'-r300');
            print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'channel_' ,channame)),'-r300');
            
            close all
        end % ichan
    end % markername
end % ipart
