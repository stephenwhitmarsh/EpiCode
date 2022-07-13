function plot_PSTH_stats(cfg, psth_data)

nrows = 1;

for ipart = 1 : size(psth_data, 2)
    for markername = string(fields(psth_data{ipart}.psth))'
        nrows = max([nrows, size(psth_data{ipart}.stat.(markername), 2)]);
    end
end

for ipart = 1 : size(psth_data, 2)
    
    ncols = numel(fieldnames(psth_data{ipart}.psth));
    
    fig = figure;
    sgtitle(sprintf('%s : Spike density %s, %d trials (Hz)', cfg.prefix(1:end-1), markername, size(psth_data{ipart}.psth.(markername).trialinfo, 1)),...
        'fontsize', 15, 'fontweight', 'bold', 'interpreter', 'none');
    orient(fig, 'portrait');
    col = 1;
    for markername = string(fields(psth_data{ipart}.psth))'
        disp(size(psth_data{ipart}.psth.(markername).label, 2));
        
        row = 1;
        for itemp = 1 : size(psth_data{ipart}.psth.(markername).label, 2)
            
            subaxis(nrows, ncols, (row-1) * ncols + col, 'SpacingVert', 0.01);
            hold on
            
            if size(psth_data{ipart}.psth.(markername).label, 2) == 1
                bar(psth_data{ipart}.psth.(markername).time, psth_data{ipart}.psth.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            else
                bar(psth_data{ipart}.psth.(markername).time, psth_data{ipart}.psth.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            end
            
            % plot positive clusters
            lag = size(psth_data{ipart}.psth.(markername).time, 2) - size(psth_data{ipart}.stat.(markername){itemp}.time, 2);
            
            if isfield(psth_data{ipart}.stat.(markername){itemp}, 'posclusters')
                for ipos = 1 : size(psth_data{ipart}.stat.(markername){itemp}.posclusters, 2)
                    if psth_data{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        sel = find(psth_data{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
                        if length(sel) == 1
                            if sel == 1
                                continue
                            end
                            w = psth_data{ipart}.stat.(markername){itemp}.time(sel) - psth_data{ipart}.stat.(markername){itemp}.time(sel-1);
                        else
                            w = 1;
                        end
                        bar(psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), w, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                    end
                end
            end
            
            % plot negative clusters
            if isfield(psth_data{ipart}.stat.(markername){itemp}, 'negclusters')
                for ineg = 1 : size(psth_data{ipart}.stat.(markername){itemp}.negclusters, 2)
                    if psth_data{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        sel = find(psth_data{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                        if length(sel) == 1
                            if sel == 1 
                                continue
                            end
                            w = psth_data{ipart}.stat.(markername){itemp}.time(sel) - psth_data{ipart}.stat.(markername){itemp}.time(sel-1);
                        else
                            w = 1;
                        end
                        bar(psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), w, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                    end
                end
            end
            axis tight
            %
            %
            %             if size(psth_data{ipart}.psth.(markername).label, 2) == 1
            %                 bar(psth_data{ipart}.psth.(markername).time, psth_data{ipart}.psth.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            %             else
            %                 bar(psth_data{ipart}.psth.(markername).time, psth_data{ipart}.psth.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            %             end
            %             if isfield(psth_data{ipart}.stat.(markername){itemp}, 'posclusters')
            %                 for ipos = 1 : size(psth_data{ipart}.stat.(markername){itemp}.posclusters, 2)
            %                     if psth_data{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
            %                         sel = find(psth_data{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
            %                         if size(psth_data{ipart}.psth.(markername).label, 2) == 1
            %                             lag = size(psth_data{ipart}.psth.(markername).avg, 1) - size(psth_data{ipart}.stat.(markername){itemp}.mask, 2);
            %                             if length(sel) == 1
            %                                 bar([psth_data{ipart}.stat.(markername){itemp}.time(sel)-0.001, psth_data{ipart}.stat.(markername){itemp}.time(sel)+0.01], [psth_data{ipart}.psth.(markername).avg( sel+lag), psth_data{ipart}.psth.(markername).avg(sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
            %                             else
            %                                 bar( psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
            %                             end
            %                         else
            %                             lag = size(psth_data{ipart}.psth.(markername).avg, 2) - size(psth_data{ipart}.stat.(markername){itemp}.mask, 2);
            %                             if length(sel) == 1
            %                                 bar([psth_data{ipart}.stat.(markername){itemp}.time(sel)-0.001, psth_data{ipart}.stat.(markername){itemp}.time(sel)+0.01], [psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
            %                             else
            %                                 bar( psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
            %                             end
            %                         end
            %                     end
            %                 end
            %             end
            %
            %             if isfield(psth_data{ipart}.stat.(markername){itemp}, 'negclusters')
            %                 for ineg = 1 : size(psth_data{ipart}.stat.(markername){itemp}.negclusters, 2)
            %                     if psth_data{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
            %                         sel = find(psth_data{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
            %                         if size(psth_data{ipart}.psth.(markername).label, 2) == 1
            %                             lag = size(psth_data{ipart}.psth.(markername).avg, 1) - size(psth_data{ipart}.stat.(markername){itemp}.mask, 2);
            %                             if length(sel) == 1
            %                                 bar([psth_data{ipart}.stat.(markername){itemp}.time(sel)-0.001, psth_data{ipart}.stat.(markername){itemp}.time(sel)+0.01], [psth_data{ipart}.psth.(markername).avg( sel+lag), psth_data{ipart}.psth.(markername).avg(sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
            %                             else
            %                                 bar( psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
            %                             end
            %                         else
            %                             lag = size(psth_data{ipart}.psth.(markername).avg, 2) - size(psth_data{ipart}.stat.(markername){itemp}.mask, 2);
            %                             if length(sel) == 1
            %                                 bar([psth_data{ipart}.stat.(markername){itemp}.time(sel)-0.001, psth_data{ipart}.stat.(markername){itemp}.time(sel)+0.01], [psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
            %                             else
            %                                 bar( psth_data{ipart}.stat.(markername){itemp}.time(sel), psth_data{ipart}.psth.(markername).avg(itemp, sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
            %                             end
            %                         end
            %                     end
            %                 end
            %             end
            
            % patch baseline
            x = cfg.stats.bl.(markername);
            y = ylim;
            patch('XData', [x(1) x(2), x(2), x(1)], 'YData', [y(1), y(1), y(2), y(2)], 'edgecolor', 'none', 'facecolor', 'k', 'facealpha', 0.1);
            
            box off
            set(gca, 'TickDir', 'out');
            %xlim(cfg.epoch.toi.(markername));
            ytick = get(gca, 'ytick');
            set(gca, 'ytick', ytick(end));
            y = split(psth_data{ipart}.psth.(markername).label{itemp}, '_');
            ylabel(y{2}, 'interpreter', 'none');
            if row ~= size(psth_data{ipart}.psth.(markername).label, 2)
                set(gca,'xtick',[]);
                set(gca,'xticklabel',[]);
                set(gca,'XColor','none')
            else
                xlabel('Time (seconds)');
            end
            set(gca,'FontSize', 7, 'tickdir', 'out');
%             if row == 1
%                 title(sprintf('%s (%d trials)', markername, size(psth_data{ipart}.psth.(markername).trialinfo, 1)), 'fontsize', 15);
%             end
            row = row + 1;
        end
        col = col + 1;
    end
    
    % print to file
    %     fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    %     set(fig,'PaperOrientation','portrait');
    %     set(fig,'PaperUnits','normalized');
    %     set(fig,'PaperPosition', [0 0 1 1]);
    fname_fig = fullfile(cfg.imagesavedir, 'psth', strcat(cfg.prefix, 'p', num2str(ipart), '_stats_bargraphs'));
    %     isdir_or_mkdir(fileparts(fname_fig));
    savefigure_own(fig, fname_fig, 'png', 'pdf', 'close');
    %     exportgraphics(fig, [fname_fig, '.jpg'], 'resolution', 150);
    %     exportgraphics(fig, [fname_fig, '.pdf']);
    %     close(fig);
    
end


