function plotstats(cfg)

nrows = 18;
ncols = 3;

%     SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg);
SpikeDensity_timelocked = spikeTrialDensity(cfg);

for ipart = 1 : size(SpikeDensity_timelocked, 2)
    
    fig = figure;
    orient(fig, 'portrait');
    col = 1;
    for markername = string(fields(SpikeDensity_timelocked{ipart}.sdf_bar))'
        
        row = 1;
        for itemp = 1 : size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).label, 2)
            
            subaxis(nrows, ncols, (row-1) * ncols + col, 'SpacingVert', 0.01);
            
            bar(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).time, SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', [127/255,127/255,127/255]);
            hold;
            if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'posclusters')
                for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        lag = size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.mask, 2);
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
                        bar(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', [252/255,187/255,62/255]);
                    end
                end
            end
            
            subaxis(nrows, ncols, (row-1) * ncols + col, 'SpacingVert', 0.01);
            if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'negclusters')
                for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters, 2)
                    if SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        lag = size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.mask, 2);
                        sel = find(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                        bar(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', [70/255,93/255,250/255]);
                    end
                end
            end
            
            box off
            set(gca, 'TickDir', 'out');
            
            if row ~= nrows
                set(gca,'xtick',[]);
                set(gca,'xticklabel',[]);
                set(gca,'XColor','none')
            else
                xlabel('Time (seconds)');
                ylabel('Count');
            end
            set(gca,'FontSize',6);
            if row == 1
                title(markername);
            end
            row = row + 1;
            
        end
        col = col + 1;
    end
    
    % print ISI to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','portrait');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    fname_fig = fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'stats_bargraphs'));
    exportgraphics(fig, [fname_fig, '.jpg'], 'resolution', 600);
    exportgraphics(fig, [fname_fig, '.pdf']);

end


