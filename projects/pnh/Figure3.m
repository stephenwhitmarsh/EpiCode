function Figure3(cfg_fig)

cfg = pnh_setparams;

cfg_fig.plot.name = {'PSW', 'FA', 'ES'};

ipart = 1;

for ipatient = 1 : 4
    SpikeTrials{ipatient}           = readSpikeTrials_MuseMarkers(cfg{ipatient});
    SpikeDensity{ipatient}          = spikeTrialDensity(cfg{ipatient});
    SpikeWaveforms{ipatient}        = readSpikeWaveforms(cfg{ipatient});
    SpikeStats_windowed{ipatient}   = spikeTrialStats(cfg{ipatient}); 
    SpikeTrials_windowed{ipatient}  = readSpikeTrials_windowed_samples(cfg{ipatient});
end

% so that FA of nodule 3 will end up in line with PSWs
SpikeTrials{3}{1}.PSW           = SpikeTrials{3}{1}.FA; 
SpikeTrials{3}{1}               = rmfield(SpikeTrials{3}{1}, 'FA');
SpikeDensity{3}{1}.stat.PSW     = SpikeDensity{3}{1}.stat.FA;
SpikeDensity{3}{1}.stat         = rmfield(SpikeDensity{3}{1}.stat, 'FA');
SpikeDensity{3}{1}.sdf_bar.PSW  = SpikeDensity{3}{1}.sdf_bar.FA;
SpikeDensity{3}{1}.sdf_bar      = rmfield(SpikeDensity{3}{1}.sdf_bar, 'FA');
SpikeDensity{3}{1}.psth.PSW     = SpikeDensity{3}{1}.psth.FA;
SpikeDensity{3}{1}.psth         = rmfield(SpikeDensity{3}{1}.psth, 'FA');

% which units to plot
units(1, :) = [1, 8];
units(2, :) = [3, 9];
units(3, :) = [1, 9];
units(4, :) = [1, 4];

% configure figure
papersize = 1000;
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [20 -60  papersize papersize*sqrt(2)]);
set(fig, 'Renderer', 'Painters');

ntemp       = 2;
nmarker     = 3;
npatient    = 4;

w           = 1/3;
hratio      = 1/3.00;
vratio      = 1/16.8;
htop        = 1/16 * 0.8;
hbottom     = 1/16 * 0.8;
spacehor    = 0.1;
rshift      = 0.015;
upshift     = -0.01;

w = w * (1-spacehor);

% first print average LFP and TFR
for ipatient = 1 : 4
    
    imarker = 1;
    for markername = string(cfg_fig.plot.name)
        
        if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername)
            imarker = imarker + 1;
            continue
        end
    
        for tempnr = 1 : ntemp
            
            row     = (tempnr-1) * 2 + (ipatient - 1) * 4 + 1;
            itemp   = units(ipatient, tempnr);
            s1      = axes('Position', [hratio*(imarker-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
            set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');
            
            hold on
            if size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2) == 1
                bar(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time, SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            else
                bar(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time, SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
            end
            
            if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'posclusters')
                for ipos = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters, 2)
                    if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg{ipatient}.stats.alpha
                        sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
                        if size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2) == 1
                            lag = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 1) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                            if length(sel) == 1
                                bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg( sel+lag), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                            else
                                bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                            end
                        else
                            lag = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                            if length(sel) == 1
                                bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                            else
                                bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                            end
                        end
                    end
                end
            end
            
            if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'negclusters')
                for ineg = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters, 2)
                    if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg{ipatient}.stats.alpha
                        sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                        if size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2) == 1
                            lag = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 1) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                            if length(sel) == 1
                                bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg( sel+lag), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                            else
                                bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                            end
                        else
                            lag = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                            if length(sel) == 1
                                bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                            else
                                bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                            end
                        end
                    end
                end
            end
            axis tight;
            xlim(cfg{ipatient}.epoch.toi.(markername));
            y = ylim;
            yticks(floor(y(end)));
            
            if tempnr == 1
                switch markername
                    case "ES"
                        th = title("Epileptic Spike", 'FontSize', 14);
                        
                    case "PSW"
                        th = title("Periodic Slow Waves", 'FontSize', 14);
                        if ipatient == 3
                            th = title("Periodic Fast Activity", 'FontSize', 14);
                        end
                    case "FA"
                        th = title("Fast Activity", 'FontSize', 14);
                end
                set(th,'position',get(th,'position')-[0 0.07 0])     
            end
            
            if ~((imarker == 1) || (ipatient == 3 && imarker == 2))
                ylabel([]);
            else
                l = ylabel('Count');
                set(l, 'Units', 'normalized');
                x = l.Position;
                x(1) = -0.025;
                set(l, 'Position', x);
            end
            
            % baseline line
            y = ylim;
            y = y(1) + (y(2) - y(1)) * 0.05;
            x = xlim;
            width = x(2) - x(1);
            line(cfg{ipatient}.stats.bl.(markername), [y, y], 'linewidth', 2, 'color', 'k');
            text(cfg{ipatient}.stats.bl.(markername)(1) + width * 0.01, y, 'Baseline', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            
            % time line
            if strcmp(markername, "ES")
                duration = 0.5;
                t = sprintf('%.1f second', duration);
            else
                duration = 1;
                t = sprintf('1 second');
            end
            line([cfg{ipatient}.spike.toi.(markername)(end) - duration, cfg{ipatient}.spike.toi.(markername)(end)], [y, y], 'linewidth', 2, 'color', 'k');
            text(cfg{ipatient}.spike.toi.(markername)(end) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
            y = ylim;
            
            % label MUA/SUA
            if imarker == 1 || (imarker == 2 && ipatient == 3)
                if strfind(SpikeTrials_windowed{ipatient}{ipart}.window.cluster_group{itemp}, 'good')
                    text(cfg{ipatient}.stats.bl.(markername)(1) + width * 0.01, y(2) * 0.95, "SUA", 'color', 'k', 'fontsize', 8);
                else
                    text(cfg{ipatient}.stats.bl.(markername)(1) + width * 0.01, y(2) * 0.95, "MUA", 'color', 'k', 'fontsize', 8); %0.9
                end
            end
            
            % add inset with waveshape
            if imarker == 1 || (imarker == 2 && ipatient == 3)
                
                y_all   = vertcat(SpikeWaveforms{ipatient}{ipart}{itemp}.trial{:})';
                r       = rms(y_all - mean(y_all, 2));
                [~, i]  = sort(r, 'ascend'); 
                y_sel   = y_all(:, i(1:100));
                
                x = s1.Position;
                x(1) = x(1) + 0.01;
                x(2) = x(2) + 0.020; % 0.015
                x(3) = x(3) / 6;
                x(4) = x(4) / 1.5;   
                inset1 = axes('position', x, 'color', 'none', 'XColor','none', 'YColor', 'none'); hold on; %left bottom width height
 
                for i = 1 : size(y_sel, 2)
                    lh = plot(y_sel(:, i), 'k'); lh.Color = [0, 0, 0, 0.2];
                end
                plot(mean(y_all, 2), 'k', 'linewidth', 2);
                plot(mean(y_all, 2), 'w', 'linewidth', 1);    
            end
      
%           % add inset with ISI histogram
%             if imarker == 1 || (imarker == 2 && ipatient == 3)
%                 x = s1.Position;
%                 x(1) = x(1) + 0.07;
%                 x(2) = x(2) + 0.017;
%                 x(3) = x(3) / 6;
%                 x(4) = x(4) / 1.5;
%                 inset2 = axes('position', x, 'color', 'none', 'YColor','none', 'Tickdir', 'out'); hold on;
%                 bar(SpikeStats_windowed{ipatient}{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats_windowed{ipatient}{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0, 0.5, 0], 'edgecolor', 'none');
%                 x = find(SpikeStats_windowed{ipatient}{ipart}.window{itemp}.isi_avg_time < 0.002, 1, 'last');                
%                 bar(SpikeStats_windowed{ipatient}{ipart}.window{itemp}.isi_avg_time(1:x) * 1000, SpikeStats_windowed{ipatient}{ipart}.window{itemp}.isi_avg(1:x), 1, 'facecolor', [0.5, 0, 0], 'edgecolor', 'none');
%                 xlim([0, 30]);
%                 xticks([0, 30]);
%                 l = xlabel('time');
%                 x = get(l,'Position');
%                 set(l, 'Position', x .* [1, 0.1, 1]);
%                 y = ylim;
%             end

            % next row
            row = (tempnr-1) * 2 + (ipatient - 1) * 4 + 1 + 1;
            s2 = axes('Position', [hratio*(imarker-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + (vratio-hbottom) + upshift, w, hbottom], 'Color', 'none', 'XColor', 'none');
            
            % plot raster
            cfg_raster              = [];
            cfg_raster.trialborders = 'no';
            cfg_raster.spikechannel = itemp;
            ft_spike_plot_raster(cfg_raster, SpikeTrials{ipatient}{ipart}.(markername));
            set(gca,'Ydir','reverse')
            
            % increased units in raster for visibility
            if strcmp(markername, "ES") && ~(ipatient == 2 && itemp == 9) && ~(ipatient == 4 && itemp == 1)
                lines = findobj(s2,'Type','Line');
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
            end
            yticks(size(SpikeDensity{ipatient}{ipart}.psth.(markername).trial, 1));
            
            if row ~= 16
                set(s2, 'box', 'off', 'xticklabel', [], 'TickDir', 'out', 'xlabel', []);
                xticks([]);
            end
            
            if ~((imarker == 1)) % || (ipatient == 3 && imarker == 2))
                ylabel([]);
            else
                l = ylabel('Trial');
                set(l, 'Units', 'normalized');
                x = l.Position;
                x(1) = -0.025;
                set(l, 'Position', x);
            end
            
        end
        imarker = imarker + 1;
    end
end

fname = fullfile(cfg{ipart}.imagesavedir, 'article', 'Figure3');
% exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));
close all
disp('done'); 