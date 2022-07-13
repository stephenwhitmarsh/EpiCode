function plot_spike_waveforms(cfg, markerlist, spikewaveformstats, spikestats, spikewaveforms)

% Plot spike waveform data. If the argument spikewaveforms is omitted, then
% it only plots average data.
%
% Use as :
%       plot_spike_waveforms(cfg, spikewaveformstats, spikestats, spikewaveforms)
%
% INPUT :
% - cfg.prefix               : prefix used in the name of the image
% - cfg.imagesavedir         : where to save images
% - cfg.spike.RPV            : time for refractory period violation
% - cfg.plotspike.plotraw    : true/false, whether to plot raw spike waveforms
%                             (default = true)
% cfg.plotspike.max_nr_of_spikes : maximum number of raw spikes plotted. If
%                              there are more spikes, then a random subselection 
%                              of spikes is done (default = 1000) 
% - cfg.plotspike.plotavg    : true/false, whether to plot the mean of the
%                              spike waveforms (default = true)
% - cfg.plotspike.plotstd    : true/false, whether to plot the std of the
%                              spike waveforms (default = true)
% - cfg.plotspike.isi_lim    : x limits of the isi plot (default = [0 0.05]
% - cfg.plotspike.adapt_y    : 'yes', 'no', wheter to adapt y scaling in
%                               case of outlier raw spike traces (dafault = yes).
% - cfg.plotspike.suffix     : suffix used in the name of the image, usefull
%                              if you need to apply this function with different
%                              parameters (default = [])
% - cfg.plotspike.img_format : format used to save the plot. See :
%                              savefigure_own.m (default = "png"). Can be a 
%                              list of formats, ie : ["png", "pdf"];
% - markerlist               : name of the period to plot, as defined in
%                              the setparams script. It should be a string
%                              (ie : "marker1"), or a list of strings (ie :
%                              ["marker1", "marker2"])
% - spikewaveformstats       : output from spikeWaveformStats.m
% - spikestats               : output from readSpikeTrials
% - spikewaveforms (optional): output from readSpikeWaveforms
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.
% 
cfg.plotspike                  = ft_getopt(cfg, 'plotspike', []);
cfg.plotspike.plotraw          = ft_getopt(cfg.plotspike, 'plotraw', true);
cfg.plotspike.max_nr_of_spikes = ft_getopt(cfg.plotspike, 'max_nr_of_spikes', 1000);
cfg.plotspike.plotavg          = ft_getopt(cfg.plotspike, 'plotavg', true);
cfg.plotspike.plotstd          = ft_getopt(cfg.plotspike, 'plotstd', true);
cfg.plotspike.isi_lim          = ft_getopt(cfg.plotspike, 'isi_lim', [0 0.05]);
cfg.plotspike.adapt_y          = ft_getopt(cfg.plotspike, 'adapt_y', 'yes');
cfg.plotspike.suffix           = ft_getopt(cfg.plotspike, 'suffix', []);
cfg.plotspike.img_format       = ft_getopt(cfg.plotspike, 'img_format', "png");

if nargin == 4
    cfg.plotspike.plotraw = false;
elseif nargin < 4
    error('not enough arguments : the first 4 arguments are mandatory');
end

cfg.plotspike.img_format = cellstr(cfg.plotspike.img_format);

for ipart = 1:size(spikewaveformstats,2)
    
    for markername = markerlist
        group = spikewaveformstats{ipart}.cluster_group;
        
        for itype = ["good", "mua"]
            
            n_units = sum(contains(group, itype));
            if n_units == 0
                continue
            end
            n_subplots = ceil(n_units/5);
            
            %one plot per type
            iplot = -1;
            fig = figure;
            sgtitle(sprintf('%s : %s', cfg.prefix(1:end-1), itype), 'interpreter', 'none');
            
            ft_progress('init', 'text');%, sprintf('%s p%d %s : %s %d/%d', cfg.prefix(1:end-1), ipart, markername, itype
            for i_unit = 1:length(group)
                ft_progress(0, '%s p%d : %s %d/%d', cfg.prefix(1:end-1), ipart, itype, i_unit, length(group));

                if ~contains(group{i_unit}, itype)
                    continue
                end
                
                iplot = iplot+2;
                rpv   = spikestats{ipart}.(markername){i_unit}.RPV * 100;
                
                subplot(ceil(n_units/n_subplots),2*n_subplots,iplot);hold on;
                
                %plot ISI
                bar(spikestats{ipart}.(markername){i_unit}.isi_avg_time * 1000, spikestats{ipart}.(markername){i_unit}.isi_avg, 1, 'facecolor', [0 0 0], 'edgecolor', 'none');
                idx = spikestats{ipart}.(markername){i_unit}.isi_avg_time <= cfg.spike.RPV;
                bar(spikestats{ipart}.(markername){i_unit}.isi_avg_time(idx) * 1000, spikestats{ipart}.(markername){i_unit}.isi_avg(idx), 1, 'facecolor', [1 0 0], 'edgecolor', 'none');
                set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out', 'FontSize', 5);
                axis tight
                xlim(cfg.plotspike.isi_lim * 1000);
                %xlabel('Time (ms)'); ylabel('Spikecount');
                
                ax = axis;
                titlepos = title(sprintf('%s, RPV=%.2f%%, %d spikes',spikestats{ipart}.(markername){i_unit}.label, rpv, length(spikestats{ipart}.(markername){i_unit}.isi)),'Interpreter', 'none', 'FontSize', 10);
                titlepos.Position(1) = ax(1);
                titlepos.HorizontalAlignment = 'left';
                
                subplot(ceil(n_units/n_subplots),2*n_subplots, iplot+1); hold on;
                if cfg.plotspike.plotraw
                    color_avg = 'r';
                    linewidth_avg = 1;
                    alpha_std = 0.4;
                else
                    color_avg = 'k';
                    linewidth_avg = 1;
                    alpha_std = 0.2;
                end
                
                if isfield(spikewaveformstats{ipart}.waveformavg{i_unit}, 'time')
                    flip  = -spikewaveformstats{ipart}.peak_direction(i_unit);
                    if istrue(cfg.plotspike.plotraw)
                        %plot only randomly selected trials if there
                        %are more than cfg.plotspike.max_nr_of_spikes trials
                        if size(spikewaveforms{ipart}{i_unit}.trial, 2) > cfg.plotspike.max_nr_of_spikes
                            trial_list = randperm(size(spikewaveforms{ipart}{i_unit}.trial,2),cfg.plotspike.max_nr_of_spikes);
                        else
                            trial_list = 1:length(spikewaveforms{ipart}{i_unit}.trial);
                        end
                        for itrial = trial_list
                            p = plot(spikewaveforms{ipart}{i_unit}.time{itrial}, spikewaveforms{ipart}{i_unit}.trial{itrial}*flip, 'k');
                            p.Color(4) = 0.2;
                        end
                    end
                    if istrue(cfg.plotspike.plotstd)
                        x = spikewaveformstats{ipart}.waveformavg{i_unit}.time;
                        y = spikewaveformstats{ipart}.waveformavg{i_unit}.avg.*flip;
                        ystd = sqrt(spikewaveformstats{ipart}.waveformavg{i_unit}.var).*flip;
                        p = patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], color_avg, 'edgecolor', 'none', 'facealpha', alpha_std);
                        %p.ZData = ones(size(p.YData)) .* -1;
                    end
                    if istrue(cfg.plotspike.plotavg)
                        plot(spikewaveformstats{ipart}.waveformavg{i_unit}.time, spikewaveformstats{ipart}.waveformavg{i_unit}.avg.*flip, 'color', color_avg, 'linewidth', linewidth_avg);
                        plot(spikewaveformstats{ipart}.halfwidth.x(i_unit,:), spikewaveformstats{ipart}.halfwidth.y(i_unit,:).*flip, '-x', 'color', color_avg);
                        plot(spikewaveformstats{ipart}.troughpeak.x(i_unit,:), [0, 0], '-x', 'color', color_avg);
                        axis tight
                    end
%                     if istrue(cfg.plotspike.plotraw)
%                         %avoid aberrant scaling due to outlier spikes
                        ylim([min(spikewaveformstats{ipart}.waveformavg{i_unit}.avg*flip)*2 max(spikewaveformstats{ipart}.waveformavg{i_unit}.avg*flip)*10]);
%                     end
                end
                xticklabels(xticks.*1000);
                %xlabel('Time (ms)');
                %ylabel('uV');
                set(gca, 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out', 'FontSize', 5);
                axis tight 
                %scale y limits if needed
                if istrue(cfg.plotspike.adapt_y)
                    yinf = min(spikewaveformstats{ipart}.waveformavg{i_unit}.avg*flip)*2;
                    ysup = max(spikewaveformstats{ipart}.waveformavg{i_unit}.avg*flip)*10;
                    y = ylim;
                    if y(1) < yinf * 3
                        ylim([yinf, y(2)]);
                    end
                    y = ylim;
                    if y(2) > ysup * 3
                        ylim([y(1), ysup]);
                    end
                end
                
            end %iunit
            ft_progress('close');
            
            %print figure to file
            figname = fullfile(cfg.imagesavedir, 'plot_spike_waveforms', sprintf('%sp%d_%s_spikewaveforms_%s%s', cfg.prefix, ipart, markername, itype, cfg.plotspike.suffix));
            savefigure_own(fig, figname, 'portrait', cfg.plotspike.img_format{:}, 'close');
            
        end %itype
    end %markernmae
end %markerlist