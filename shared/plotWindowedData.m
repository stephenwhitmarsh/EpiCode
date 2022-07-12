function plotWindowedData(cfg, MuseStruct, markername, varargin)

% PLOTWINDOWEFDATA plots windowed data from different data types a common 
% axis.
%
% use as
%   plotWindowedData(cfg, data1, data2, data3);
%
% required fields:
%     data can be any of the following in any order: 
%        SpikeTrials
%        SpikeStats
%        FFT
%        hypnogram
%
%     cfg                 = [];
%     cfg.ipart           = ipart;
%     cfg.filename        = fullfile(config{ipatient}.imagesavedir, sprintf('windowed_s%sp%d.jpg', config{ipatient}.prefix, ipart));
%     cfg.xlim            = seconds([MuseStruct{1}{1}.starttime, MuseStruct{1}{2}.endtime] - MuseStruct{1}{2}.markers.CriseStart.clock);
%     
% optional fields:
%     cfg.orientation     = 'landscape' or 'portrait'
%     cfg.offset          = e.g.: MuseStruct{1}{2}.markers.CriseStart.clock; 
%     
% followed by the following fields, one per row, e.g.
%     cfg.type{1}         = 'power'; % or: 'relpower', 'spike', 'trialinfo' 
%     cfg.frequency{1}    = [1, 7];
%     cfg.channel{1}      = 'all';
%     cfg.title{1}        = sprintf('Power %d-%dHz', cfg.frequency{1}(1), cfg.frequency{1}(2));
%     cfg.plotart(1)      = false; % plot artefacts as shaded areas
%     cfg.log(1)          = true;  % yaxis log
%     cfg.hideart(1)      = false; % hide data when artefacted
%     cfg.marker_indx{1}  = FFT{ipart}.(markername).trialinfo.BAD_cnt>0; % plot
%     markers on the same axis as the data
%     cfg.marker_sign{1}  = '.r'; % how to plot markers
%     cfg.marker_label{1} = 'artefact'; % their label
%     
%     cfg.type{2}         = 'power';
%     cfg.frequency{2}    = [8, 14];
%     cfg.channel{2}      = 'all';
%     cfg.title{2}        = sprintf('Power %d-%dHz', cfg.frequency{2}(1), cfg.frequency{2}(2));
%     cfg.plotart(2)      = false;
%     cfg.log(2)          = true;
%     cfg.hideart(2)      = true;
%     
%     cfg.type{3}         = 'relpower';
%     cfg.frequency1{3}   = [1, 7];
%     cfg.frequency2{3}   = [8, 14];
%     cfg.channel{3}      = 'all';
%     cfg.title{3}        = sprintf('Power (%d-%dHz)/(%d-%dHz)', cfg.frequency1{3}(1), cfg.frequency1{3}(2), cfg.frequency2{3}(1), cfg.frequency2{3}(2));
%     cfg.plotart(3)      = true;
%     cfg.log(3)          = false;
%     
   
% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
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

ncols           = 1;
nrows           = size(cfg.type, 2);
ipart           = ft_getopt(cfg, 'ipart', 1); % which part
cfg.plotart     = ft_getopt(cfg, 'plotart', false(1, nrows)); % plot artefacts as shaded regions
cfg.log         = ft_getopt(cfg, 'log', false(1, nrows)); % plot artefacts as shaded regions
cfg.unit        = ft_getopt(cfg, 'unit', ones(1, nrows));
cfg.offset      = ft_getopt(cfg, 'offset', []);
cfg.marker      = ft_getopt(cfg, 'marker', []);
cfg.orientation = ft_getopt(cfg, 'orientation', 'landscape');
cfg.minbadlength = ft_getopt(cfg, 'minbadlength', 0);

for i = 1 : nargin - 2
    try
        if isfield(varargin{i}{ipart}.(markername), 'sample')
            SpikeTrials = varargin{i};
            fprintf('Data input %d is of type "SpikeTrials"\n', i);
        end
    catch
    end
    
    try
        if isfield(varargin{i}{ipart}.(markername), 'powspctrm')
            FFT = varargin{i}{ipart}.(markername);
            fprintf('Data input %d is of type "Powerspctrm"\n', i);
        end
    catch
    end
    
    try
        if istable(varargin{i})
            hypnogram = varargin{i};
            fprintf('Data input %d is of type "Hypnogram"\n', i);
        end
    catch
    end
    
    try
        if isfield(varargin{i}{ipart}.(markername){1}, 'isi')
            SpikeStats = varargin{i};
            fprintf('Data input %d is of type "SpikeStats"\n', i);
        end
    catch
    end
    
end

i = 1;
for idir = 1 : size(MuseStruct{ipart}, 2)
    if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
        warning('No BAD marker field found in part %d, dir %d!\n', ipart, idir)
        continue
    end
    if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'clock')
        fprintf('No BAD markers found in part %d, dir %d (lucky you!)\n', ipart, idir)
        continue
    end    
    for ibad = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.clock, 2)
        if ~isempty(cfg.offset)
            bad_t1(i) = seconds(MuseStruct{ipart}{idir}.markers.BAD__START__.clock(ibad) - cfg.offset);
            bad_t2(i) = seconds(MuseStruct{ipart}{idir}.markers.BAD__END__.clock(ibad) - cfg.offset);
        else
            bad_t1(i) = MuseStruct{ipart}{idir}.markers.BAD__START__.clock(ibad);
            bad_t2(i) = MuseStruct{ipart}{idir}.markers.BAD__END__.clock(ibad);
        end
        i = i + 1;
    end
end

% if cfg.minbadlength > 0
%     bad_diff = bad_t2 - bad_t1;
%     sel = bad_diff < cfg.minbadlength;
%     bad_t1(sel) = [];
%     bad_t2(sel) = [];
% end

% prepare subplot dimensions
w           = 1/ncols - 0.12; % add some space for legend
hratio      = 1/ncols; 
vratio      = 1/nrows * 0.95; % to make space for x-labels
htop        = 1/nrows * 0.7; %*0.8
spacehor    = 0.1;
rshift      = 0.02; %18
upshift     = 0.0; %0.01
w           = w * (1-spacehor);
legend_x    = 0.88;

% hypnogram color scheme - with added gray for non-scored 
cm_hyp      = [cool(5); 0.8, 0.8, 0.8];

% create figure
papersize = 800;
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
if strcmp(cfg.orientation, 'landscape')
    set(fig, 'position', [20 -60 papersize*sqrt(2) papersize]);
else
    set(fig, 'position', [20 -60 papersize papersize*sqrt(2)]);
end
set(fig, 'Renderer', 'Painters');

for irow = 1 : nrows
    
    if isempty(cfg.type{irow})
        continue
    end
    
    icol    = 1;
    row     = irow;
    s1      = axes('Position', [hratio*(icol-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
    
    % some defaults
    try cfg.hideart(irow);  catch; cfg.hideart(irow)    = false; end
    try cfg.log(irow);      catch; cfg.log(irow)        = false; end
    try cfg.plotart(irow);  catch; cfg.plotart(irow)    = false; end
    try cfg.index{irow};    catch; cfg.index{irow}      = 1; end
    
    % only add x labels on bottom row
    if irow < nrows
        set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');
    else
        set(s1, 'XGrid', 'off', 'box', 'off', 'tickdir', 'out', 'Color', 'none', 'XMinorTick','on');
    end
    hold on
    
    title(cfg.title{irow});
    l = []; % legend

    switch cfg.type{irow}
        
        case 'hypnogram'
            
            % select part
            hypnogram = hypnogram(hypnogram.part == ipart, :);
            
            % apply offset
            if ~isempty(cfg.offset)
                hypnogram.starttime = seconds(hypnogram.starttime - cfg.offset);
                hypnogram.endtime   = seconds(hypnogram.endtime - cfg.offset);
            end
            
            % plot hypnogram
            plot_hyp_lines(hypnogram);
            axis tight;
            plot_hyp_colors(hypnogram, cm_hyp, ylim);
            fill([cfg.xlim(1), cfg.xlim(2), cfg.xlim(2), cfg.xlim(1)], [1 1 5 5], [1, 1, 1] ,'edgecolor', 'none');
            set(gca, 'children', flipud(get(gca, 'children'))); % flip order of images

            % legend
            clear p
            for ii = 1:size(cm_hyp,1)
                p(ii) = patch(NaN, NaN, cm_hyp(ii,:));
            end  
            hl              = legend(p, ["REM", "Wake", "Stage 1", "Stage 2", "Stage 3", "Unscored"], 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
            pos_legend      = get(hl, 'position');
            pos_legend(1)   = legend_x;
            set(hl, 'position', pos_legend);
            
        case 'power'
            
            % select and average data
            cfgtemp             = [];
            cfgtemp.frequency   = cfg.frequency{irow};
            cfgtemp.channel     = cfg.channel{irow};
            cfgtemp.avgoverfreq = 'yes';
            cfgtemp.avgoverchan = 'yes';
            data                = ft_selectdata(cfgtemp, FFT);
            time                = FFT.trialinfo.starttime + (FFT.trialinfo.endtime - FFT.trialinfo.starttime) / 2;
            
            % apply offset
            if ~isempty(cfg.offset)
                time = seconds(time - cfg.offset);
            end

            if cfg.hideart(irow)
                datatemp = data.powspctrm;
                datatemp(FFT.trialinfo.BAD_sec>cfg.minbadlength) = nan;
                bar(time, datatemp, 1, 'k');
            else
                bar(time, data.powspctrm, 1, 'k');
            end
                 
            % plot marker
            mplot = false;
            try
                plot(time(cfg.marker_indx{irow}), data.powspctrm(cfg.marker_indx{irow}), cfg.marker_sign{irow});
                mplot = true;
            catch
            end       
            
            ylim([0, std(data.powspctrm*4)]);

            % plot artefacts
            if cfg.plotart(irow)
                y = ylim;
                y2 = y;
                if istrue(cfg.log(irow)) && y2(1) == 0
                    h = findobj(gca,'Type','bar');
                    y2(1) = min(h.YData);
                end
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y2(1) y2(1) y2(2) y2(2)], 'r', 'edgecolor', 'r', 'facealpha', 0.4, 'edgealpha', 0.4);
                end
                ylim(y);
                hl = legend(cfg.channel{irow}, "artefact", 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
            else 
                hl = legend(cfg.channel{irow}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
            end
            
%             % legend
%             if mplot
%                 hl = legend({cfg.channel{irow}, cfg.marker_label{irow}}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
%             else
%                 hl = legend(cfg.channel{irow}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
%             end
            
            pos_legend = get(hl, 'position');
            pos_legend(1) = legend_x;
            set(hl, 'position', pos_legend);
            clear data

        case 'relpower'
            
            % select and average data
            cfgtemp             = [];
            cfgtemp.frequency   = cfg.frequency1{irow};
            cfgtemp.channel     = cfg.channel{irow};
            cfgtemp.avgoverfreq = 'yes';
            cfgtemp.avgoverchan = 'yes';
            data1               = ft_selectdata(cfgtemp, FFT);
            
            % select and average data
            cfgtemp             = [];
            cfgtemp.frequency   = cfg.frequency2{irow};
            cfgtemp.channel     = cfg.channel{irow};
            cfgtemp.avgoverfreq = 'yes';
            cfgtemp.avgoverchan = 'yes';
            data2               = ft_selectdata(cfgtemp, FFT); 
            
            % take ratio
            data                = data1;
            data.powspctrm      = data1.powspctrm ./ data2.powspctrm;
            clear data1 data2
            
            time                = FFT.trialinfo.starttime + (FFT.trialinfo.endtime - FFT.trialinfo.starttime) / 2;
            
            if ~isempty(cfg.offset)
                time = seconds(time - cfg.offset);
            end

            if cfg.hideart(irow)
                datatemp = data.powspctrm;
                datatemp(FFT.trialinfo.BAD_sec>cfg.minbadlength) = nan;
                bar(time, datatemp, 1, 'k');
            else
                bar(time, data.powspctrm, 1, 'k');
            end
            ylim([0, std(data.powspctrm*4)]);
            
            if cfg.plotart(irow)
                y = ylim;
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'r', 'edgecolor', 'r', 'facealpha', 0.4, 'edgealpha', 0.4);
                end
                ylim(y);
                hl = legend(cfg.channel{irow}, "artefact", 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');           
            else 
                hl = legend(cfg.channel{irow}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
            end
            pos_legend = get(hl, 'position');
            pos_legend(1) = legend_x;
            set(hl, 'position', pos_legend);
            clear data
            
        case 'spike'
            
            if isempty(cfg.unit{irow})
                if isfield(cfg, 'xlim')
                    xlim(cfg.xlim);
                end
                continue
            end
            
            % colormap
            clear l p
            cm  = lines(size(SpikeStats{ipart}.(markername), 2));
            for i = 1 : size(SpikeStats{ipart}.(markername), 2)
                l{i} = SpikeStats{ipart}.(markername){i}.label;
            end
            i = 0;
            if isempty(cfg.index{irow})
                cfg.index{irow} = 1;
            end
            for iunit = cfg.unit{irow}
                
                i = i + 1;
                time = SpikeStats{ipart}.(markername){iunit}.trialinfo.starttime + (SpikeStats{ipart}.(markername){iunit}.trialinfo.endtime - SpikeStats{ipart}.(markername){iunit}.trialinfo.starttime) / 2;
                
                % offset time
                if ~isempty(cfg.offset)
                    time = seconds(time - cfg.offset);
                end
                
                % in case an extra index is given, e.g. with Spiky distance

                                
                % plot data
                y = SpikeStats{ipart}.(markername){iunit}.(cfg.metric{irow})(cfg.index{irow}, :);
    
                if cfg.hideart(irow)
                    y(cfg.index{irow}, SpikeStats{ipart}.(markername){iunit}.trialinfo.BAD_sec>cfg.minbadlength) = nan;
                end
                
                if numel(cfg.index{irow}) == 1
                    %                     if numel(cfg.unit{irow}) == 1
                    %                         p(iunit) = bar(time, y, 1, 'k');
                    %                     else
                    p(iunit) = plot(time, y, '.-', 'color', cm(iunit, :));
                    %                     end
                    hl = legend(l{cfg.unit{irow}}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
                else
                    l2 = [];
                    for itarget = cfg.index{irow}
                        ci = strcmp(SpikeStats{ipart}.(markername){iunit}.dist_label{itarget}, SpikeTrials{ipart}.(markername).label);
                        p(iunit) = plot(time, y(itarget,:), '.-', 'color', cm(ci, :));
                        l2 = [l2; string(SpikeStats{ipart}.(markername){iunit}.dist_label{itarget})];
                    end
            
                    hl = legend(l2, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
                end
                
                % plot marker
                try
                    plot(time(cfg.marker_indx{irow}), SpikeStats{ipart}.(markername){iunit}.(cfg.metric{irow})(cfg.marker_indx{irow}), cfg.marker_sign{irow});
                catch
                end
                
            end % iunit
            
            pos_legend = get(hl, 'position');
            pos_legend(1) = legend_x;
            set(hl, 'position', pos_legend);
            
            if cfg.plotart(irow)
                y = [inf, -inf];
                for iunit = cfg.unit{irow}
                    for index = cfg.index{irow}
                        y(1) = min([SpikeStats{ipart}.(markername){iunit}.(cfg.metric{irow})(index, :), y(1)]);
                        y(2) = max([y(2), SpikeStats{ipart}.(markername){iunit}.(cfg.metric{irow})(index, :)]);
                    end
                end
                y2 = y;
                %find first non zero data because of the log
                if istrue(cfg.log(irow)) && y2(1) == 0
                    h = findobj(gca,'Type','line');
                    clear ytemp
                    for iline = 1:size(h, 1)
                        ytemp{iline} = min(h(iline).YData(h(iline).YData>0));
                    end
                    if isempty([ytemp{:}]) %means that all values are equal to zero
                        y = ylim;
                        y2 = y;
                    else
                        y2(1) = min([ytemp{:}]);
                    end
                end
                if y(1) == inf
                    y = [-inf inf];
                end
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y2(1) y2(1) y2(2) y2(2)], 'r', 'edgecolor', 'r', 'facealpha', 0.4, 'edgealpha', 0.4);
                end
                ylim(y);
                
                % replace legend with added artefact
                if strcmp(cfg.metric{irow}, 'dist')
                    if exist('l2', 'var')
                        delete(hl);
                        hl = legend([l2; "artefact"], 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
                    end
                else
                    delete(hl);
                    hl = legend({l{cfg.unit{irow}}, 'artefact'}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
                end
                pos_legend = get(hl, 'position');
                pos_legend(1) = legend_x;
                set(hl, 'position', pos_legend);
            end
            
            if cfg.log(irow)
                set(gca, 'YScale', 'log')
            end
            
        case 'trialinfo'
            
            
            time = SpikeTrials{ipart}.(markername).trialinfo.starttime + (SpikeTrials{ipart}.(markername).trialinfo.endtime - SpikeTrials{ipart}.(markername).trialinfo.starttime) / 2;
            if ~isempty(cfg.offset)
                time = seconds(time - cfg.offset);
            end
            
            % in case a second index is given, e.g. with Spiky distance
            if isempty(cfg.index{irow}); cfg.index{irow} = 1; end
            
            if cfg.index{irow} == 0; cfg.index{irow} = 1; end
            
            if cfg.index{irow} == 1
                bar(time, SpikeTrials{ipart}.(markername).trialinfo.(cfg.metric{irow}), 1, 'k');
                try
                    plot(time(cfg.marker_indx{irow}), SpikeTrials{ipart}.(markername){cfg.unit{irow}}.(cfg.metric{irow})(cfg.marker_indx{irow}), cfg.marker_sign{irow});
                catch
                end
            else
                bar(time, SpikeTrials{ipart}.(markername).trialinfo.(cfg.metric{irow})(cfg.index{irow}, :), 1, 'k');
                try
                    plot(time(cfg.marker_indx{irow}), SpikeTrials{ipart}.(markername).(cfg.metric{irow})(indx, cfg.marker_indx{irow}), cfg.marker_sign{irow});
                catch
                end
            end
            
            if cfg.plotart(irow)
                y = [min(SpikeTrials{ipart}.(markername).trialinfo.(cfg.metric{irow})), max(SpikeTrials{ipart}.(markername).trialinfo.(cfg.metric{irow}))];
                for ibad = 1 : size(bad_t1, 2)
                    fill([bad_t1(ibad), bad_t2(ibad), bad_t2(ibad), bad_t1(ibad)],[y(1) y(1) y(2) y(2)], 'r', 'edgecolor', 'r', 'facealpha', 0.4, 'edgealpha', 0.4);
                end
                ylim(y);
            end
            
            hl = legend(cfg.metric{irow}, 'location', 'eastoutside', 'Interpreter', 'none', 'AutoUpdate','off');
            pos_legend = get(hl, 'position');
            pos_legend(1) = legend_x;
            set(hl, 'position', pos_legend);
            
        otherwise
            error('Cannot recognize type %s', cfg.type{irow});
            %             continue
    end
    
    if isfield(cfg, 'xlim')
        xlim(cfg.xlim);
    end
    
    if irow == nrows
        s1.XAxis.MinorTick = 'on';
        if isdatetime(cfg.xlim)
            s1.XAxis.MinorTickValues = dateshift(cfg.xlim(1) : minutes(60) : cfg.xlim(2), 'start', 'hour');
        end
    end
    
    try
        if cfg.PSGcolors(irow) && ~strcmp(cfg.type{irow}, 'hypnogram')
            plot_hyp_colors(hypnogram, cm_hyp, ylim);
            set(gca, 'children', flipud(get(gca, 'children'))); % flip order of images
        end
    catch
    end
    
    if cfg.log(irow)
        set(gca, 'YScale', 'log')
    end
    
    
%     if irow == 1
%         xleft  = time(1);
%         xright = time(end);
%     else
%         xleft   = max(time(1), xleft);
%         xright  = min(time(end), xright);
%     end
%     xlim([xleft, xright]);
%     
end

isdir_or_mkdir(fileparts(cfg.filename));
exportgraphics(fig, cfg.filename);

function plot_hyp_colors(h, cm, y)

% plot hypnogram in colors
for im = 1 : height(h)
    switch h.hyplabel{im}
        case 'NO_SCORE'
            ci = 6;
        case 'PRESLEEP'
            ci = 6;
        case 'POSTSLEEP'
            ci = 6;            
        case 'REM'
            ci = 1;
        case 'AWAKE'
            ci = 2;
        case 'PHASE_1'
            ci = 3;
        case 'PHASE_2'
            ci = 4;
        case 'PHASE_3'
            ci = 5;
    end
    fill([h.starttime(im), h.endtime(im), h.endtime(im), h.starttime(im)],[y(1) y(1) y(2) y(2)], cm(ci, :) ,'edgecolor', 'none', 'facealpha', 0.4);
end
