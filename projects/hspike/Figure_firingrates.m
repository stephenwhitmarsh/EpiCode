function Figure_firingrates


restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'))
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end

ft_defaults

config = hspike_setparams;

% load data
for ipatient = 1 : 8
    config{ipatient}.spike.name  = ["template1", "template2", "template3", "template4", "template5", "template6"];
    SpikeTrials{ipatient}        = readSpikeTrials(config{ipatient});
    config{ipatient}.spike.name  = ["window"];
    SpikeStats{ipatient}         = spikeTrialStats(config{ipatient});
    config{ipatient}.LFP.name    = ["template1", "template2", "template3", "template4", "template5", "template6"];
    SpikeDensity{ipatient}       = spikeTrialDensity(config{ipatient});
    LFPavg{ipatient}             = readLFPavg(config{ipatient});
    SpikeWaveforms{ipatient}     = readSpikeWaveforms(config{ipatient});
end

%% rereference template average
for ipatient = 1:8
    for ipart = 1 : size(LFPavg{ipatient}, 2)
        for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
            
            if ~isfield(LFPavg{ipatient}{ipart}, markername)
                continue
            end
            if isempty(LFPavg{ipatient}{ipart}.(markername))
                continue
            end
            LFPavg_reref{ipatient}{ipart}.(markername) = LFPavg{ipatient}{ipart}.(markername);
            
            if strcmp(config{ipatient}.template.reref, 'yes') & ~contains(LFPavg{ipatient}{ipart}.(markername).label{1}, '-')
                
                labels_nonum    = regexprep(LFPavg{ipatient}{ipart}.(markername).label, '[0-9_]', '');
                [~, ~, indx]    = unique(labels_nonum);
                
                % average per part (day) then reref
                clear group
                for i = 1 : max(indx)
                    cfgtemp             = [];
                    cfgtemp.reref       = 'yes';
                    cfgtemp.refmethod   = 'bipolar';
                    cfgtemp.channel     = LFPavg{ipatient}{ipart}.(markername).label(indx==i);
                    group{i}            = ft_preprocessing(cfgtemp, LFPavg{ipatient}{ipart}.(markername));
                end
                LFPavg_reref{ipatient}{ipart}.(markername) = ft_appenddata([], group{:});
            end
        end
    end
end
clear LFPavg

% configure figure
papersize       = 1000;
nLFPtemplate    = 6;

for ipatient = 1 : 8
    
    for ipart = 1  : 3
        
        fig = figure('visible', true);
        set(fig, 'PaperPositionMode', 'auto');
        %     set(fig, 'position', get(0,'ScreenSize'));
        %     set(fig, 'position', [20 -60  papersize papersize*sqrt(2)]);
        set(fig, 'Renderer', 'Painters');
        
        for iLFPtemplate = 1 : 6 % LFP clusters
            
            markername = sprintf('template%d', iLFPtemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername) || isempty(SpikeTrials{ipatient}{ipart}.(markername))
                continue
            end
            
            % configure number and dimensions of plots
            ntemplates  = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2);
            w           = 1/nLFPtemplate * 0.9;
            hratio      = 1/nLFPtemplate;
            vratio      = 1/(ntemplates+1) ;
            htop        = 1/(ntemplates+1) * 0.8;
            hbottom     = 1/(ntemplates+1) * 0.8;
            spacehor    = 0.1;
            rshift      = 0.015;
            upshift     = 0.01;
            w           = w * (1-spacehor);
            col         = 1;
            
            % plot LFP
            row     = 1;
            s0      = axes('Position', [hratio*(iLFPtemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
            set(s0, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');
            maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{1}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');
            plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), 'k')
            axis tight;
            xlim(config{ipatient}.spike.toi.(markername));
            y = ylim;
            yticks(floor(y(end)));
            xticks([]);
            
            for itemp = 1 : ntemplates % spike clusters
                
                row     = itemp + 1;
                s1      = axes('Position', [hratio*(iLFPtemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
                set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');
                
                hold on
                if size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2) == 1
                    bar(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time, SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                else
                    bar(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time, SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                end
                
                if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'posclusters')
                    for ipos = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters, 2)
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < config{ipatient}.stats.alpha
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
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < config{ipatient}.stats.alpha
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
                xlim(config{ipatient}.spike.toi.(markername));
                y = ylim;
                yticks(floor(y(end)));
                
                % baseline line
                if itemp == ntemplates
                    y = ylim;
                    y = y(1) + (y(2) - y(1)) * 0.05;
                    x = xlim;
                    width = x(2) - x(1);
                    line(config{ipatient}.stats.bl.(markername), [y, y], 'linewidth', 2, 'color', 'k');
                    text(config{ipatient}.stats.bl.(markername)(1) + width * 0.01, y, 'Baseline', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                    
                    % time line
                    duration = 0.5;
                    t = sprintf('500 ms');
                    line([config{ipatient}.spike.toi.(markername)(end) - duration, config{ipatient}.spike.toi.(markername)(end)], [y, y], 'linewidth', 2, 'color', 'k');
                    text(config{ipatient}.spike.toi.(markername)(end) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
                    y = ylim;
                end
                
                % label MUA/SUA
                if contains(SpikeTrials{ipatient}{ipart}.(markername).cluster_group{itemp}, 'good')
                    text(config{ipatient}.stats.bl.(markername)(1) + width * 0.91, y(2) * 0.95, "SUA", 'color', 'k', 'fontsize', 8);
                else
                    text(config{ipatient}.stats.bl.(markername)(1) + width * 0.91, y(2) * 0.95, "MUA", 'color', 'k', 'fontsize', 8); %0.9
                end
                
                % add inset with waveshape
                if iLFPtemplate == 1
                    
                    y_all   = vertcat(SpikeWaveforms{ipatient}{ipart}{itemp}.trial{:})';
                    r       = rms(y_all - mean(y_all, 2));
                    [~, i]  = sort(r, 'ascend');
                    y_sel   = y_all(:, i(1:100));
                    
                    pos = s1.Position;
                    pos(1) = 1 / nLFPtemplate -rshift - 0.02;
                    pos(2) = pos(2) + 0.020; % 0.015
                    pos(3) = pos(3) / 6;
                    pos(4) = pos(4) / 1.5;
                    inset1 = axes('position', pos, 'color', 'none', 'XColor','none', 'YColor', 'none'); hold on; %left bottom width height
                    
                    for i = 1 : size(y_sel, 2)
                        lh = plot(y_sel(:, i), 'k'); lh.Color = [0, 0, 0, 0.2];
                    end
                    plot(mean(y_all, 2), 'k', 'linewidth', 2);
                    plot(mean(y_all, 2), 'w', 'linewidth', 1);
                end
                
                % add inset with ISI histogram
                if iLFPtemplate == 1
                    pos = s1.Position;
                    pos(1) = 1 / nLFPtemplate -rshift - 0.05;
                    pos(2) = pos(2) + 0.020; % 0.015
                    pos(3) = pos(3) / 6;
                    pos(4) = pos(4) / 1.5;
                    inset2 = axes('position', pos, 'color', 'none', 'YColor','none', 'Tickdir', 'out'); hold on;
                    bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0, 0.5, 0], 'edgecolor', 'none');
                    x = find(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time < 0.002, 1, 'last');
                    bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time(1:x) * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg(1:x), 1, 'facecolor', [0.5, 0, 0], 'edgecolor', 'none');
                    xlim([0, 30]);
                    xticks([0, 30]);
                    l = xlabel('time');
                    x = get(l,'Position');
                    set(l, 'Position', x .* [1, 0.1, 1]);
                end
            end
        end % iLFPtemplate
        
        fname = fullfile(config{ipatient}.imagesavedir, sprintf('firingrates_patient%d_part%d', ipatient, ipart));
        % exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
        exportgraphics(fig, strcat(fname, '.pdf'));
    end % ipart
end % ipatient


%% normalize mean
for ipatient = 1 : 8
    for ipart = 1  : 3
        for iLFPtemplate = 1 : 6 % LFP clusters
            
            markername = sprintf('template%d', iLFPtemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername)
                continue
            end
            
            toinclude{ipatient}{ipart}.(markername) = [];
            
            % configure number and dimensions of plots
            ntemplates  = size(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).label, 2);
            SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg_norm = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg;
            for itemp = 1 : ntemplates
                
                if isempty(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername))
                    continue
                end
                if isempty(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername))
                    continue
                end
                t =  SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time >= config{ipatient}.stats.bl(1).(markername)(1) & SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time <= config{ipatient}.stats.bl(1).(markername)(2);
                if ntemplates > 1
                    SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg_norm(itemp, :) = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, :) ./ mean(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(itemp, t));
                else
                    SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg_norm = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg ./ mean(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(t));
                    SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg_norm = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg_norm';
                end
                
                if isempty(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp})
                    continue
                end
                if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.responsive_pos || SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.responsive_neg
                    toinclude{ipatient}{ipart}.(markername) = [toinclude{ipatient}{ipart}.(markername), itemp];
                end
            end
        end
    end
end

for ipatient = 1 : 8
    for ipart = 1  : 3
        SpikeDensity{ipatient}{ipart}.sdf_bar_avg = SpikeDensity{ipatient}{ipart}.sdf_bar;
        for iLFPtemplate = 1 : 6 % LFP clusters
            markername = sprintf('template%d', iLFPtemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername)
                continue
            end
            try
                SpikeDensity{ipatient}{ipart}.psth.(markername) = rmfield(SpikeDensity{ipatient}{ipart}.psth.(markername),'cfg');
            catch
            end
            try
                SpikeDensity{ipatient}{ipart}.sdf_lin.(markername) = rmfield(SpikeDensity{ipatient}{ipart}.sdf_lin.(markername),'cfg');
            catch
            end
            try
                SpikeDensity{ipatient}{ipart}.sdf_bar.(markername) = rmfield(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername),'cfg');
            catch
            end     
            try
                for itemp = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername), 2)
                    SpikeDensity{ipatient}{ipart}.stat.(markername){itemp} = rmfield( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp},'cfg');
                end
            catch
            end
            
            SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername)          = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername);
            SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername).label    = {'active'};
            if size(toinclude{ipatient}{ipart}.(markername), 2) > 1
                SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername).avg      = mean(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(toinclude{ipatient}{ipart}.(markername), :));
            else
                SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername).avg      = SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).avg(toinclude{ipatient}{ipart}.(markername), :);
            end
        end
    end
end




%% plot normalized

fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
%     set(fig, 'position', get(0,'ScreenSize'));
%     set(fig, 'position', [20 -60  papersize papersize*sqrt(2)]);
set(fig, 'Renderer', 'Painters');

% configure figure
papersize       = 1000;
nLFPtemplate    = 6;

for ipatient = 1 : 8
    
    %     for ipart = 1  : 3
    
    
    for iLFPtemplate = 1 : 6 % LFP clusters
        
        markername = sprintf('template%d', iLFPtemplate);
        if ~isfield(SpikeDensity{ipatient}{ipart}.sdf_bar, markername) || isempty(SpikeTrials{ipatient}{ipart}.(markername)) || isempty(SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername).avg)
            continue
        end
                
        subplot(8, 6, iLFPtemplate + (ipatient-1) * 6);

        yyaxis left
        bar(SpikeDensity{ipatient}{ipart}.sdf_bar.(markername).time, SpikeDensity{ipatient}{ipart}.sdf_bar_avg.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
        xticks([]);
        y = ylim;
        yticks(floor(y(end)));

        yyaxis right
        hold;
        % plot LFP
        try
            maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');
            plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), 'k')
            axis tight;
            xlim(config{ipatient}.spike.toi.(markername));
            y = ylim;
            yticks(floor(y(end)));
            xticks([]);
        catch
        end

        axis tight;
        xlim(config{ipatient}.spike.toi.(markername));
        y = ylim;
        yticks(floor(y(end)));
        
        % baseline line
        if ipatient == 8
            y = ylim;
            y = y(1) + (y(2) - y(1)) * 0.05;
            x = xlim;
            width = x(2) - x(1);
            line(config{ipatient}.stats.bl.(markername), [y, y], 'linewidth', 2, 'color', 'k');
            text(config{ipatient}.stats.bl.(markername)(1) + width * 0.01, y, 'Baseline', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            
            % time line
            duration = 0.5;
            t = sprintf('500 ms');
            line([config{ipatient}.spike.toi.(markername)(end) - duration, config{ipatient}.spike.toi.(markername)(end)], [y, y], 'linewidth', 2, 'color', 'k');
            text(config{ipatient}.spike.toi.(markername)(end) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
            y = ylim;
        end
        
        
        
        %         end % iLFPtemplate
        
        %         fname = fullfile(config{ipatient}.imagesavedir, sprintf('firingrates_patient%d_part%d', ipatient, ipart));
        % exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
        %     exportgraphics(fig, strcat(fname, '.pdf'));
    end % ipart
end % ipatient





