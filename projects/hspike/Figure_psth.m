function Figure_psth


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
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer
    
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
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer
    
end

ft_defaults
config = hspike_setparams;

%% load data
for ipatient = 1 : 8
    config{ipatient}.spike.name  = ["template1", "template2", "template3", "template4", "template5", "template6"];
    SpikeTrials{ipatient}        = readSpikeTrials(config{ipatient});
    
    config{ipatient}.spike.name  = ["window"];
    SpikeStats{ipatient}         = spikeTrialStats(config{ipatient}); %% formally: spikeTrialDensity(config{ipatient});
    
    config{ipatient}.LFP.name    = ["template1", "template2", "template3", "template4", "template5", "template6"];
    config{ipatient}.LFP.postfix = {'_all'};
    LFPavg{ipatient}             = readLFPavg(config{ipatient});
    SpikeDensity{ipatient}       = spikePSTH(config{ipatient});    
    SpikeWaveforms{ipatient}     = readSpikeWaveforms(config{ipatient});
end

%% rereference LFP average
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

%% demean/baseline correct LFP average
for ipatient = 1 : 8
    for ipart = 1 : size(LFPavg_reref{ipatient}, 2)
        for template = ["template1","template2","template3","template4","template5","template6"]
            if ~isfield(LFPavg_reref{ipatient}{ipart}, template)
                continue
            end
            if isempty(LFPavg_reref{ipatient}{ipart}.(template))
                continue
            end
            cfgtemp                 = [];
            cfgtemp.demean          = 'yes';
            cfgtemp.baselinewindow  = [-0.15, -0.05];
            LFPavg_reref{ipatient}{ipart}.(template) = ft_preprocessing(cfgtemp, LFPavg_reref{ipatient}{ipart}.(template));
        end
    end
end

%% select latency LFP
latency = [-0.2 0.5];
cfgtemp = [];
cfgtemp.latency = latency;
for ipatient = 1 : 8
    for ipart = 1 : size(LFPavg_reref{ipatient}, 2)
        for template = ["template1","template2","template3","template4","template5","template6"]
            if ~isfield(LFPavg_reref{ipatient}{ipart}, template)
                continue
            end
            if isempty(LFPavg_reref{ipatient}{ipart}.(template))
                continue
            end
            LFPavg_reref{ipatient}{ipart}.(template) = ft_selectdata(cfgtemp, LFPavg_reref{ipatient}{ipart}.(template));
        end
    end
end

%% flip polarity of rereferenced LFP patient 6
ipatient = 6;
for ipart = 1 : size(LFPavg_reref{ipatient}, 2)
    for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
        
        if ~isfield(LFPavg_reref{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFPavg_reref{ipatient}{ipart}.(markername))
            continue
        end
        LFPavg_reref{ipatient}{ipart}.(markername).trial{1} = -LFPavg_reref{ipatient}{ipart}.(markername).trial{1};
        
    end
end

%% plot each unit

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
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername) || isempty(SpikeTrials{ipatient}{ipart}.(markername))
                continue
            end
            
            % configure number and dimensions of plots
            ntemplates  = size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2);
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
            try
                
                maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{1}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');
                plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), 'k')
            catch
            end
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
                if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                    bar(SpikeDensity{ipatient}{ipart}.psth.(markername).time, SpikeDensity{ipatient}{ipart}.psth.(markername).avg, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                else
                    bar(SpikeDensity{ipatient}{ipart}.psth.(markername).time, SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                end
                
                if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'posclusters')
                    for ipos = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters, 2)
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < config{ipatient}.stats.alpha
                            sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
                            if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg( sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                end
                            else
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                end
                            end
                        end
                    end
                end
                
                if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'negclusters')
                    for ineg = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters, 2)
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < config{ipatient}.stats.alpha
                            sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                            if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg( sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                end
                            else
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                end
                            end
                        end
                    end
                end
                
                axis tight;
%                 xlim(config{ipatient}.spike.toi.(markername));
                xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);

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
                if iLFPtemplate == 1
                    
                    x = xlim;
                    width = x(2) - x(1);
                    if contains(SpikeTrials{ipatient}{ipart}.(markername).cluster_group{itemp}, 'good')
                        text(config{ipatient}.stats.bl.(markername)(1) + width * 0.81, y(2) * 0.95, "SUA", 'color', 'k', 'fontsize', 8);
                    else
                        text(config{ipatient}.stats.bl.(markername)(1) + width * 0.81, y(2) * 0.95, "MUA", 'color', 'k', 'fontsize', 8); %0.9
                    end
                end
                
%                 % add inset with waveshape
%                 if iLFPtemplate == 1
%                     
%                     y_all   = vertcat(SpikeWaveforms{ipatient}{ipart}{itemp}.trial{:})';
%                     r       = rms(y_all - mean(y_all, 2));
%                     [~, i]  = sort(r, 'ascend');
%                     y_sel   = y_all(:, i(1:100));
%                     
%                     pos = s1.Position;
%                     pos(1) = 1 / nLFPtemplate -rshift - 0.02;
%                     pos(2) = pos(2) + 0.020; % 0.015
%                     pos(3) = pos(3) / 6;
%                     pos(4) = pos(4) / 1.5;
%                     inset1 = axes('position', pos, 'color', 'none', 'XColor','none', 'YColor', 'none'); hold on; %left bottom width height
%                     
%                     for i = 1 : size(y_sel, 2)
%                         lh = plot(y_sel(:, i), 'k'); lh.Color = [0, 0, 0, 0.2];
%                     end
%                     plot(mean(y_all, 2), 'k', 'linewidth', 2);
%                     plot(mean(y_all, 2), 'w', 'linewidth', 1);
%                 end
%                 
%                 % add inset with ISI histogram
%                 if iLFPtemplate == 1
%                     pos = s1.Position;
%                     pos(1) = 1 / nLFPtemplate -rshift - 0.05;
%                     pos(2) = pos(2) + 0.020; % 0.015
%                     pos(3) = pos(3) / 6;
%                     pos(4) = pos(4) / 1.5;
%                     inset2 = axes('position', pos, 'color', 'none', 'YColor','none', 'Tickdir', 'out'); hold on;
%                     bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0, 0.5, 0], 'edgecolor', 'none');
%                     x = find(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time < 0.002, 1, 'last');
%                     bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time(1:x) * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg(1:x), 1, 'facecolor', [0.5, 0, 0], 'edgecolor', 'none');
%                     xlim([0, 30]);
%                     xticks([0, 30]);
%                     l = xlabel('time');
%                     x = get(l,'Position');
%                     set(l, 'Position', x .* [1, 0.1, 1]);
%                 end
            end
        end % iLFPtemplate
        
        fname = fullfile(config{ipatient}.imagesavedir, sprintf('firingrates_patient%d_part%d', ipatient, ipart));
        % exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
        exportgraphics(fig, strcat(fname, '.pdf'));
        close all
    end % ipart
end % ipatient

%% standardize PSTH & find responsive neurons
for ipatient = 1 : 8
    for ipart = 1 : 3
        for iLFPtemplate = 1 : 6 % LFP clusters
            
            markername = sprintf('template%d', iLFPtemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername)
                continue
            end
            
            toinclude{ipatient}{ipart}.(markername) = [];
            
            % configure number and dimensions of plots
            ntemplates  = size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2);
            SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm = SpikeDensity{ipatient}{ipart}.psth.(markername).avg;
            for itemp = 1 : ntemplates

                if isempty(SpikeDensity{ipatient}{ipart}.psth.(markername))
                    continue
                end
                
                % select baseline
                t = SpikeDensity{ipatient}{ipart}.psth.(markername).time >= config{ipatient}.stats.bl(1).(markername)(1) & SpikeDensity{ipatient}{ipart}.psth.(markername).time <= config{ipatient}.stats.bl(1).(markername)(2);
                if ntemplates == 1
                    temp    = SpikeDensity{ipatient}{ipart}.psth.(markername).avg(t);
                else
                    temp    = SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, t);
                end
                temp(temp == 0) = nan;
                bl = nanmean(temp);
                
                if ntemplates > 1
                    SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, :) = SpikeDensity{ipatient}{ipart}.psth.(markername).avg(itemp, :) ./ bl;
                else
                    SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm = (SpikeDensity{ipatient}{ipart}.psth.(markername).avg ./ bl)';
                end
                if isempty(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp})
                    continue
                end
                if (SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.responsive_pos || SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.responsive_neg)
                    toinclude{ipatient}{ipart}.(markername) = [toinclude{ipatient}{ipart}.(markername), itemp];
                end
            end
        end
    end
end

%% plot normalized PSTH for each patient/part

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
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername) || isempty(SpikeTrials{ipatient}{ipart}.(markername))
                continue
            end
            
            % configure number and dimensions of plots
            ntemplates  = size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2);
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
            try
                maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{1}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');
                plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), 'k')
            catch
            end
            axis tight;
            xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);
            
            y = ylim;
            yticks(floor(y(end)));
            xticks([]);
            
            for itemp = 1 : ntemplates % spike clusters
                
                row     = itemp + 1;
                s1      = axes('Position', [hratio*(iLFPtemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
                set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');
                
                hold on
                if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                    bar(SpikeDensity{ipatient}{ipart}.psth.(markername).time, SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                else
                    bar(SpikeDensity{ipatient}{ipart}.psth.(markername).time, SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, :), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
                end
                
                if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'posclusters')
                    for ipos = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters, 2)
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < config{ipatient}.stats.alpha
                            sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos);
                            if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm, 1) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm( sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                end
                            else
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag)], 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag), 1, 'facecolor', [252/255,187/255,62/255], 'edgecolor', 'none');
                                end
                            end
                        end
                    end
                end
                
                if isfield(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}, 'negclusters')
                    for ineg = 1 : size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters, 2)
                        if SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < config{ipatient}.stats.alpha
                            sel = find(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg);
                            if size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2) == 1
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm, 1) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm( sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                end
                            else
                                lag = size(SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm, 2) - size(SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.mask, 2);
                                if length(sel) == 1
                                    bar([SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)-0.001, SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel)+0.01], [SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag)], 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                else
                                    bar( SpikeDensity{ipatient}{ipart}.stat.(markername){itemp}.time(sel), SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(itemp, sel+lag), 1, 'facecolor', [70/255,93/255,250/255], 'edgecolor', 'none');
                                end
                            end
                        end
                    end
                end
                
                axis tight;
                xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);

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
                if iLFPtemplate == 1
                    
                    x = xlim;
                    width = x(2) - x(1);
                    if contains(SpikeTrials{ipatient}{ipart}.(markername).cluster_group{itemp}, 'good')
                        text(config{ipatient}.stats.bl.(markername)(1) + width * 0.81, y(2) * 0.95, "SUA", 'color', 'k', 'fontsize', 8);
                    else
                        text(config{ipatient}.stats.bl.(markername)(1) + width * 0.81, y(2) * 0.95, "MUA", 'color', 'k', 'fontsize', 8); %0.9
                    end
                end
                
%                 % add inset with waveshape
%                 if iLFPtemplate == 1
%                     
%                     y_all   = vertcat(SpikeWaveforms{ipatient}{ipart}{itemp}.trial{:})';
%                     r       = rms(y_all - mean(y_all, 2));
%                     [~, i]  = sort(r, 'ascend');
%                     y_sel   = y_all(:, i(1:100));
%                     
%                     pos = s1.Position;
%                     pos(1) = 1 / nLFPtemplate -rshift - 0.02;
%                     pos(2) = pos(2) + 0.020; % 0.015
%                     pos(3) = pos(3) / 6;
%                     pos(4) = pos(4) / 1.5;
%                     inset1 = axes('position', pos, 'color', 'none', 'XColor','none', 'YColor', 'none'); hold on; %left bottom width height
%                     
%                     for i = 1 : size(y_sel, 2)
%                         lh = plot(y_sel(:, i), 'k'); lh.Color = [0, 0, 0, 0.2];
%                     end
%                     plot(mean(y_all, 2), 'k', 'linewidth', 2);
%                     plot(mean(y_all, 2), 'w', 'linewidth', 1);
%                 end
%                 
%                 % add inset with ISI histogram
%                 if iLFPtemplate == 1
%                     pos = s1.Position;
%                     pos(1) = 1 / nLFPtemplate -rshift - 0.05;
%                     pos(2) = pos(2) + 0.020; % 0.015
%                     pos(3) = pos(3) / 6;
%                     pos(4) = pos(4) / 1.5;
%                     inset2 = axes('position', pos, 'color', 'none', 'YColor','none', 'Tickdir', 'out'); hold on;
%                     bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg, 1, 'facecolor', [0, 0.5, 0], 'edgecolor', 'none');
%                     x = find(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time < 0.002, 1, 'last');
%                     bar(SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg_time(1:x) * 1000, SpikeStats{ipatient}{ipart}.window{itemp}.isi_avg(1:x), 1, 'facecolor', [0.5, 0, 0], 'edgecolor', 'none');
%                     xlim([0, 30]);
%                     xticks([0, 30]);
%                     l = xlabel('time');
%                     x = get(l,'Position');
%                     set(l, 'Position', x .* [1, 0.1, 1]);
%                 end
            end
        end % iLFPtemplate
        
        fname = fullfile(config{ipatient}.imagesavedir, sprintf('firingrates_patient%d_part%d_normalized', ipatient, ipart));
        exportgraphics(fig, strcat(fname, '.pdf'));
    end % ipart
end % ipatient

%% average standardized PSTH of responsive neurons 
for ipatient = 1 : 8
    for ipart = 1  : 3
        SpikeDensity{ipatient}{ipart}.psth_avg = [];
        for iLFPtemplate = 1 : 6 
            markername = sprintf('template%d', iLFPtemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername)
                continue
            end
            if isempty(toinclude{ipatient}{ipart}.(markername))
                continue
            end
            if size(toinclude{ipatient}{ipart}.(markername), 2) > 1
                SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg = nanmean(SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(toinclude{ipatient}{ipart}.(markername), :));
            elseif toinclude{ipatient}{ipart}.(markername) ~= 1
                SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg = SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm(toinclude{ipatient}{ipart}.(markername), :);
            elseif toinclude{ipatient}{ipart}.(markername) == 1
                SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg = SpikeDensity{ipatient}{ipart}.psth.(markername).avg_norm;
            end
        end
    end
end

%% average LFP over nights 
SpikeDensity_avg = [];
clear temp_time temp_avg
for ipatient = 1 : 8
    for iLFPtemplate = 1 : 6
        markername = sprintf('template%d', iLFPtemplate);
        temp_avg = [];
        for ipart = 1 : 3
            
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth_avg, markername)
                continue
            end
            
            d = size(SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg);
            if d(1) > d(2)
                temp_avg = [temp_avg; SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg'];
            else
                temp_avg = [temp_avg; SpikeDensity{ipatient}{ipart}.psth_avg.(markername).avg];
            end
        end
        
        if size(temp_avg, 1) > 1
            SpikeDensity_avg{ipatient}.(markername) = nanmean(temp_avg);
        end
    end
end

%% plot averaged normalized PSTH of responsive neurons for each template
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
%     set(fig, 'position', get(0,'ScreenSize'));
%     set(fig, 'position', [20 -60  papersize papersize*sqrt(2)]);
set(fig, 'Renderer', 'Painters');

% configure figure
papersize       = 1000;
nLFPtemplate    = 6;

for ipatient = 1 : 8
   
    for iLFPtemplate = 1 : 6 % LFP clusters
        
        markername = sprintf('template%d', iLFPtemplate);
        if ~isfield(SpikeDensity_avg{ipatient}, markername) || isempty(SpikeDensity_avg{ipatient}.(markername))
            continue
        end
        
        subplot(8, 6, iLFPtemplate + (ipatient-1) * 6);
        
        yyaxis left
        bar(SpikeDensity{ipatient}{1}.psth.(markername).time, SpikeDensity_avg{ipatient}.(markername), 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
        y = ylim;
        yticks(sort(unique([-floor(y(1)), 1, floor(y(2))])));
        xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);

        % labels
        if iLFPtemplate == 1 && ipatient ==8
            xlabel('Time');
            ylabel('Change vs. baseline');
        end
        
        yyaxis right
        hold;
        
        % plot LFP
        if isfield(LFPavg_reref{ipatient}{ipart}, markername)
            maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');
            plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), 'k')
            axis tight;
            y = ylim;
            yticks(floor(y(end)));

            axis tight;
            xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);
            y = ylim;
            yticks(sort([floor(y(1)), 0, floor(y(2))]));
        end
        
        % labels
        if iLFPtemplate == 1 && ipatient ==8
            ylabel('Amplitude (mV)');
        end
        
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
        end
    end
end % ipatient

fname = fullfile(config{ipatient}.imagesavedir, 'firingrates_patients');
% exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));

%% average PSTH over templates

SpikeDensity_avg_avg = [];
for ipatient = 1 : 8
    temp_avg = [];
    for iLFPtemplate = 1 : 6
        markername = sprintf('template%d', iLFPtemplate);
        if isfield(SpikeDensity_avg{ipatient}, markername)
            temp_avg = [temp_avg; SpikeDensity_avg{ipatient}.(markername)];
        end
    end
    if size(temp_avg, 1) > 1
        SpikeDensity_avg_avg{ipatient} = nanmean(temp_avg);
    end
end

%% plot average normalized averages per template of responsive neurons
latency = [-0.2 0.5];

fig = figure('visible', true);
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'Renderer', 'Painters');
set(0,'DefaultAxesTitleFontWeight','normal');
ylims = [600, 1000, 200, 700, 1700, 150, 250, 400] / 1000; % go to mV

for ipatient = 1 : 8
    
    s1 = subplot(2, 4, ipatient);
%     set(s, 'XGrid', 'off', 'box', 'off', 'tickdir', 'out', 'Color', 'none');
    
%     yyaxis left
%     bar(SpikeDensity{ipatient}{1}.psth.(markername).time, SpikeDensity_avg_avg{ipatient}, 1, 'facecolor', [127/255,127/255,127/255], 'edgecolor', 'none');
    bar(SpikeDensity{ipatient}{1}.psth.template3.time, SpikeDensity_avg_avg{ipatient}, 1, 'facecolor', 'k', 'edgecolor', 'none');
    set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none');

    % zero line
    hold on
    
    % smooth line
    s = smooth(SpikeDensity{ipatient}{1}.psth.template3.time, smooth(SpikeDensity_avg_avg{ipatient}), 3, 'moving')';
%     plot(temp_time, s, 'k-');
   
    axis tight
%     xlim([config{ipatient}.stats.bl.(markername)(1), config{ipatient}.spike.toi.(markername)(2)]);
    xlim(latency);

    y = ylim;
    if y(2) < 2
        yticks([round(y(2)+0.1,1)]);
        ylim([0, round(y(2)+0.1,1)]);
    else
        yticks([round(y(2))]);
        ylim([0, ceil(y(2))]);
    end
    y = ylim;
    ylim([y(1), y(2)*1.5]);
    if ipatient == 7
        ylim([0, 5]);
        yticks([1]);        
    end
    
    % baseline line
    if ipatient == 4
        y = ylim;
        y = 0.1;
        x = xlim;
        line([-0.3, -0.1], [y, y], 'linewidth', 2, 'color', 'w');
        text(-0.2, y, 'BL', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        
        % time line
        duration = 0.300;
        t = sprintf('300 ms');
        line([x(2) - duration, config{ipatient}.spike.toi.(markername)(end)], [y, y], 'linewidth', 2, 'color', 'w');
        text(x(2) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    
    % labels
    if ipatient == 1 || ipatient == 5
        % xlabel('Time');
        l = ylabel('Change vs. baseline');
        pos = l.Position;
        pos(1) = -0.4;
        l.Position = pos;
    end
    
    % plot average LFP
    yyaxis right
    temp_avg = [];
    for iLFPtemplate = 1:6
        markername = sprintf('template%d', iLFPtemplate);
        if ~isfield(LFPavg_reref{ipatient}{ipart}, markername)
            continue
        end
        maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'last');
        temp_avg(iLFPtemplate,:) = LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:);
    end    
  
    ylim([-ylims(ipatient), ylims(ipatient)]);
    yticks([-ylims(ipatient), 0, ylims(ipatient)]);
    yticklabels([-ylims(ipatient), "", ylims(ipatient)]);
    set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0.03]);
    y = ylim;
    ylim([y(1), y(2)*3]);
    
    plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, nanmean(temp_avg)/1000, 'k', 'linewidth', 1)
    set(gca, 'ydir', 'reverse'); 
    
    if ipatient ~= 7
        % positive area
        temp_time = SpikeDensity{ipatient}{1}.psth.template3.time;
        indx = temp_time > -0.1 & temp_time < 1;
        [Ypk,Xpk,Wpk,Ppk] = findpeaks(s(indx), temp_time(indx), 'SortStr', 'descend', 'NPeaks',1);
        toi_pos{ipatient}(1) = Xpk-Wpk/2;
        toi_pos{ipatient}(2) = Xpk+Wpk/2;
        y = ylim;
        fill([toi_pos{ipatient}(1), toi_pos{ipatient}(2), toi_pos{ipatient}(2), toi_pos{ipatient}(1)],[y(1), y(1), y(2), y(2)], 'r','facealpha',0.1, 'edgecolor', 'none');
        
        % negative area
        indx = temp_time > -0.1 & temp_time < 1;
        [Ypk,Xpk,Wpk,Ppk] = findpeaks(-s(indx), temp_time(indx), 'SortStr', 'descend', 'NPeaks',1);
        toi_neg{ipatient}(1) = max(Xpk-Wpk/2, toi_pos{ipatient}(2));
        toi_neg{ipatient}(2) = Xpk+Wpk/2;
        y = ylim;
        fill([toi_neg{ipatient}(1), toi_neg{ipatient}(2), toi_neg{ipatient}(2), toi_neg{ipatient}(1)],[y(1), y(1), y(2), y(2)], 'b','facealpha',0.1, 'edgecolor', 'none');
    end
    
    if ipatient == 4 || ipatient == 8
        ylabel('Amplitude (mV)');
    end
    
    title(sprintf('Patient %d', ipatient));
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

end % ipatient
set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 24);
fname = fullfile(config{ipatient}.imagesavedir, 'firingrates_all');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'), 'ContentType', 'vector');

%% Make table for R - PSTH table

t = table;
i = 1;
hyplabels = ["PHASE_3", "PHASE_2", "PHASE_1", "AWAKE", "PRESLEEP", "POSTSLEEP", "REM", "NO_SCORE"];

for ipatient = [1:8]
    for ipart = 1 : 3
        for itemplate = 1 : 6
            markername = sprintf('template%d', itemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername)
                continue
            end
            

            
            for iunit = 1 : size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2)
                
                if isempty(SpikeDensity{ipatient}{ipart}.stat.(markername){iunit})
                    continue
                end
                
                for hyplabel = hyplabels

                    t.Patient(i)  = ipatient;
                    t.part(i)     = ipart;
                    t.template(i) = itemplate+ipatient*100;
                    t.unit(i)     = iunit;
                    t.responsive(i) = SpikeDensity{ipatient}{ipart}.stat.(markername){iunit}.responsive_neg || SpikeDensity{ipatient}{ipart}.stat.(markername){iunit}.responsive_pos;
                    if strcmp(strtrim(SpikeTrials{ipatient}{ipart}.(markername).cluster_group{iunit}), 'good')
                        t.SUA(i) = 1;
                    else
                        t.SUA(i) = 0;
                    end
                    
                    
                    if ~isempty(toi_pos{ipatient}) || ~isempty(toi_neg{ipatient})
    
                    indx_pos = SpikeDensity{ipatient}{ipart}.psth.(markername).time > toi_pos{ipatient}(1) & SpikeDensity{ipatient}{ipart}.psth.(markername).time < toi_pos{ipatient}(2);
                    indx_neg = SpikeDensity{ipatient}{ipart}.psth.(markername).time > toi_neg{ipatient}(1) & SpikeDensity{ipatient}{ipart}.psth.(markername).time < toi_neg{ipatient}(2);
                    
                    
                    indx          = SpikeDensity{ipatient}{ipart}.psth.(markername).trialinfo.hyplabel == hyplabel;
                    
                    trlcnt_pos    = permute(SpikeDensity{ipatient}{ipart}.psth.(markername).trial(indx, iunit, indx_pos), [1, 3, 2]);
                    trlavg_pos    = nanmean(trlcnt_pos);
                    avg_pos       = nanmean(trlavg_pos);
                    t.poscnt(i)   = avg_pos;
                    t.posrate(i)  = avg_pos / config{ipatient}.spike.psthbin.(markername);
                    
                    trlcnt_neg    = permute(SpikeDensity{ipatient}{ipart}.psth.(markername).trial(indx, iunit, indx_neg), [1, 3, 2]);
                    trlavg_neg    = nanmean(trlcnt_neg);
                    avg_neg       = nanmean(trlavg_neg);
                    t.negcnt(i)   = avg_neg;
                    t.negrate(i)  = avg_neg / config{ipatient}.spike.psthbin.(markername);        
                    else
                        t.poscnt(i)  = nan;
                        t.posrate(i) = nan;
                        t.negcnt(i)  = nan;
                        t.negrate(i) = nan;
                    end
                    
                    switch hyplabel
                        case "PHASE_1"
                            t.hyplabel(i) = "S1";
                        case "PHASE_2"
                            t.hyplabel(i) = "S2";
                        case "PHASE_3"
                            t.hyplabel(i) = "S3";
                        case "AWAKE"
                            t.hyplabel(i) = "Wake";
                        case "REM"
                            t.hyplabel(i) = "REM";                            
                        case "POSTSLEEP"
                            t.hyplabel(i) = "Post";
                        case "PRESLEEP"
                            t.hyplabel(i) = "Pre";
                        case "NO_SCORE"
                            t.hyplabel(i) = "NO_SCORE";
                    end
                    
                    i = i + 1;
                end
            end
        end
    end
end

% save data to table for R
fname   = fullfile(config{ipatient}.datasavedir, 'psth_table');
writetable(t, fname);


%% Make table for R (correlation between firing rates and LFP)
t = table;
i = 1;

for ipatient = 1:8
    for ipart = 1 : 3
        for itemplate = 1 : 6
            
            markername = sprintf('template%d', itemplate);
            if ~isfield(SpikeDensity{ipatient}{ipart}.psth, markername)
                continue
            end
            
            temp_avg = [];
            for iLFPtemplate = 1:6
                markername = sprintf('template%d', iLFPtemplate);
                if ~isfield(LFPavg_reref{ipatient}{ipart}, markername)
                    continue
                end
                maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'last');
                temp_avg(iLFPtemplate,:) = LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:);
            end
            
            if ~isfield(LFPavg_reref{ipatient}{ipart}, markername)
                continue
            end
            
            % resample to same time-axis
            maxchan         = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'first');

            cfgtemp = [];
%             cfgtemp.latency = [-0.1, 0.5];
            temp = ft_selectdata(cfgtemp, LFPavg_reref{ipatient}{ipart}.(markername));
            
            cfgtemp         = [];
            cfgtemp.time{1} = SpikeDensity{ipatient}{ipart}.psth.(markername).time;
            LFP_ds          = ft_resampledata(cfgtemp, temp);
            
            for iunit = 1 : size(SpikeDensity{ipatient}{ipart}.psth.(markername).label, 2)
                if isempty(SpikeDensity{ipatient}{ipart}.stat.(markername){iunit})
                    continue
                end
                
                t.Patient(i)  = ipatient;
                t.part(i)     = ipart;
                t.template(i) = itemplate+ipatient*100;
                t.unit(i)     = iunit;
                t.responsive(i) = SpikeDensity{ipatient}{ipart}.stat.(markername){iunit}.responsive_neg || SpikeDensity{ipatient}{ipart}.stat.(markername){iunit}.responsive_pos;
                
                if strcmp(strtrim(SpikeTrials{ipatient}{ipart}.(markername).cluster_group{iunit}), 'good')
                    t.SUA(i) = 1;
                else
                    t.SUA(i) = 0;
                end
                
                if ~isempty(LFPavg_reref{ipatient}{ipart}.(markername))
                    [t.corr_rho(i), t.corr_p(i)]   = corr(LFP_ds.trial{1}(maxchan, :)', SpikeDensity{ipatient}{ipart}.psth.(markername).avg(iunit, :)', 'type', 'pearson');                    
                end
                i = i + 1;
            end
        end
    end
end
                
% save data to table for R
fname   = fullfile(config{ipatient}.datasavedir, 'rho_table');
writetable(t, fname);

%% Combine trial-by-trial LFP and firing rates
% for ipatient = 1 : 8
%     LFP{ipatient} = readLFP(config{ipatient});
% end


disp('done');

%%
% 
%             if isfield(SpikeDensity{ipatient}{ipart}.psth.(markername), "PHASE_1")
%                 LFPstageavg{ipatient}.(template).S1 = LFPstageavg{ipatient}.(template).PHASE_1;
%                 LFPstageavg{ipatient}.(template)    = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_1');
%             end
%             
%             if isfield(LFPstageavg{ipatient}.(template), "PHASE_2")
%                 LFPstageavg{ipatient}.(template).S2      = LFPstageavg{ipatient}.(template).PHASE_2;
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_2');
%             end
%             if isfield(LFPstageavg{ipatient}.(template), "PHASE_3")
%                 LFPstageavg{ipatient}.(template).S3      = LFPstageavg{ipatient}.(template).PHASE_3;
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_3');
%             end
%             
%             if isfield(LFPstageavg{ipatient}.(template), "AWAKE")
%                 LFPstageavg{ipatient}.(template).Wake    = LFPstageavg{ipatient}.(template).AWAKE;
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'AWAKE');
%             end
%             
%             if isfield(LFPstageavg{ipatient}.(template), "POSTSLEEP")
%                 LFPstageavg{ipatient}.(template).Post    = LFPstageavg{ipatient}.(template).POSTSLEEP;
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'POSTSLEEP');
%             end
%             
%             if isfield(LFPstageavg{ipatient}.(template), "PRESLEEP")
%                 LFPstageavg{ipatient}.(template).Pre    = LFPstageavg{ipatient}.(template).PRESLEEP;
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PRESLEEP');
%             end
%             
%             if isfield(LFPstageavg{ipatient}.(template), "NO_SCORE")
%                 LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'NO_SCORE');
%             end
%             
%             

%% average over patients
% SpikeDensity_avg_avg_avg = mean(vertcat(SpikeDensity_avg_avg{:}));
% 
% clear temp
% for ipatient = 1 : 8
%     i = 1;
%     
%     for ipart = 1 : 3
%         for iLFPtemplate = 1:6
%             
%             markername = sprintf('template%d', iLFPtemplate);
%             
%             if isfield(LFPavg_reref{ipatient}{ipart}, markername)
%                 temp{ipatient}{i} = LFPavg_reref{ipatient}{ipart}.(markername);
%                 i = i + 1;
%             end
%         end
%     end
%     cfg = [];
%     LFPavg_avg{ipatient} = ft_appenddata(cfg, temp{ipatient}{:});
%     cfg = [];
%     cfg.avgoverrpt = 'yes';
%     LFPavg_avg{ipatient} = ft_selectdata(cfg, LFPavg_avg{ipatient});
%     
% end
% 
% %%% ONLY AVERAGE MAXCHAN %%%%
% 
% cm = cbrewer('qual', 'Set2', 8);
% 
% % plot average over patients
% fig = figure('visible', true);
% set(fig, 'PaperPositionMode', 'auto');
% %     set(fig, 'position', get(0,'ScreenSize'));
% %     set(fig, 'position', [20 -60  papersize papersize*sqrt(2)]);
% set(fig, 'Renderer', 'Painters');
% 
% 
% % configure figure
% papersize       = 1000;
% nLFPtemplate    = 6;
% 
% hold
% set(gca, 'YScale', 'log')
% t = SpikeDensity{ipatient}{1}.psth.template3.time;
% 
% for ipatient = 1 : 8      
%     
% %         bar(temp_time{1}, SpikeDensity_avg_avg{ipatient}, 1, 'facecolor', cm(ipatient, :), 'edgecolor', 'none', 'facealpha', 0.8);
%         plot(t, SpikeDensity_avg_avg{ipatient}, 'color', cm(ipatient, :));
% end
%         yyaxis left
% 
%         xticks([]);
%         y = ylim;
%         yticks(floor(y(end)));
%         
%         % labels
%         if iLFPtemplate == 1 && ipatient ==8
%             xlabel('Time');
%             ylabel('Change vs. baseline');
%         end
%         
%         yyaxis right
%         hold;
%         % plot LFP
%         %         try
%         for iLFPtemplate = 1:6
%             markername = sprintf('template%d', iLFPtemplate);
%             if ~isfield(LFPavg_reref{ipatient}{ipart}, markername)
%                 continue
%             end 
%             
%             %             if ~((ipatient == 1 && iLFPtemplate == 4) || (ipatient == 4 && iLFPtemplate == 3) || (ipatient == 5 && iLFPtemplate == 4) || (ipatient == 5 && iLFPtemplate == 6))
%             %             plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}', '-k')
%             
%             maxchan = find(~cellfun(@isempty, strfind(LFPavg_reref{ipatient}{ipart}.(markername).label, config{ipatient}.align.zerochannel)), 1, 'last');
%             plot(LFPavg_reref{ipatient}{ipart}.(markername).time{1}, LFPavg_reref{ipatient}{ipart}.(markername).trial{1}(maxchan,:), '-k')
%             %             end
%             axis tight;
%             xlim(config{ipatient}.spike.toi.(markername));
%             y = ylim;
%             yticks(floor(y(end)));
%             xticks([]);
%             axis tight;
%             y = ylim;
%             yticks(sort([floor(y(1)), 0, floor(y(2))]));
%             
%         end
%         
%         xlim(config{ipatient}.spike.toi.(markername));
%         % labels
%         if iLFPtemplate == 1 && ipatient ==8
%             ylabel('Amplitude (mV)');
%         end
%         
%         % baseline line
%         if ipatient == 8
%             y = ylim;
%             y = y(1) + (y(2) - y(1)) * 0.05;
%             x = xlim;
%             width = x(2) - x(1);
%             line(config{ipatient}.stats.bl.(markername), [y, y], 'linewidth', 2, 'color', 'k');
%             text(config{ipatient}.stats.bl.(markername)(1) + width * 0.01, y, 'Baseline', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
%             
%             % time line
%             duration = 0.5;
%             t = sprintf('500 ms');
%             line([config{ipatient}.spike.toi.(markername)(end) - duration, config{ipatient}.spike.toi.(markername)(end)], [y, y], 'linewidth', 2, 'color', 'k');
%             text(config{ipatient}.spike.toi.(markername)(end) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
%         end
%         
%         %         fname = fullfile(config{ipatient}.imagesavedir, sprintf('firingrates_patient%d_part%d', ipatient, ipart));
%         % exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
%         %     exportgraphics(fig, strcat(fname, '.pdf'));
% end % ipatient
