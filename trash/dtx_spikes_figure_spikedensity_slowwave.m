
%% set parameters
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    ft_defaults
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/development'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    ft_defaults
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/Nlx2Mat_release-v7_Dec2015/binaries
end

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info on
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config = dtx_spikes_setparams;
ipart = 1;


%% load precomputed data
for irat = 1:size(config,2)
    
    neurons_table{irat} = summarized_neurons_table(config{irat}, ipart, false);
    %sanity check
    for igroup = string(unique(neurons_table{irat}.group))'
        if ~ismember(igroup, ["mua", "sua", "noise"])
            error('wrong group name : %s');
        end
    end
        
        waveformstats{irat}          = spikeWaveformStats(config{irat}, [], false);
        SpikeStats_timelocked{irat}  = spikeTrialDensity(config{irat}, [], false);
        SpikeTrials_timelocked{irat} = readSpikeTrials_MuseMarkers(config{irat}, [], [], false);
        SpikeTrials_timelocked{irat} = removeArtefactedTrials(config{irat}, SpikeTrials_timelocked{irat});
%         LFP{irat}                    = dtx_readLFP(config{irat}, [], false);
      
end


%% trouver le groupe de comportement pendant l'onde lente
for irat = 1:size(config,2)
    if isempty(SpikeStats_timelocked{irat})
        continue
    end
    for i_unit = 1:size(neurons_table{irat}, 1)
        if contains('noise', neurons_table{irat}.group{i_unit})
            continue
        end
        t = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time;
        
        % find positive clusters
        %lag = size(SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.time, 2) - size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time, 2);
        sel_pos = false(size(t));
        if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'posclusters')
            for ipos = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters, 2)
                if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters(ipos).prob < config{irat}.stats.alpha
                    %sel = find(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusterslabelmat == ipos);
                    sel_pos(ipos, :) = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusterslabelmat == ipos;
                    %bar(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time(sel), SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.avg(i_unit, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
                end
            end
        end
        %concatenate all clusters
        if size(sel_pos, 1) > 1
            sel_pos = any(sel_pos);
        end
        %keep only time of interest [0 1.8]
        sel_pos(t<0 | t>1.8) = [];
        
        % find negative clusters
        sel_neg = false(size(t));
        if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'negclusters')
            for ineg = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters, 2)
                if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters(ineg).prob < config{irat}.stats.alpha
                    sel_neg(ineg, :) = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusterslabelmat == ineg;
                    %bar(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time(sel), SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.avg(i_unit, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
                end
            end
        end
        %concatenate all clusters
        if size(sel_neg, 1) > 1
            sel_neg = any(sel_neg);
        end
        %keep only time of interest [0 1.8]
        sel_neg(t<0 | t>1.8)  = [];
        
        %groupe 1 : décharge augmentée et continue
        if sum(sel_pos) > 0 && sum(sel_neg) == 0 && sum(sel_pos)/length(sel_pos) > 0.5
            neurons_table{irat}.sw_code(i_unit) = 1;
            %groupe 2 : décharge augmentée puis diminuée
        elseif sum(sel_pos) > 0 && sum(sel_neg) > 0
            neurons_table{irat}.sw_code(i_unit) = 2;
            %groupe 3 : décharge augmentée puis retour à baseline
        elseif sum(sel_pos) > 0 && sum(sel_neg) == 0 && sum(sel_pos)/length(sel_pos) < 0.5
            neurons_table{irat}.sw_code(i_unit) = 3;
            %groupe 4 : décharge diminuée sans augmentation
        elseif  sum(sel_pos) == 0 && sum(sel_neg) > 0
            neurons_table{irat}.sw_code(i_unit) = 4;
            %groupe 5 : pas d'effet
        elseif  sum(sel_pos) == 0 && sum(sel_neg) == 0
            neurons_table{irat}.sw_code(i_unit) = 5;
        end
    end
end

%% plot tous les neurones : sdf line SlowWave de chaque neurone en fonction du groupe (une figure par groupe)
% sua seuls, sua + mua

for igroup = ["sua"]%, "sua_mua"]    
    behaviour_count.(igroup) = [0 0 0 0 0];
    behaviour_check.(igroup) = cell(1,5);
    for behav_list = {1,[2 3],[4 5]}
        data.label  = {'dummy'};
        data.time   = {};
        data.trial  = {};
        fig = figure; hold on;
        for irat = 1:size(config,2)
            if isempty(SpikeStats_timelocked{irat})
                continue
            end
            for i_unit = 1:size(neurons_table{irat}, 1)
                if ~contains(igroup, neurons_table{irat}.group{i_unit})
                    continue
                end
                if ismember(neurons_table{irat}.sw_code(i_unit), behav_list{1})
                    ibehaviour = neurons_table{irat}.sw_code(i_unit);
                    behaviour_count.(igroup)(ibehaviour) = behaviour_count.(igroup)(ibehaviour) +1;
                    behaviour_check.(igroup){ibehaviour}(end+1) = irat + i_unit/100;
                    x = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.time;
                    y = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.avg(i_unit,:);
                    
                    if max(y) > 150
                        error('stop here')
                    end
                    
                    p = plot(x,y,'k');
                    %p.Color(4) = 0.3;
                    %store to average : 
                    data.time{end+1} = x;
                    data.trial{end+1} = y;
                end
            end
        end
        dataavg = ft_timelockanalysis([], data);
        x = dataavg.time(~isnan(dataavg.avg));
        y = dataavg.avg(~isnan(dataavg.avg));
        ystd = sqrt(dataavg.var(~isnan(dataavg.avg)));
        patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], 'r', 'edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.2);
        plot(x, y, 'r','LineWidth',5);
        set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 25);
        ylim([-10 200]);
        xlabel('time (s)');
        ylabel('Hz');
        title(sprintf('%d units', sum(behaviour_count.(igroup)(behav_list{1}))));
        fig_name = fullfile(config{irat}.imagesavedir, '..', 'sdf_OL', sprintf('allrats_%s_%d%d', igroup, behav_list{1}));
        dtx_savefigure(fig, fig_name, 'pdf', 'png', 'close');
    end
end

%% plot average : sdf line SlowWave de chaque neurone en fonction du groupe (une figure par groupe)
% sua seuls, sua + mua
% 1 : decharge continue
% 2 : decharge provisoire
% 3 : arrêt sans augmentation
for igroup = ["sua", "sua_mua"]    
    behaviour_count.(igroup) = [0 0 0 0 0 0];
    behaviour_check.(igroup) = cell(1,6);
    for behav_list = {1,2,3,[2 3],4,5, [4 5]}
        data.label  = {'dummy'};
        data.time   = {};
        data.trial  = {};
        fig = figure; hold on;
        for irat = 1:size(config,2)
            if isempty(SpikeStats_timelocked{irat})
                continue
            end
            for i_unit = 1:size(neurons_table{irat}, 1)
                if ~contains(igroup, neurons_table{irat}.group{i_unit})
                    continue
                end
                if ismember(neurons_table{irat}.sw_code(i_unit), behav_list{1})
                    ibehaviour = neurons_table{irat}.sw_code(i_unit);
                    behaviour_count.(igroup)(ibehaviour) = behaviour_count.(igroup)(ibehaviour) +1;
                    behaviour_check.(igroup){ibehaviour}(end+1) = irat + i_unit/100;
                    x = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.time;
                    y = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.avg(i_unit,:);
                    
                    %p = plot(x,y,'k');
                    %p.Color(4) = 0.3;
                    %store to average : 
                    data.time{end+1} = x;
                    data.trial{end+1} = y;
                end
            end
        end
        dataavg = ft_timelockanalysis([], data);
        x = dataavg.time(~isnan(dataavg.avg));
        y = dataavg.avg(~isnan(dataavg.avg));
        ystd = sqrt(dataavg.var(~isnan(dataavg.avg)));
        patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], 'k', 'edgecolor', 'none', 'facecolor', 'k', 'facealpha', 0.2);
        plot(x, y, 'k','LineWidth',2);
        set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 25);
        xlabel('time (s)');
        ylabel('Hz');
        ylim([-5 100]);
        title(sprintf('%d units', sum(behaviour_count.(igroup)(behav_list{1}))));
        fig_name = fullfile(config{irat}.imagesavedir, '..', 'sdf_OL', sprintf('allrats_average_%s_%d%d', igroup, behav_list{1}));
        dtx_savefigure(fig, fig_name, 'pdf', 'png', 'close');
    end
end

%% raster example for each group
% Cluster 6 E07 : groupe 1 : décharge continue
% Cluster 7 E14 : groupe 2 : décharge transitoire
% Cluster 0 E10 : groupe 3 : pas d'augmentation significative de décharge

irat = 1;
for unitname = ["cluster_6_E07", "cluster_7_E14", "cluster_0_E10"]
    i_unit = find(strcmp(unitname, SpikeTrials_timelocked{irat}{ipart}.SlowWave.label));
    fig = figure;subplot(2,1,1);
    cfgtemp                 = [];
    cfgtemp.spikechannel    = i_unit;
    cfgtemp.trialborders    = 'yes';
    cfgtemp.linewidth       = 4;
    ft_spike_plot_raster(cfgtemp, SpikeTrials_timelocked{irat}{ipart}.SlowWave);
    set(gca, 'TickDir', 'out');
    title('Rasterplot');
    ylabel('Trial nr');
    set(gca, 'TickDir', 'out', 'FontSize', 25);
    xlabel('Time (s)');
    
    
    subplot(2,1,2); hold on;
    title('Spike Density (Hz)');
    bar(SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.time, squeeze(nanmean(SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.trial(:, i_unit, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
    set(gca, 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
    %ylabel('Firing rate (Hz)');
    xlabel('Time (s)');
    
    % plot positive clusters
    lag = size(SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.time, 2) - size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time, 2);
    
    if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'posclusters')
        for ipos = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters, 2)
            if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters(ipos).prob < config{irat}.stats.alpha
                sel = find(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusterslabelmat == ipos);
                bar(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time(sel), SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.avg(i_unit, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
            end
        end
    end
    
    % plot negative clusters
    if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'negclusters')
        for ineg = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters, 2)
            if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters(ineg).prob < config{irat}.stats.alpha
                sel = find(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusterslabelmat == ineg);
                bar(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time(sel), SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.avg(i_unit, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
            end
        end
    end
    
    % plot smoothed
    plot(SpikeStats_timelocked{irat}{ipart}.sdf_bar.SlowWave.time, squeeze(nanmean(SpikeStats_timelocked{irat}{ipart}.sdf_bar.SlowWave.trial(:, i_unit, :), 1))', 'Color', [0 0 0], 'LineWidth', 2);
    
%     % plot baseline patch
%     x = [config{irat}.stats.bl.SlowWave(1) config{irat}.stats.bl.SlowWave(2) config{irat}.stats.bl.SlowWave(2) config{irat}.stats.bl.SlowWave(1)];
%     ax = axis;
%     y = [ax(3) ax(3) ax(4) ax(4)];
%     patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
%     set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out', 'FontSize', 25);
%     
%     % plot baseline indicator line
%     y = mean(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.baseline(:, i_unit));
%     plot([config{irat}.epoch.toi.SlowWave(1), config{irat}.epoch.toi.SlowWave(2)], [y, y], ':k');
%     clear y
    
%     ylim([-10 200]);

    figname = fullfile(config{irat}.imagesavedir, '..', 'raster_examples', sprintf('%s%s_raster', config{irat}.prefix, unitname));
    dtx_savefigure(fig, figname, 'png', 'pdf', 'close');
end


