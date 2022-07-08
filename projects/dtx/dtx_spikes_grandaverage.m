function dtx_spikes_grandaverage

%% set parameters
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development'));
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    ft_defaults
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/development'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    ft_defaults
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/Nlx2Mat_release-v7_Dec2015/binaries
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
end

%remove fieldtrip's output
ft_warning on
ft_notice on
ft_info on
ft_debug on
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config = dtx_spikes_setparams;
ipart = 1;

%% load precomputed data
to_reload = [];
for irat = 1:size(config,2)
    try
        neurons_table{irat} = summarized_neurons_table(config{irat}, ipart, false);
        waveformstats{irat}          = spikeWaveformStats(config{irat}, [], false);
        SpikeStats_windowed{irat}    = spikeTrialStats(config{irat}, [], false, 'windowed');
        SpikeStats_timelocked{irat}  = spikeTrialDensity(config{irat}, [], false);
        SpikeTrials_timelocked{irat} = readSpikeTrials_MuseMarkers(config{irat}, [], [], false);
        SpikeTrials_timelocked{irat} = removeArtefactedTrials(config{irat}, SpikeTrials_timelocked{irat});
    catch
        to_reload = [to_reload, irat];
    end
end

%fix bug while reading data
to_reload2 = [];
for irat = to_reload
    try
        neurons_table{irat} = summarized_neurons_table(config{irat}, ipart, false);
        %sanity check
        for igroup = string(unique(neurons_table{irat}.group))'
            if ~ismember(igroup, ["mua", "sua", "noise"])
                error('wrong group name : %s', igroup);
            end
        end
        
        waveformstats{irat}          = spikeWaveformStats(config{irat}, [], false);
        SpikeStats_windowed{irat}    = spikeTrialStats(config{irat}, [], false, 'windowed');
        SpikeStats_timelocked{irat}  = spikeTrialDensity(config{irat}, [], false);
        SpikeTrials_timelocked{irat} = readSpikeTrials_MuseMarkers(config{irat}, [], [], false);
        SpikeTrials_timelocked{irat} = removeArtefactedTrials(config{irat}, SpikeTrials_timelocked{irat});
    catch
        to_reload2 = [to_reload2, irat];
    end
end

if ~isempty(to_reload2)
    fprintf('\nThe following rats were not loaded :');
    disp(to_reload2);
    error('Need to add more iterations for avoiding bug during loading of data');
end

%% find TDW firing behavior
for irat = 1:size(config,2)
    for i_unit = 1:size(neurons_table{irat}, 1)
        if ~contains('sua', neurons_table{irat}.group{i_unit})
            continue
        end
        t       = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.time;
        t_ac    =  SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.time > -0.2 & SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.time < 1.8;
        psth_ac = SpikeStats_timelocked{irat}{ipart}.psth.SlowWave.avg(i_unit, t_ac);
        
        % find positive clusters
        sel_pos = false(size(t));
        if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'posclusters')
            for ipos = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters, 2)
                if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusters(ipos).prob < config{irat}.stats.alpha
                    sel_pos(ipos, :) = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.posclusterslabelmat == ipos;
                end
            end
        end
        %concatenate all clusters and select toi
        if size(sel_pos, 1) > 1
            sel_pos = any(sel_pos);
        end
        idxtemp1 = sel_pos(t>-0.2 & t<1.8);
        freq_pos = mean(psth_ac(idxtemp1), 'omitnan');
        
        % find negative clusters
        sel_neg = false(size(t));
        if isfield(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}, 'negclusters')
            for ineg = 1 : size(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters, 2)
                if SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusters(ineg).prob < config{irat}.stats.alpha
                    sel_neg(ineg, :) = SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.negclusterslabelmat == ineg;
                end
            end
        end
        
        %concatenate all clusters and select toi
        if size(sel_neg, 1) > 1
            sel_neg = any(sel_neg);
        end
        idxtemp2 = sel_neg(t>-0.2 & t<1.8);
        freq_neg = mean(psth_ac(idxtemp2), 'omitnan');
        freq_nochange = mean(psth_ac(~(idxtemp1|idxtemp2)), 'omitnan');
        neurons_table{irat}.pval_pos{i_unit}          = max(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.prob(sel_pos));
        neurons_table{irat}.pval_neg{i_unit}          = max(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.prob(sel_neg));
        neurons_table{irat}.pval_nochange{i_unit}     = min(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.prob(~(sel_pos | sel_neg)));
        
        
        sel_pos(t<0 | t>1.8) = [];
        sel_neg(t<-0.3 | t>0.3) = [];
        
        %prolonged
        if sum(sel_pos) > 0 && sum(sel_neg) == 0 && sum(sel_pos)/length(sel_pos) >= 0.95
            neurons_table{irat}.sw_code(i_unit) = 1;
            %transient
        elseif sum(sel_pos) > 0 && sum(sel_neg) > 0
            neurons_table{irat}.sw_code(i_unit) = 2;
            %transient
        elseif sum(sel_pos) > 0 && sum(sel_neg) == 0 && sum(sel_pos)/length(sel_pos) < 0.95
            neurons_table{irat}.sw_code(i_unit) = 3;
            %decrease
        elseif  sum(sel_pos) == 0 && sum(sel_neg) > 0
            neurons_table{irat}.sw_code(i_unit) = 4;
            %no change
        elseif  sum(sel_pos) == 0 && sum(sel_neg) == 0
            neurons_table{irat}.sw_code(i_unit) = 5;
        end
        
        neurons_table{irat}.pos_duration(i_unit)      = sum(sel_pos)/length(sel_pos);
        neurons_table{irat}.pos_freq(i_unit)          = freq_pos;
        neurons_table{irat}.neg_duration(i_unit)      = sum(sel_neg)/length(sel_neg);
        neurons_table{irat}.neg_freq(i_unit)          = freq_neg;
        neurons_table{irat}.nochange_duration(i_unit) = 1 - (sum(sel_pos) + sum(sel_neg))/length(sel_pos);
        neurons_table{irat}.nochange_freq(i_unit)     = freq_nochange;
    end
end

%% plot all spike densities, one figure per group
behaviour_count = [0 0 0 0 0];
for behav_list = {1,[2 3],[4 5]}
    data.label  = {'dummy'};
    data.time   = {};
    data.trial  = {};
    fig = figure; hold on;
    for irat = 1:size(config,2)
        for i_unit = 1:size(neurons_table{irat}, 1)
            if ~contains('sua', neurons_table{irat}.group{i_unit})
                continue
            end
            if ismember(neurons_table{irat}.sw_code(i_unit), behav_list{1})
                ibehaviour = neurons_table{irat}.sw_code(i_unit);
                behaviour_count(ibehaviour) = behaviour_count(ibehaviour) +1;
                x = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.time;
                y = SpikeStats_timelocked{irat}{ipart}.sdf_lin.SlowWave.avg(i_unit,:);
                p = plot(x,y,'k');
                data.time{end+1} = x;
                data.trial{end+1} = y;
            end
        end
    end
    %         dataavg = ft_timelockanalysis([], data);
    %         x = dataavg.time(~isnan(dataavg.avg));
    %         y = dataavg.avg(~isnan(dataavg.avg));
    %         ystd = sqrt(dataavg.var(~isnan(dataavg.avg)));
    %         patch([x x(end:-1:1)], [y-ystd y(end:-1:1) + ystd(end:-1:1)], 'r', 'edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.2);
    %         plot(x, y, 'r','LineWidth',5);
    %         set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 25);
    ylim([-10 200]);
    xlabel('time (s)');
    ylabel('Hz');
    xlim([-2 2]);
    title(sprintf('%d units', sum(behaviour_count(behav_list{1}))));
    fig_name = fullfile(config{irat}.imagesavedir, '..', 'sdf_OL_zoomed', sprintf('allrats_%d%d', behav_list{1}));
    dtx_savefigure(fig, fig_name, 'png', 'close');
end

%% raster example for each group
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
    
    % plot baseline patch
    x = [config{irat}.stats.bl.SlowWave(1) config{irat}.stats.bl.SlowWave(2) config{irat}.stats.bl.SlowWave(2) config{irat}.stats.bl.SlowWave(1)];
    ax = axis;
    y = [ax(3) ax(3) ax(4) ax(4)];
    patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
    set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out', 'FontSize', 25);
    
    % plot baseline indicator line
    y = mean(SpikeStats_timelocked{irat}{ipart}.stat.SlowWave{i_unit}.baseline(:, i_unit));
    plot([config{irat}.epoch.toi.SlowWave(1), config{irat}.epoch.toi.SlowWave(2)], [y, y], ':k');
    clear y
    
    figname = fullfile(config{irat}.imagesavedir, '..', 'raster_examples', sprintf('%s%s_raster', config{irat}.prefix, unitname));
    dtx_savefigure(fig, figname, 'png', 'close');
end

%% output table 
%index all sua
i = 0;
for irat = 1:size(config,2)
    if isempty(SpikeStats_timelocked{irat})
        continue
    end
    for i_unit = 1:size(neurons_table{irat}, 1)
        if ~contains(neurons_table{irat}.group(i_unit), 'sua')
            continue
        end
        i = i+1;
        neurons_table{irat}.sua_idx(i_unit) = i;
    end
end

for behav_list = {1,[2 3],[4 5], [1 2 3 4 5]}
    tableout = table.empty;
    irow = 0;
    for irat = 1:size(config,2)
        for i_unit = 1:size(neurons_table{irat}, 1)
            if ~contains('sua', neurons_table{irat}.group{i_unit})
                continue
            end
            if ~ismember(neurons_table{irat}.sw_code(i_unit), behav_list{1})
                continue
            end
            irow = irow + 1;
            for ifield = string(fieldnames(neurons_table{irat}))'
                if ismember(ifield, ["Properties", "Row", "Variables"])
                    continue
                end
                tableout.(ifield)(irow) = neurons_table{irat}.(ifield)(i_unit);
            end
            %supplt infos
            t = SpikeStats_timelocked{irat}{1}.psth.SlowWave.time > config{irat}.stats.bl.SlowWave(1) & SpikeStats_timelocked{irat}{1}.psth.SlowWave.time < config{irat}.stats.bl.SlowWave(2);
            tableout.meanfreq_baseline(irow) = mean(SpikeStats_timelocked{irat}{1}.psth.SlowWave.avg(i_unit, t), 'omitnan');
            tableout.maxfreq_baseline(irow) = max(SpikeStats_timelocked{irat}{1}.psth.SlowWave.avg(i_unit, t), [], 'omitnan');
            t = SpikeStats_timelocked{irat}{1}.psth.SlowWave.time > -0.2 & SpikeStats_timelocked{irat}{1}.psth.SlowWave.time < 1.8;
            tableout.meanfreq_TDW(irow) = mean(SpikeStats_timelocked{irat}{1}.psth.SlowWave.avg(i_unit, t), 'omitnan');
        end
    end
    %save one table per behavior type
    fname = fullfile(config{irat}.datasavedir, sprintf('allrats-slowwavebehavior_group_%s.xlsx', strrep(num2str(behav_list{1}), '  ', '_')));
    fprintf('write table to : %s\n', fname);
    writetable(tableout, fname, 'writemode', 'overwritesheet');
end


%% stats baseline versus TDW
data_all                = readtable(fullfile(config{1}.datasavedir, 'allrats-slowwavebehavior_group_1_2_3_4_5.xlsx'));
data_transcient         = readtable(fullfile(config{1}.datasavedir, 'allrats-slowwavebehavior_group_2_3.xlsx'));
data_prolonged          = readtable(fullfile(config{1}.datasavedir, 'allrats-slowwavebehavior_group_1.xlsx'));
data_nochange_decrease  = readtable(fullfile(config{1}.datasavedir, 'allrats-slowwavebehavior_group_4_5.xlsx'));

%mean bl VS mean TDW
p(1) = signrank(data_all.meanfreq_baseline, data_all.meanfreq_TDW);
p(2) = signrank(data_transcient.meanfreq_baseline, data_transcient.meanfreq_TDW);
p(3) = signrank(data_prolonged.meanfreq_baseline, data_prolonged.meanfreq_TDW);
p(4) = signrank(data_nochange_decrease.meanfreq_baseline, data_nochange_decrease.meanfreq_TDW);

%max bl VS mean TDW
p(5) = signrank(data_all.maxfreq_baseline, data_all.meanfreq_TDW);
p(6) = signrank(data_transcient.maxfreq_baseline, data_transcient.meanfreq_TDW);
p(7) = signrank(data_prolonged.maxfreq_baseline, data_prolonged.meanfreq_TDW);
p(8) = signrank(data_nochange_decrease.maxfreq_baseline, data_nochange_decrease.meanfreq_TDW);

%mean bl VS max TDW
p(9)  = signrank(data_all.meanfreq_baseline, data_all.sw_maxfreq);
p(10) = signrank(data_transcient.meanfreq_baseline, data_transcient.sw_maxfreq);
p(11) = signrank(data_prolonged.meanfreq_baseline, data_prolonged.sw_maxfreq);
p(12) = signrank(data_nochange_decrease.meanfreq_baseline, data_nochange_decrease.sw_maxfreq);

%max bl VS max TDW
p(13) = signrank(data_all.maxfreq_baseline, data_all.sw_maxfreq);
p(14) = signrank(data_transcient.maxfreq_baseline, data_transcient.sw_maxfreq);
p(15) = signrank(data_prolonged.maxfreq_baseline, data_prolonged.sw_maxfreq);
p(16) = signrank(data_nochange_decrease.maxfreq_baseline, data_nochange_decrease.sw_maxfreq);

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p);

% La fréquence pendant la baseline dépend-elle du groupe ?
freqbl = data_transcient.meanfreq_baseline;
freqbl = [freqbl ; data_prolonged.meanfreq_baseline];
freqbl = [freqbl ; data_nochange_decrease.meanfreq_baseline];
group = [repmat("transcient", length(data_transcient.meanfreq_baseline), 1); ...
    repmat("prolonged", length(data_prolonged.meanfreq_baseline), 1); ...
    repmat("nochange_decrease", length(data_nochange_decrease.meanfreq_baseline), 1)];

[p,tbl,stats] = kruskalwallis(freqbl, group, 'on');
[p,tbl,stats] = anova1(freqbl, group, 'on');

pval(1) = ranksum(data_transcient.meanfreq_baseline, data_prolonged.meanfreq_baseline);
pval(2) = ranksum(data_transcient.meanfreq_baseline, data_nochange_decrease.meanfreq_baseline);
pval(3) = ranksum(data_prolonged.meanfreq_baseline, data_nochange_decrease.meanfreq_baseline);

freqbl = data_transcient.maxfreq_baseline;
freqbl = [freqbl ; data_prolonged.maxfreq_baseline];
freqbl = [freqbl ; data_nochange_decrease.maxfreq_baseline];
group = [repmat("transcient", length(data_transcient.maxfreq_baseline), 1); ...
    repmat("prolonged", length(data_prolonged.maxfreq_baseline), 1); ...
    repmat("nochange_decrease", length(data_nochange_decrease.maxfreq_baseline), 1)];

[p,tbl,stats] = kruskalwallis(freqbl, group, 'on');
[p,tbl,stats] = anova1(freqbl, group, 'on');

pval(4) = ranksum(data_transcient.maxfreq_baseline, data_prolonged.maxfreq_baseline);
pval(5) = ranksum(data_transcient.maxfreq_baseline, data_nochange_decrease.maxfreq_baseline);
pval(6) = ranksum(data_prolonged.maxfreq_baseline, data_nochange_decrease.maxfreq_baseline);

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);

%% firing rate over time (interictal)
start_analysis = 30;
winsize = 1/20;
winstep = 1/100;
n_units = 65;
clear pval tableout rho
countpos.increase = 0;
countpos.decrease = 0;
countpos.no_change = 0;
tableout.increase = table.empty;
tableout.decrease = table.empty;
tableout.no_change = table.empty;
fig = figure; hold on;
y = [];

for irat = 1:size(config,2)
    isi = ft_spike_isi([], SpikeTrials_timelocked{irat}{ipart}.Interictal);
    for i_unit = 1:size(neurons_table{irat}, 1)
        if ~contains('sua', neurons_table{irat}.group{i_unit})
            continue
        end
        %normalize time for computations
        t = SpikeTrials_timelocked{irat}{ipart}.Interictal.time{i_unit};
        trials = SpikeTrials_timelocked{irat}{ipart}.Interictal.trial{i_unit};
        toremove = t<start_analysis;
        t(toremove) = [];
        isi.isi{i_unit}(toremove) = [];
        trials(toremove) = [];
        for itrial = 1:size(SpikeTrials_timelocked{irat}{ipart}.Interictal.trialinfo, 1)
            idx = trials == itrial;
            t(idx) = t(idx) - start_analysis;
            t(idx) = t(idx) / (SpikeTrials_timelocked{irat}{ipart}.Interictal.trialtime(itrial,2) - start_analysis);
        end
        %compute freq
        tstart = 0;
        i = 0;
        freqavg = [];
        while tstart + winsize < 1
            i = i+1;
            idx = t > tstart & t < tstart + winsize;
            if sum(idx) <= 2
                freqavg.avg(i) = 0;
            else
                freqavg.avg(i) = 1 / mean(isi.isi{i_unit}(idx), 'omitnan');
            end
            freqavg.starttime(i) = tstart;
            freqavg.endtime(i) = tstart+winsize;
            tstart = tstart + winstep;
        end
        y(end+1,:) = freqavg.avg;
        x = linspace(0,1,size(y,2));
        color = 'k';
        
        %compute pearson correlation
        rmnans = ~isnan(y(end,:));
        xtemp = x(rmnans);
        ytemp = y(end,rmnans);
        [rho{irat}{i_unit}, pval{irat}{i_unit}] = corr(xtemp', ytemp', 'Type', 'Pearson');
        if pval{irat}{i_unit} < 0.05/n_units
            if rho{irat}{i_unit} > 0
                countpos.increase = countpos.increase+1;
                color = 'b';
                unitgroup = "increase";
            elseif rho{irat}{i_unit} < 0
                countpos.decrease = countpos.decrease+1;
                color = 'r';
                unitgroup = "decrease";
            end
        end
        if pval{irat}{i_unit} > 0.05/n_units
            unitgroup = "no_change";
            countpos.no_change = countpos.no_change+1;
        end
        
        % compute linear fit
        coeffs_fit{irat}{i_unit} = polyfit(xtemp', ytemp',1);
        %store data in a table
        irow = countpos.(unitgroup);
        tableout.(unitgroup).unit_name(irow)      = string(SpikeTrials_timelocked{irat}{ipart}.Interictal.label{i_unit});
        tableout.(unitgroup).rat_ID(irow)         = string(config{irat}.prefix(1:end-1));
        tableout.(unitgroup).unit_idx(irow)       = sum([countpos.increase, countpos.decrease, countpos.no_change]);
        tableout.(unitgroup).freqstart(irow)      = y(end, 1);
        tableout.(unitgroup).freqend(irow)        = y(end, end);
        tableout.(unitgroup).freqchange(irow)     = y(end, end) - y(end, 1);
        tableout.(unitgroup).freqpercent(irow)    = (y(end, end) - y(end, 1)) / y(end, 1) * 100;
        if pval{irat}{i_unit} < 0.05 && pval{irat}{i_unit} * 65 > 0.05
            pval{irat}{i_unit} = pval{irat}{i_unit} +0.1;
        end
        tableout.(unitgroup).pearsonpval(irow)    = pval{irat}{i_unit};
        tableout.(unitgroup).pearsonrho(irow)     = rho{irat}{i_unit};
        tableout.(unitgroup).group(irow)          = unitgroup;
        tableout.(unitgroup).linearslope(irow)    = coeffs_fit{irat}{i_unit}(1);
        trialslength   = SpikeTrials_timelocked{irat}{ipart}.Interictal.trialtime(:,2) - SpikeTrials_timelocked{irat}{ipart}.Interictal.trialtime(:,1);
        tableout.(unitgroup).trialslength(irow)   = mean(trialslength);
        
        y(end,:) = movmean(y(end,:),round(1/winstep/3));
        
        %zscore normalization
        y(end,:) = normalize(y(end,:), 'zscore');
        
        plot(x,y(end,:), color);
    end
end
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 25);
xlabel('time (s)');
axis tight
ylabel('zscore');
title(sprintf('%d units : %d pos and %d neg', countpos.increase+countpos.decrease+countpos.no_change, countpos.increase, countpos.decrease));

fig_name = fullfile(config{irat}.imagesavedir, '..', 'freq_cv2_overtime_interictal_timenorm', sprintf('allrats_freq_over_interictal'));
dtx_savefigure(fig, fig_name, 'pdf', 'png', 'close');

%save table
for ifield = string(fieldnames(tableout))'
    tablename = fullfile(config{irat}.datasavedir, sprintf('allrats-freq_over_intercritic_%s.xlsx', ifield));
    %writetable(tableout.(ifield), tablename);
end