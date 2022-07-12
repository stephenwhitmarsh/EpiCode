function neurons_table = PET_neurons_table_timelocked(cfg, ipart, statsPSTH, spiketrials, markername, force)

fname = fullfile(cfg.tablesavedir,[cfg.prefix, 'p', num2str(ipart), '-neurons_table_timelocked_', char(markername), '.xlsx']);

if exist(fname, 'file') && force == false
    fprintf('Reading %s\n', fname);
    neurons_table = readtable(fname);
    return
end

neurons_table = table.empty;

for itemp = 1:size(statsPSTH{ipart}.psth.(markername).label, 2)
    t       = statsPSTH{ipart}.stat.(markername){itemp}.time;
    t_bl    = statsPSTH{ipart}.psth.(markername).time >= cfg.stats.bl.(markername)(1) & statsPSTH{ipart}.psth.(markername).time <= cfg.stats.bl.(markername)(2);
    t_ac    = statsPSTH{ipart}.psth.(markername).time > cfg.stats.bl.(markername)(2);
    psth_ac = statsPSTH{ipart}.psth.(markername).avg(itemp, t_ac);
    psth_bl = statsPSTH{ipart}.psth.(markername).avg(itemp, t_bl);
    
    % find positive clusters
    %lag = size(spikestats{ipart}.psth.(markername).time, 2) - size(spikestats{ipart}.stat.(markername){i_unit}.time, 2);
    sel_pos = false(size(t));
    if isfield(statsPSTH{ipart}.stat.(markername){itemp}, 'posclusters')
        for ipos = 1 : size(statsPSTH{ipart}.stat.(markername){itemp}.posclusters, 2)
            if statsPSTH{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                %sel = find(spikestats{ipart}.stat.(markername){i_unit}.posclusterslabelmat == ipos);
                sel_pos(ipos, :) = statsPSTH{ipart}.stat.(markername){itemp}.posclusterslabelmat == ipos;
                %bar(spikestats{ipart}.stat.(markername){i_unit}.time(sel), spikestats{ipart}.psth.(markername).avg(i_unit, sel+lag), 1, 'facecolor', 'r', 'edgecolor', 'none');
            end
        end
    end
    %concatenate all clusters
    if size(sel_pos, 1) > 1
        sel_pos = any(sel_pos);
    end
    data_pos = psth_ac(sel_pos);
    
    % find negative clusters
    sel_neg = false(size(t));
    if isfield(statsPSTH{ipart}.stat.(markername){itemp}, 'negclusters')
        for ineg = 1 : size(statsPSTH{ipart}.stat.(markername){itemp}.negclusters, 2)
            if statsPSTH{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                sel_neg(ineg, :) = statsPSTH{ipart}.stat.(markername){itemp}.negclusterslabelmat == ineg;
                %bar(spikestats{ipart}.stat.(markername){i_unit}.time(sel), spikestats{ipart}.psth.(markername).avg(i_unit, sel+lag), 1, 'facecolor', 'b', 'edgecolor', 'none');
            end
        end
    end
    %concatenate all clusters
    if size(sel_neg, 1) > 1
        sel_neg = any(sel_neg);
    end
    %keep only time of interest [0 1.8]
    data_neg = psth_ac(sel_neg);
    data_nochange = psth_ac(~(sel_pos|sel_neg));
    
    steps = zeros(1, length(sel_pos));
    steps(sel_pos) = 1;
    steps(sel_neg) = 2;
    temp = find(diff(steps)); %detect changes
    changes = steps([1, temp+1]);
    %convert to sting
    stepstrings = [];
    for i = changes
        switch i
            case 0
                %stepstrings = [stepstrings, "bl."];
            case 1
                stepstrings = [stepstrings, "increase"];
            case 2
                stepstrings = [stepstrings, "decrease"];
        end
    end
    %         if stepstrings(1) == " bl." && length(stepstrings) > 1
    %             stepstrings = stepstrings(2:end);
    %         end
    stepstrings_out = [];
    if isempty(stepstrings)
        stepstrings_out = 'no_significative_change';
    else
        for i = 1 : size(stepstrings,2)
            if i == 1
                stepstrings_out = stepstrings(i);
            else
                if stepstrings(i) == stepstrings(i-1)
                    continue
                end
                stepstrings_out = sprintf('%s_then_%s', stepstrings_out, stepstrings(i));
            end
        end
    end
            
    %check that the different structure have the same order :
    assert(strcmp(statsPSTH{ipart}.psth.(markername).label{itemp}, spiketrials{ipart}.(markername).label{itemp}));
    
    neurons_table.patient_ID{itemp}                 = cfg.prefix(1:end-1);
    neurons_table.clusterID{itemp}                  = statsPSTH{ipart}.psth.(markername).label{itemp};
    neurons_table.cluster_group{itemp}              = strrep(spiketrials{ipart}.(markername).cluster_group{itemp}, ' ', '');
    maxchan       = spiketrials{ipart}.(markername).template_maxchan(itemp) + 1; %+1 because starts at zero
    neurons_table.channel{itemp}                    = cfg.circus.channel{maxchan};
    if ~isfield(cfg.circus, 'channelname')
        a = split(neurons_table.clusterID{itemp}, '_');
        neurons_table.electrode_bundle{itemp}       = a{3};
    end
    
    neurons_table.n_trials{itemp}               = size(spiketrials{ipart}.(markername).trialinfo, 1);
    neurons_table.n_PAs{itemp}               = length(spiketrials{ipart}.(markername).time{itemp});
    neurons_table.unit_behavior{itemp}          = stepstrings_out;
    neurons_table.freq_baseline_mean{itemp}     = mean(psth_bl);
    neurons_table.freq_baseline_std{itemp}      = std(psth_bl);
    neurons_table.freq_baseline_max{itemp}      = max(psth_bl);
    neurons_table.freq_baseline_min{itemp}      = min(psth_bl);
    neurons_table.freq_active_mean{itemp}       = mean(psth_ac);
    neurons_table.freq_active_std{itemp}        = std(psth_ac);
    neurons_table.freq_active_max{itemp}        = max(psth_ac);
    neurons_table.freq_active_min{itemp}        = min(psth_ac);
    neurons_table.pos_freq_mean(itemp)          = mean(data_pos);
    neurons_table.pos_freq_std(itemp)           = std(data_pos);
    neurons_table.neg_freq_mean(itemp)          = mean(data_neg);
    neurons_table.neg_freq_freq_std(itemp)      = std(data_neg);
    neurons_table.nochange_freq_mean(itemp)     = mean(data_nochange);
    neurons_table.nochange_freq_freq_std(itemp) = std(data_nochange);
    
end

delete(fname)
writetable(neurons_table, fname, 'WriteMode', 'overwritesheet');