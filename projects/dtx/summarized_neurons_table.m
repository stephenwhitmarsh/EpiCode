function neurons_table = summarized_neurons_table(cfg, ipart, force, SpikeStats_windowed, SpikeDensity_timelocked, SpikeWaveforms_stats)

fname = fullfile(cfg.datasavedir,'..',[cfg.prefix, 'p', num2str(ipart), '-neurons_table.xlsx']);

if exist(fname, 'file') && force == false
    fprintf('Reading %s\n', fname);
    neurons_table = readtable(fname);
    return
end

if isempty(SpikeDensity_timelocked)
    istimelocked = false;
else
    istimelocked = true;
end

neurons_table = table.empty;

i = 0;
for itemp = 1:size(SpikeStats_windowed{ipart}.window, 2)
    i = i+1;
    neurons_table.ratID{itemp}                      = cfg.prefix(1:end-1);
    neurons_table.clusterID{itemp}                  = SpikeStats_windowed{ipart}.window{itemp}.label;
    neurons_table.idx_orig{itemp}                   = i;
    neurons_table.group{itemp}                      = "";
    neurons_table.RPV{itemp}                        = SpikeStats_windowed{ipart}.window{itemp}.RPV * 100;
    neurons_table.freq{itemp}                       = nanmean(SpikeStats_windowed{ipart}.window{itemp}.trialfreq(SpikeStats_windowed{ipart}.window{itemp}.trialfreq < 100));
    neurons_table.CV2{itemp}                        = nanmean(SpikeStats_windowed{ipart}.window{itemp}.CV2_trial);
    neurons_table.burst_per_min{itemp}              = nanmean(SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum);
    neurons_table.amplitude{itemp}                  = SpikeWaveforms_stats{ipart}.window.amplitude.val(itemp);
    neurons_table.halfwidth{itemp}                  = SpikeWaveforms_stats{ipart}.window.halfwidth.val(itemp) * 1000;
    neurons_table.peaktrough{itemp}                 = SpikeWaveforms_stats{ipart}.window.peaktrough.val(itemp) * 1000;
    neurons_table.troughpeak{itemp}                 = SpikeWaveforms_stats{ipart}.window.troughpeak.val(itemp) * 1000;
    neurons_table.select_time{itemp}                = "";
    if istimelocked
        t_sel = SpikeDensity_timelocked{ipart}.sdf_lin.SlowWave.time > 0;
        neurons_table.sw_maxfreq{itemp}             = nanmax(SpikeDensity_timelocked{ipart}.sdf_lin.SlowWave.avg(itemp, :));
        neurons_table.sw_minfreq{itemp}             = nanmin(SpikeDensity_timelocked{ipart}.sdf_lin.SlowWave.avg(itemp, t_sel));
    end
end

%do not overwrite
if exist(fname, 'file')
    fname_orig = fname(1:end-5);
    i = 0;
    while exist(fname, 'file')
        i = i + 1;
        fname = sprintf('%s_%.3d.xlsx', fname_orig, i);
    end
end

writetable(neurons_table, fname);
