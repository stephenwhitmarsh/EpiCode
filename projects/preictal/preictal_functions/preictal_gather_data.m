function data = preictal_gather_data(cfg, force, MuseStruct, LFP, FFT, waveformstats, spikestats, spiketrials)

fname = fullfile(cfg.datasavedir, 'data_gathered', sprintf('%sall_data.mat', cfg.prefix));

%% load precomputed data if required 
if nargin <= 2 && exist(fname, 'file')
    fprintf('Reading %s\n', fname);
    load(fname, 'data');
    return
elseif nargin <= 2
    error('Cannot find %s \nNot enough arguments to create the data structure', fname);
end

if force == false && exist(fname, 'file')
    fprintf('Reading %s\n', fname);
    load(fname, 'data');
    return
end

%% create the data structure

data.pat_ID             = cfg.prefix(1:end-1);
data.units_label        = waveformstats{1}.label;
data.units_group        = waveformstats{1}.cluster_group;
for itemp = 1:length(data.units_label)
    maxchan             = spiketrials{1}.window_60.template_maxchan(itemp) + 1; %+1 because starts at zero
    data.units_channels{itemp} = cfg.circus.channel{maxchan};
end

data.info1               = "Artefacted trials are removed and not present in this structure.";
data.info2               = "One trial = one window. The times of the windows are stored for each structure in trialinfo.starttime and trialinfo.endtime.";
if strcmp(cfg.seizure_index, 'last')
    data.Crise_Start = MuseStruct{1}{2}.markers.CriseStart.clock(end);
else
    data.Crise_Start = MuseStruct{1}{2}.markers.CriseStart.clock(cfg.seizure_index);
end
data.spikewaveforms      = waveformstats{1};
data.win_60s.LFP        = LFP{1}.window_60;
data.win_60s.FFT        = FFT{1}.window_60;
data.win_60s.spiketrials = spiketrials{1}.window_60;
data.win_60s.spiketrials_info = "Details about the spike trial structure organisation : http://old.fieldtriptoolbox.org/reference/ft_datatype_spike";
data.win_60s.spikestats = spikestats{1}.window_60;
data.win_60s.spikestats_info  = "One cellule corresponds to one unit";
data.win_3s.LFP         = LFP{1}.window_3;
data.win_3s.FFT         = FFT{1}.window_3;
data.win_3s.spiketrials = spiketrials{1}.window_3;
data.win_3s.spiketrials_info = "Details about the spike trial structure organisation : http://old.fieldtriptoolbox.org/reference/ft_datatype_spike";
data.win_3s.spikestats  = spikestats{1}.window_3;
data.win_3s.spikestats_info  = "One cellule corresponds to one unit";

%% save to disk
isdir_or_mkdir(fileparts(fname));
fprintf('Writing data to %s\n', fname);
save(fname, 'data');
