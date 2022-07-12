function preictal_project(ielec)

restoredefaultpath
if ispc
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\shared'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\external'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\templates'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\development\modified_fieldtrip_functions'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))
    addpath \\l2export\iss02.charpier\analyses\vn_preictal\scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/development/modified_fieldtrip_functions'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
    addpath /network/lustre/iss02/charpier/analyses/vn_preictal/scripts/fieldtrip-20200607
end
ft_defaults

config = preictal_setparams;
ipart  = 1;
redo_analysis = true;

%% analysis
% read muse markers
MuseStruct = readMuseMarkers(config{ielec}, redo_analysis);

% remove post ictal from the whole analysis,
% according to config (some 'patients' will have
% shorter postictal kept because of noise, see setparams)
MuseStruct                  = addMuseBAD(config{ielec},MuseStruct);

% add (sliding) timewindow
[config{ielec}, MuseStruct] = addSlidingWindows(config{ielec}, MuseStruct);

% %FIXME TOREMOVE : pour debuguer quand il n'y pas pas le bon nom d'electrode dans le setparams
% for ichan = 1:size(config{ielec}.LFP.channel, 2)
%     if ~ismember(config{ielec}.LFP.channel{ichan}(end-1), '_')
%         config{ielec}.LFP.channel{ichan} = [config{ielec}.LFP.channel{ichan}, '_1'];
%     end
%     if config{ielec}.LFP.channel{ichan}(end) == 'X'
%         config{ielec}.LFP.channel{ichan} =  [config{ielec}.LFP.channel{ichan}(1:end-1), '1'];
%     end
% end

% template LFP
LFP = readLFP(config{ielec}, MuseStruct, redo_analysis);
LFP = remove_artefacted_trials(config{ielec}, LFP);

% calculate FFT on sliding timewindow
FFT = FFTtrials(config{ielec}, redo_analysis, LFP);

%read spike data
SpikeRaw = readSpikeRaw_Phy(config{ielec}, redo_analysis);

%read spike waveforms
SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, redo_analysis);
waveformstats  = spikeWaveformStats(config{ielec}, SpikeWaveforms, redo_analysis);

% epoch data into windows
SpikeTrials = readSpikeTrials(config{ielec}, MuseStruct, SpikeRaw, redo_analysis);
SpikeTrials_art_removed = remove_artefacted_trials(config{ielec}, SpikeTrials);

% calculate statistics per window
SpikeStats = spikeTrialStats(config{ielec}, SpikeTrials, true);
SpikeStats_art_removed = spikeTrialStats(config{ielec}, SpikeTrials_art_removed, true);

% create one structure with all data for Mario
preictal_gather_data(config{ielec}, redo_analysis, MuseStruct, LFP, FFT, waveformstats, SpikeStats_art_removed, SpikeTrials_art_removed);

%% plot data over time
% you can plot any metric shown here: SpikeStats{1}.window{1}

for markername = ["window_60", "window_3"]
    iunit = 1:size(SpikeStats{ipart}.(markername), 2); % spike to plot
    ipart = 1;
    
    cfg                 = [];
    cfg.ipart           = ipart;
    cfg.filename        = fullfile(config{ielec}.imagesavedir, 'data_overview', sprintf('%sp%d_data_overview_%s.jpg', config{ielec}.prefix, ipart, markername));
    cfg.minbadlength    = config{ielec}.minbadtime.(markername); %in seconds
    if strcmp(config{ielec}.seizure_index, 'last')
        offset = MuseStruct{1}{2}.markers.CriseStart.clock(end);
    else
        offset = MuseStruct{1}{2}.markers.CriseStart.clock(config{ielec}.seizure_index);
    end
    if markername == "window_60"
        x1 = SpikeTrials_art_removed{1}.(markername).trialinfo.starttime(1);
        x2 = SpikeTrials_art_removed{1}.(markername).trialinfo.endtime(end);
    elseif markername == "window_3"
        x1 = offset - seconds(300);
        x2 = offset + seconds(120);
    end
    
    cfg.xlim            = seconds([x1, x2] - offset);
    % cfg.xlim            = seconds([MuseStruct{1}{1}.starttime, MuseStruct{1}{2}.endtime] - MuseStruct{1}{2}.markers.CriseStart.clock);
    cfg.orientation     = 'landscape';
    cfg.offset          = offset;
    
    cfg.type{1}         = 'trialinfo';
    cfg.title{1}        = sprintf('BAD duration (s) (trial removed if artefacts > %g s)',cfg.minbadlength);
    cfg.metric{1}       = 'BAD_sec';
    
    cfg.type{2}         = 'power';
    cfg.frequency{2}    = [1, 7];
    cfg.channel{2}      = config{ielec}.LFP.channel{1};
    cfg.title{2}        = sprintf('Power %d-%dHz (log)', cfg.frequency{2}(1), cfg.frequency{2}(2));
    cfg.plotart(2)      = true;
    cfg.log(2)          = true;
    cfg.hideart(2)      = true;
    % cfg.marker_indx{1}  = FFT{ipart}.(markername).trialinfo.BAD_cnt>0;
    % cfg.marker_sign{1}  = '.r';
    % cfg.marker_label{1} = 'artefact';
    
    cfg.type{3}         = 'power';
    cfg.frequency{3}    = [8, 14];
    cfg.channel{3}      = config{ielec}.LFP.channel{1};
    cfg.title{3}        = sprintf('Power %d-%dHz (log)', cfg.frequency{3}(1), cfg.frequency{3}(2));
    cfg.plotart(3)      = true;
    cfg.log(3)          = true;
    cfg.hideart(3)      = true;
    
    cfg.type{4}         = 'relpower';
    cfg.frequency1{4}   = [1, 7];
    cfg.frequency2{4}   = [8, 14];
    cfg.channel{4}      = config{ielec}.LFP.channel{1};
    cfg.title{4}        = sprintf('Power (%d-%dHz)/(%d-%dHz)', cfg.frequency1{4}(1), cfg.frequency1{4}(2), cfg.frequency2{4}(1), cfg.frequency2{4}(2));
    cfg.plotart(4)      = true;
    cfg.log(4)          = false;
    cfg.hideart(4)      = true;
    
    cfg.type{5}         = 'spike';
    cfg.title{5}        = sprintf('Firing rate unit %d-%d (log(Hz))', iunit(1), iunit(end));
    cfg.metric{5}       = 'trialfreq';
    cfg.unit{5}         = iunit;
    cfg.plotart(5)      = true;
    cfg.log(5)          = true;
    cfg.hideart(5)      = true;
    
    cfg.type{6}         = 'spike';
    cfg.title{6}        = sprintf('Nr. of bursts per min, unit %d-%d', iunit(1), iunit(end));
    cfg.metric{6}       = 'burst_trialsum';
    cfg.unit{6}         = iunit;
    cfg.plotart(6)      = true;
    cfg.hideart(6)      = true;
    
    cfg.type{7}        = 'spike';
    cfg.title{7}       = sprintf('CV2 unit %d-%d', iunit(1), iunit(end));
    cfg.metric{7}      = 'CV2_trial';
    cfg.unit{7}        = iunit;
    cfg.plotart(7)     = true;
    cfg.log(7)         = false;
    cfg.hideart(7)      = true;
    
    cfg.type{8}         = 'spike';
    cfg.title{8}        = sprintf('CV unit %d-%d', iunit(1), iunit(end));
    cfg.metric{8}       = 'CV_trial';
    cfg.unit{8}         = iunit;
    cfg.plotart(8)      = true;
    cfg.hideart(8)      = true;
    
    plotWindowedData(cfg, MuseStruct, markername, SpikeTrials, SpikeStats, FFT);
    close all;
end

%% plot spike distance over time
% you can plot any metric shown here: SpikeStats{1}.window{1}

for markername = ["window_60", "window_3"]
    iunit = 1:size(SpikeStats{ipart}.(markername), 2); % spike to plot
    ipart = 1;
    
    cfg                 = [];
    cfg.ipart           = ipart;
    cfg.filename        = fullfile(config{ielec}.imagesavedir, 'spike_distance', sprintf('%sp%d_spike_distance_%s.jpg', config{ielec}.prefix, ipart, markername));
    cfg.minbadlength    = config{ielec}.minbadtime.(markername); %in seconds
    if strcmp(config{ielec}.seizure_index, 'last')
        offset = MuseStruct{1}{2}.markers.CriseStart.clock(end);
    else
        offset = MuseStruct{1}{2}.markers.CriseStart.clock(config{ielec}.seizure_index);
    end
    if markername == "window_60"
        x1 = SpikeTrials_art_removed{1}.(markername).trialinfo.starttime(1);
        x2 = SpikeTrials_art_removed{1}.(markername).trialinfo.endtime(end);
    elseif markername == "window_3"
        x1 = offset - seconds(300);
        x2 = offset + seconds(120);
    end
    cfg.xlim            = seconds([x1, x2] - offset);
    cfg.orientation     = 'landscape';
    cfg.offset          = offset;
    
    cfg.type{1}         = 'trialinfo';
    cfg.title{1}        = sprintf('BAD duration (s) (trial removed if artefacts > %g s)',cfg.minbadlength);
    cfg.metric{1}       = 'BAD_sec';
    
    cfg.type{2}         = 'spike';
    cfg.title{2}        = sprintf('Firing rate unit %d-%d (log(Hz))', iunit(1), iunit(end));
    cfg.metric{2}       = 'trialfreq';
    cfg.unit{2}         = iunit;
    cfg.plotart(2)      = true;
    cfg.log(2)          = true;
    cfg.hideart(2)      = true;
    
    i = 2;
    if size(SpikeStats{ipart}.(markername), 2) > 1
        for i_dist = 1:size(SpikeStats{ipart}.(markername), 2)
            i = i+1;
            cfg.type{i}         = 'spike';
            unitlist = iunit;
            unitlist(i_dist) = [];
            cfg.title{i}        = sprintf('Distance %d to %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', iunit(i_dist), unitlist);
            cfg.metric{i}       = 'dist';
            cfg.unit{i}         = iunit(i_dist);
            cfg.index{i}        = 1:size(iunit,2)-1;
            cfg.plotart(i)      = true;
            cfg.hideart(i)      = true;
        end
    else
        cfg.type{3}         = 'spike';
        cfg.title{3}        = sprintf('Cannot compute spike distance with only one spike train');
        cfg.metric{3}       = 'dist';
        cfg.plotart(3)      = false;
        cfg.unit{3}         = [];
    end
    
    plotWindowedData(cfg, MuseStruct, markername, SpikeTrials, SpikeStats, FFT);
    close all;
    
end

%% plot spike waveform
config{ielec}.plotspike.plotraw    = true;
config{ielec}.plotspike.suffix     = '_raw';
config{ielec}.plotspike.isi_lim    = [0 0.05];
config{ielec}.plotspike.img_format = ["png", "pdf"];
plot_spike_waveforms(config{ielec}, "window_60", waveformstats, SpikeStats_art_removed, SpikeWaveforms);

config{ielec}.plotspike.plotraw    = false;
config{ielec}.plotspike.suffix     = '_avg';
config{ielec}.plotspike.img_format = ["png", "pdf"];
plot_spike_waveforms(config{ielec}, "window_60", waveformstats, SpikeStats_art_removed);
