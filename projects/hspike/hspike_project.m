%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses, looping over patients

for ipatient = 2
    
    % load settings
    config = hspike_setparams;
    
    % export hypnogram to muse
    %% REDO AND CHECK 02718; 15_04-31
    exportHypnogram(config{ipatient});
    
    % read muse markers
    [MuseStruct] = readMuseMarkers(config{ipatient}, true);
    
    % align Muse markers according to peaks
    [MuseStruct] = alignMuseMarkers(config{ipatient},MuseStruct, true);
    
    % automatically detect spikes
    detectSpikes(config{ipatient}, MuseStruct, true, true);
    
    % read hypnogram as table
    [PSGtable] = PSG2table(config{ipatient}, MuseStruct, true);
    
    % plot hypnogram
    plotHypnogram(config{ipatient},MuseStruct);
    
    % events vs. hypnogram statistics and plots
    [MuseStruct, marker, hypnogram] = hypnogramStats(config{ipatient}, MuseStruct, false);
     
    % calculate TFR over all files, per part
    TFR = doTFRcontinuous(config{ipatient}, MuseStruct, true);
    
    % plot TFR 
    plotTFRcontinuous(config{ipatient},TFR);
    
    % write data concatinated for SC, artefacts, and output sampleinfo per file
    writeSpykingCircus(config{ipatient}, MuseStruct, true, true);
    
    % write parameters for spyking circus
    writeSpykingCircusParameters(config{ipatient})
    
    % read spike-clustering results, and epoch around events
    [SpikeRaw, SpikeTrials] = readSpykingCircus(config{ipatient}, MuseStruct, false, 'all');
    
    % compute event-related changes of spike rates, and other stats
    [stats_smooth, stats_binned] = spikeratestatsEvents(config{ipatient}, SpikeRaw, SpikeTrials, false);
    
    % read spike-clustering results, and label according to polysomnograpy
    [SpikeRawPSG, SpikeTrialsPSG] = readSpykingCircusPSG(config{ipatient}, MuseStruct, true, 'all');
    
    % stratify according to nr. of windows
    SpikeTrialsPSG = rectifywindownr(config{ipatient},SpikeTrialsPSG);

    % computer spike stats, and label according to polysomnography
    [SpikeStatsPSG] = spikeratestatsPSG(config{ipatient}, SpikeRawPSG, SpikeTrialsPSG, hypnogram, true);
    
    % read LFP data
    [LFP] = readLFP(config{ipatient}, MuseStruct, false, false);
end
     
figure; hold;
for ichan = 1 : size(LFP{1}{1}.label,1)
    subplot(size(LFP{1}{1}.label,1),1,ichan);
    cfg = [];
    cfg.channel = ichan;
    ft_singleplotER(cfg,LFP{1}{1});
end

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'bipolar';
LFP_bipolar = ft_preprocessing(cfg,LFP{1}{1});

figure; hold;
for ichan = 1 : size(LFP_bipolar.label,1)
    subplot(size(LFP_bipolar.label,1),1,ichan);
    cfg = [];
    cfg.channel = ichan;
    ft_singleplotER(cfg,LFP_bipolar);
end
    
    figure; hold;
    for ichan = 1 : size(LFP{1}{1}.label,1)
        plot(LFP{1}{1}.time{1},mean(LFP{1}{1}.trial{1}));
    end
    
    
    
    % append nights
    cfg = [];
    cfg.keepsampleinfo = 'no';
    dat_micro_append{1} = ft_appenddata(cfg,dat_micro{1}{1},dat_micro{2}{1},dat_micro{3}{1});
    dat_macro_append{1} = ft_appenddata(cfg,dat_macro{1}{1},dat_macro{2}{1},dat_macro{3}{1});
    dat_micro_append{2} = ft_appenddata(cfg,dat_micro{1}{2},dat_micro{2}{2},dat_micro{3}{2});
    dat_macro_append{2} = ft_appenddata(cfg,dat_macro{1}{2},dat_macro{2}{2},dat_macro{3}{2});
    
    % average
    dat_micro_avg{1} = ft_timelockanalysis([],dat_micro_append{1});
    dat_macro_avg{1} = ft_timelockanalysis([],dat_macro_append{1});
    dat_micro_avg{2} = ft