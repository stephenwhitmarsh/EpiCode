%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\  
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\      
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses, looping over patients
config = hspike_setparams;
for ipatient = 3 : 4   
    [MuseStruct_orig]                       = readMuseMarkers(config{ipatient},                                                             true);
    [MuseStruct_aligned]                    = alignMuseMarkersXcorr(config{ipatient},   MuseStruct_orig,                                    true);
    [clusterindx, LFP_cluster]              = clusterLFP(config{ipatient},              MuseStruct_aligned,                                 true);
    [MuseStruct_template, indx, LFP_avg]    = detectTemplate(config{ipatient},          MuseStruct_aligned, LFP_cluster{1}{1}.kmedoids{6},  true);
    [t]                                     = plotHypnogram(config{ipatient},           MuseStruct_template);
    [marker, hypnogram]                     = hypnogramStats(config{ipatient},          MuseStruct_template,                                true);  
end





%% REDO AND CHECK HYPNOGRAM 02718-15_04-31 and 02680_2019-01-16_01-31

% read LFP data
config{ipatient}.LFP = rmfield(config{ipatient}.LFP, 'resamplefs');

[LFP] = readLFP(config{ipatient}, MuseStruct_orig, true);
     
     % average for 'template'
     cfg = [];
     temp = ft_timelockanalysis(cfg, LFP{1}{1});
     
%      labels_nonum    = regexprep(temp.label, '[0-9_]', '');
%      [~,~,indx]      = unique(labels_nonum);
%      clear group
%      for i = 1 : max(indx)
%          cfgtemp             = [];
%          cfgtemp.reref       = 'yes';
%          cfgtemp.refmethod   = 'bipolar';
%          cfgtemp.channel     = temp.label(indx==i);
%          group{i}            = ft_preprocessing(cfgtemp,temp);
%      end
%      template = ft_appenddata([],group{:});
%      clear group
     
     C = detectTemplate(config{ipatient}, MuseStruct, template, true);
     
     % loop over parts within subject
     for ipart = 1 : size(config{ipatient}.directorylist,2)
         % loop over directories
         C2{ipart} = [];
         for idir = 1 : size(config{ipatient}.directorylist{ipart}, 2)
             C2{ipart} = cat(1,C2{ipart},C{ipart}{idir});
         end
     end
     figure; plot(C2{1}(:,4))
     
%     end
%     % align Muse markers according to peaks
%     [MuseStruct] = alignMuseMarkers(config{ipatient},MuseStruct, false);

    % automatically detect spikes
    detectSpikes(config{ipatient}, MuseStruct, true, true);
    
    % read hypnogram as table
    
    % plot hypnogram
    plotHypnogram(config{ipatient},MuseStruct);
    
    % events vs. hypnogram statistics and plots
     
    % calculate TFR over all files, per part
    TFR = doTFRcontinuous(config{ipatient}, MuseStruct, true);
    
    % plot TFR 
    plotTFRcontinuous(config{ipatient},TFR);
    
    % write data concatinated for SC, artefacts, and output sampleinfo per file
    writeSpykingCircus(config{ipatient}, MuseStruct, true, true);
    
    % write parameters for spyking circus
    writeSpykingCircusParameters(config{ipatient})
    
    % read spike-clustering results, and epoch around events
    [SpikeRaw, SpikeTrials] = readSpykingCircus(config{ipatient}, MuseStruct, true, 1);
    
    % compute event-related changes of spike rates, and other stats
    [stats_smooth, stats_binned] = spikeratestatsEvents(config{ipatient}, SpikeRaw, SpikeTrials, true);
    
    % read spike-clustering results, and label according to polysomnograpy
    [SpikeRawPSG, SpikeTrialsPSG] = readSpykingCircusPSG(config{ipatient}, MuseStruct, true, 'all');
    
    % stratify according to nr. of windows
    SpikeTrialsPSG = rectifywindownr(config{ipatient},SpikeTrialsPSG);

    % computer spike stats, and label according to polysomnography
    [SpikeStatsPSG] = spikeratestatsPSG(config{ipatient}, SpikeRawPSG, SpikeTrialsPSG, hypnogram, true);
    
    % read LFP data
    [LFP] = readLFP(config{ipatient}, MuseStruct, false);
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