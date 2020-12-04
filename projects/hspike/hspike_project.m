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
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/altmany-export_fig-8b0ba13
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/binaries'));
    %     addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    %     addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = hspike_setparams;

%% General analyses, looping over patients
%      exportHypnogram(config{ipatient})

for ipatient = 1:7
    
    [MuseStruct_orig{ipatient}]                                                                     = readMuseMarkers(config{ipatient}, false);
    [MuseStruct_aligned{ipatient}]                                                                  = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
    %     [clusterindx{ipatient}, LFP_cluster{ipatient}]                                                  = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
    % [MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]                            = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
    [MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]                            = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, [], false);
    
    % update to any new artefacts
    [MuseStruct_orig{ipatient}]                                                                     = readMuseMarkers(config{ipatient}, true);
    
    for ipart = 1 : 3
        for idir = 1 : size(MuseStruct_template{ipatient}{ipart}, 2)
            try
                MuseStruct_template{ipatient}{ipart}{idir}.markers.BAD__START__ = MuseStruct_orig{ipatient}{ipart}{idir}.markers.BAD__START__;
                MuseStruct_template{ipatient}{ipart}{idir}.markers.BAD__END__ = MuseStruct_orig{ipatient}{ipart}{idir}.markers.BAD__END__;
                
            catch
            end
            if MuseStruct_template{ipatient}{ipart}{idir}.starttime ~= MuseStruct_orig{ipatient}{ipart}{idir}.starttime
                disp('OD');
            end
            MuseStruct_template{ipatient}{ipart}{idir}.nSamples = MuseStruct_orig{ipatient}{ipart}{idir}.nSamples;
            MuseStruct_template{ipatient}{ipart}{idir}.Fs = MuseStruct_orig{ipatient}{ipart}{idir}.Fs;
        end
    end
    
    switch ipatient
        case 1
            markernames = {'combined1', 'combined2'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined1'; 'template5', 'combined2'; 'template6', 'combined2'};
        case 2
            markernames = {'combined1'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template4', 'combined1'; 'template6', 'combined1'};
        case 3
            markernames = {'combined1', 'combined2'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined2'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined1'};
        case 4
            markernames = {'combined1', 'combined2', 'combined3'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined2'; 'template4', 'combined2'; 'template5', 'combined2'; 'template6', 'combined3'};
        case 5
            markernames = {'combined1', 'combined2', 'combined3'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined3'; 'template5', 'combined2'; };
        case 6
            markernames = {'combined1'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined1'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined1'};
        case 7
            markernames = {'combined1', 'combined2'};
            config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined2'; 'template4', 'combined2'; 'template5', 'combined1'; 'template6', 'combined2'};
    end
    
    MuseStruct_combined{ipatient} = editMuseMarkers(config{ipatient}, MuseStruct_template{ipatient});

    % focus time period a bit more
    config{ipatient}.epoch.toi.Hspike           = [-0.2  0.8];
    config{ipatient}.epoch.pad.Hspike           = 0.5;
    config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
    config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
    
    % from now on work on manual and combined templates
    itemp = 1;
    for markername = string(markernames)
        config{ipatient}.name{itemp}                      = markername;
        config{ipatient}.muse.startmarker.(markername)    = markername;
        config{ipatient}.muse.endmarker.(markername)      = markername;   
        config{ipatient}.epoch.toi.(markername)           = [-0.2  0.8];
        config{ipatient}.epoch.pad.(markername)           = 0.5;
        config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
        config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
        config{ipatient}.LFP.name{itemp}                  = markername;
        config{ipatient}.hyp.markers{itemp}               = markername;
        itemp = itemp + 1;
    end

    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_combined{ipatient}, true);    
    
    [LFP{ipatient}]                                                = readLFP(config{ipatient}, MuseStruct_combined{ipatient}, false);
    [TFR{ipatient}]                                                = TFRtrials(config{ipatient}, LFP{ipatient}, false);

    % trim files to only those within a hypnogram
    MuseStruct_trimmed  = MuseStruct_combined;
%     MuseStruct_trimmed  = MuseStruct_orig;
    config_trimmed      = config;
    for ipart = 1 : 3
        sel     = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
        first   = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
        last    = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
        config_trimmed{ipatient}.directorylist{ipart}   = config{ipatient}.directorylist{ipart}(first:last);
        MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct_combined{ipatient}{ipart}(first:last);
       
        % if still more than 7, cut off the beginning
        if size(config_trimmed{ipatient}.directorylist{ipart}, 2) > 7
            config_trimmed{ipatient}.directorylist{ipart}   = config_trimmed{ipatient}.directorylist{ipart}(end-6:end);
            MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct_trimmed{ipatient}{ipart}(end-6:end);           
        end
    end 
    
    % write data concatinated for SC, artefacts, and output sampleinfo per file
    cfg = writeSpykingCircusMultichannel(config_trimmed{ipatient}, true);
    

    sampleinfo = writeSpykingCircus(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, true, true);
    % write parameters for spyking circus
%     writeSpykingCircusParameters(config_trimmed{ipatient})
%     writeSpykingCircusParameters_new(config_trimmed{ipatient})
%     [artefacts] = writeSpykingCircus_deadfiles(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, true);

    % read spike data from Phy as one continuous trial
    SpikeRaw{ipatient}                    = readSpikeRaw_Phy(config_trimmed{ipatient}, false);
    
    % segment into equal periods
    SpikeTrials_windowed{ipatient}        = readSpikeTrials_windowed(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, true);
    SpikeStats_windowed{ipatient}         = spikeTrialStats(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true, 'windowed');
    
    
    % segment into trials based on IED markers
    SpikeTrials_timelocked{ipatient}      = readSpikeTrials_MuseMarkers(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, true);
    SpikeDensity_timelocked{ipatient}     = spikeTrialDensity(config_trimmed{ipatient}, SpikeTrials_timelocked{ipatient}, true);

    SpikeWaveforms{ipatient}              = readSpikeWaveforms_new(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);
    
    plotOverviewHspike(config{ipatient}, marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}, ...
        SpikeTrials_timelocked{ipatient}, SpikeTrials_windowed{ipatient}, SpikeStats_windowed{ipatient}, ...
        SpikeDensity_timelocked{ipatient}, LFP{ipatient}, TFR{ipatient}, SpikeWaveforms{ipatient});
end

%% collect behavioural for overview

% hypnogram labels to use
hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

% order of plots and colors
labelorder  = ["REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

for ipatient = 1:7
    [MuseStruct_orig{ipatient}]                                                                     = readMuseMarkers(config{ipatient}, false);
    [MuseStruct_aligned{ipatient}]                                                                  = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
    [MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]                            = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, [], false);
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_aligned{ipatient}, false);
end

y       = zeros(7,5);
name    = [];
for ipatient = 1:7 
    ihyp = 1;
    for hyplabel = labelorder
        isum = 0;
        for ipart = 1
            for markername = string(fields(hypmusestat{ipatient}{ipart}))'
                y(ipatient, ihyp) = y(ipatient, ihyp) + hypmusestat{ipatient}{ipart}.(markername).IEDrateNorm.(hyplabel);
                isum = isum + 1;                
            end
        end
        y(ipatient, ihyp) =  y(ipatient, ihyp) / isum;    
        ihyp = ihyp + 1;
    end
end

fig = figure;
set(gcf,'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);

subplot(2,2,2); hold;
title('IED rate vs. wake');
boxplot(y);
set(gca,'Xticklabels',labelorder);
xtickangle(90);

subplot(2,2,1); hold;
title('IED rate vs. wake');
cm = cool(5);
m = -inf;
    
for hyplabel = labelorder
    barloc      = find(hyplabel == labelorder);
    stdev       = nanstd(y(:, barloc));
    hb          = bar(barloc, mean(y(:, barloc)), 1);
    he          = errorbar(barloc, mean(y(:, barloc)), stdev, 'clipping','off','color', 'k');
    l{barloc}   = sprintf('%s=%0.1f(%0.1f)', hyplabel, mean(y(:, barloc)), nanstd(y(:, barloc)));
    m           = max(m,  mean(y(:, barloc)));
    set(hb, 'FaceColor', cm(barloc,:));
end
xlim([0.5, 5.5]); ylim([0, m * 1.1]);
set(gca,'TickLength',[0 0], 'Xticklabels', [], 'TickLabelInterpreter', 'none', 'XGrid', 'off', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');

clear p
for ii = 1:size(cm,1)
    p(ii) = patch(NaN, NaN, cm(ii,:));
end
hl = legend(p, l, 'location', 'eastoutside', 'Interpreter', 'none');
set(findall(fig, '-property', 'fontsize'), 'fontsize', 8);
fname = fullfile(config{1}.imagesavedir, 'overview_relative_IEDrate');
export_fig(fname, '-png'); % need to install https://www.ghostscript.com/download/gsdnld.html
close all


%% collect spikes for overview

config = hspike_setparams;
% 3 = 2660;

for ipatient = [1 5 6]
%     SpikeTrials_timelocked{ipatient}      = readSpikeTrials_MuseMarkers2(config{ipatient}, [], [], false);
%     SpikeDensity_timelocked{ipatient}     = spikeTrialDensity(config{ipatient}, [], false);
    SpikeStats_windowed{ipatient}         = spikeTrialStats(config_trimmed{ipatient}, [], false, 'windowed');
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, [], false);    

end
config{1}.epoch.toi.combined1 = [-0.2  0.8];
config{1}.epoch.toi.combined2 = [-0.2  0.8];
config{1}.epoch.toi.combined2 = [-0.2  0.8];

[GA, stats] = spikeTrialStats_GrandAverage(config, SpikeTrials_timelocked, SpikeDensity_timelocked);
  


%% Create slurm job list
config              = hspike_setparams;
fname_slurm_joblist = fullfile(config{1}.datasavedir, 'slurm_job_list3.txt');
delete(fname_slurm_joblist);
for ipatient = 1:7
    for ipart = 1 : size(config{ipatient}.directorylist, 2)
        subjdir     = config{ipatient}.prefix(1:end-1);
        partdir     = ['p', num2str(ipart)];
        if ~isfield(config{ipatient}.circus, 'channelname')
            fid = fopen(fname_slurm_joblist, 'a');
            if fid == -1
                error('Could not create/open %s', fname_slurm_joblist);
            end
            filename    = strcat(config{ipatient}.prefix, 'p', num2str(ipart), '-multifile-', config{ipatient}.circus.channel{1} ,'.ncs');
            dirname     = strcat('//network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/', subjdir, '/', partdir);
            fprintf(fid,'module load spyking-circus/0.9.9 ; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28; echo DONE!!! \n', dirname, filename, filename);
            fprintf('module load spyking-circus/0.9.9 ; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28; echo DONE!!! \n', dirname, filename, filename);
            fclose(fid);  
        else
            for chandir = unique(config{ipatient}.circus.channelname)
                fid = fopen(fname_slurm_joblist, 'a');
                if fid == -1
                    error('Could not create/open %s', fname_slurm_joblist);
                end
                temp        = strcmp(config{ipatient}.circus.channelname, chandir);
                firstchan   = string(config{ipatient}.circus.channel(find(temp,1,'first')));
                filename    = strcat(config{ipatient}.prefix, 'p', num2str(ipart), '-multifile-', firstchan ,'.ncs');
                dirname     = strcat('//network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/hspike/', subjdir, '/', partdir, '/', string(chandir));
                fprintf(fid,'module load spyking-circus/0.9.9 ; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28; echo DONE!!! \n', dirname, filename, filename);
                fprintf('module load spyking-circus/0.9.9 ; cd %s; spyking-circus %s -c 28; spyking-circus %s -m converting -c 28; echo DONE!!! \n', dirname, filename, filename);
                fclose(fid);
            end
        end
    end
end

%% General analyses
config                                         = hspike_setparams;
[MuseStruct_orig{ipatient}]                    = readMuseMarkers(config{ipatient}, false);
[MuseStruct_aligned{ipatient}]                 = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, true);
[clusterindx{ipatient}, LFP_cluster{ipatient}] = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, true);
[MuseStruct_template{ipatient}, ~,~, ~]        = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, true);

% add templates to config
for itemp = 1 : 6
    markername = sprintf("template%d", itemp);
    config{ipatient}.muse.startmarker.(markername)   = markername;
    config{ipatient}.muse.endmarker.(markername)     = markername;
    config{ipatient}.epoch.toi.(markername)          = [-0.5  1];
    config{ipatient}.epoch.pad.(markername)          = 0.5;
    config{ipatient}.LFP.baselinewindow.(markername) = [-0.5  1];
    config{ipatient}.LFP.baselinewindow.(markername) = [-0.5  1];
    config{ipatient}.LFP.name{itemp}                 = markername;
    config{ipatient}.hyp.markers{itemp}              = markername;
end

%% TODO
% Authomatically determine template matching threshold by
% iteration/minimizing/maximizing sensitivity and selectivity compared to
% manual annotation
% demean before clustering
% use Amygdala as well? often micro + nice SWDs, e.g. in patient 7
% REDO AND CHECK HYPNOGRAM 02718-15_04-31 and 02680_2019-01-16_01-31
% extract spike-by-spike paramewters: EOC, amplitude, durations, etc. to
% put in full model

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
    plotTFRcontinuous(config{ipatient}, TFR);
    
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
    SpikeTrialsPSG = rectifywindownr(config{ipatient}, SpikeTrialsPSG);

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