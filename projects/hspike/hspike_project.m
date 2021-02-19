%% Analysis script for Spike x Sleep x SUA project
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end
ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = hspike_setparams;
labels = ["BAD__START__", "BAD__END__", "PHASE_1__START__", "PHASE_1__END__", "PHASE_2__START__", "PHASE_2__END__", "PHASE_3__START__", "PHASE_3__END__", "REM__START__", "REM__END__", "AWAKE__START__", "AWAKE__END__", "NO_SCORE__START__", "NO_SCORE__END__"];

%% General analyses, looping over patients
%      exportHypnogram(config{ipatient})

for ipatient = 1:7
    
    [MuseStruct_orig{ipatient}]                                             = readMuseMarkers(config{ipatient}, false);
    [MuseStruct_aligned{ipatient}]                                          = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);    
    [clusterindx{ipatient}, LFP_cluster{ipatient}]                          = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
    [MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]    = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
    
    % update - add any new artefacts
    MuseStruct_template{ipatient}                                          = updateMarkers(config{ipatient}, MuseStruct_template{ipatient}, labels);
        
    % rename markers to combined markers (see config)
    MuseStruct_combined{ipatient}                                          = editMuseMarkers(config{ipatient}, MuseStruct_template{ipatient});
    
    % focus time period a bit more
    config{ipatient}.epoch.toi.Hspike           = [-0.2  0.8];
    config{ipatient}.epoch.pad.Hspike           = 0.5;
    config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
    config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
    
    % from now on work on manual and combined templates
    itemp = 1;
    config{ipatient}.name = [];
    for markername = string(unique(config{ipatient}.editmarkerfile.torename(:, 2))')
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
    
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_combined{ipatient}, false);
    [LFP{ipatient}]                                                = readLFP(config{ipatient}, MuseStruct_combined{ipatient}, false);
    [TFR{ipatient}]                                                = TFRtrials(config{ipatient}, LFP{ipatient}, false);
    
    % trim files to only those within a hypnogram
    MuseStruct_trimmed  = MuseStruct_combined;
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
    
    % read spike data from Phy as one continuous trial
    SpikeRaw{ipatient}                    = readSpikeRaw_Phy(config_trimmed{ipatient}, false);
    
    if ipatient == 3
        SpikeRaw{ipatient}{1}.template_maxchan(1) = 1;
    end
    
    % SpikeWaveforms_fromRaw{ipatient}              = readSpikeWaveforms_fromRaw(config_trimmed{ipatient}, SpikeRaw{ipatient}, true);
    
    % segment into trials based on IED markers
    SpikeTrials_timelocked{ipatient}      = readSpikeTrials_MuseMarkers(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, false);
    SpikeDensity_timelocked{ipatient}     = spikeTrialDensity(config_trimmed{ipatient}, SpikeTrials_timelocked{ipatient}, false);
    
    % segment into equal periods
    SpikeTrials_windowed{ipatient}        = readSpikeTrials_windowed(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, true);
    SpikeStats_windowed{ipatient}         = spikeTrialStats(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true, '-windowed');
    SpikeWaveforms{ipatient}              = readSpikeWaveforms(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);
    
    
    
    % plotOverviewHspike(config{ipatient}, marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}, ...
    %     SpikeTrials_timelocked{ipatient}, SpikeTrials_windowed{ipatient}, SpikeStats_windowed{ipatient}, ...
    %     SpikeDensity_timelocked{ipatient}, LFP{ipatient}, TFR{ipatient}, SpikeWaveforms{ipatient});
    
    plotOverviewHspike(config{ipatient}, marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}, ...
        SpikeTrials_timelocked{ipatient}, SpikeTrials_windowed{ipatient}, SpikeStats_windowed{ipatient}, ...
        SpikeDensity_timelocked{ipatient}, LFP{ipatient}, TFR{ipatient}, SpikeWaveforms{ipatient});
    
    plotGrandAverageAll(config);
end

%% collect behavioural for overview
config = hspike_setparams;

% hypnogram labels to use
hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

% order of plots and colors
labelorder  = ["REM", "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

for ipatient = 5:7
    [MuseStruct_orig{ipatient}]                                                                     = readMuseMarkers(config{ipatient}, false);
    [MuseStruct_aligned{ipatient}]                                                                  = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
    [MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]                            = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, [], false);
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_aligned{ipatient}, false);
end
fname = fullfile(config{1}.datasavedir, 'compilation');
save(fname, '-v7.3')

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

for ipatient = 1 : 7
    itemp = 1;
    config{ipatient}.name = [];
    for markername = string(unique(config{ipatient}.editmarkerfile.torename(:, 2))')
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
end
config = config(1 : 7);

% load average LFPs
for ipatient = 1 : 7
    
    for markername = string(unique(config{ipatient}.editmarkerfile.torename(:, 2))')
        
        fname = fullfile(config{ipatient}.datasavedir, strcat(config{ipatient}.prefix, 'LFP_', markername, '.mat'));
        
        if ~exist(fname,'file')
            continue
        end
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            fprintf('Loading %s\n', fname);
            try
                temp = load(fname);
                for ipart = 1 : 3
                    LFP{ipatient}{ipart}.(markername) = temp.LFP{ipart}.(markername);
                end
                clear temp
            catch ME
                err_count = err_count + 1;
            end
            count = count + 1;
        end
    end
end


[dat, dat_hyp, trialinfo] = Trials2GrandAverage(config, true);

plotGrandAverage_timelocked(config, LFP, dat, dat_hyp, trialinfo);



%% load windowed data
config = hspike_setparams;

for ipatient = 1 : 7
    
    fname = fullfile(config{ipatient}.datasavedir, [config{ipatient}.prefix, 'SpikeTrials_Windowed.mat']);
    load(fname);
    SpikeTrials_windowed{ipatient}        = readSpikeTrials_windowed(config{ipatient});
    
    config{ipatient}.postfix              = '-windowed';
    SpikeStats_windowed{ipatient}         = spikeTrialStats(config{ipatient});
    
    SpikeWaveforms{ipatient}              = readSpikeWaveforms(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);
    
end




plotGrandAverage_windowed(config, []);

%% Create slurm job list
config              = hspike_setparams;
fname_slurm_joblist = fullfile(config{1}.datasavedir, 'slurm_job_list.txt');
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


%% BrainNetViewer
addpath '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\BrainNetViewer_20191031'

config  = hspike_setparams;
nodes_all = [];

% Separate per patient
for ipatient = 1 : 7
    
    cfg     = config{ipatient};
    fname   = cfg.locfile;
    fid     = fopen(fname);
    dat     = textscan(fid, '%s%d%s%f%f%f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', 'Headerlines', 4, 'delimiter', ';');
    fclose(fid);
    
    clear elec
    toclear = [];
    
    for i = 1 : size(dat{1}, 1)
        elec(i).name = string(dat{1}{i});
        if isempty(dat{1}{i})
            elec(i).name = elec(i-1).name;
        end
        elec(i).plotnr = dat{2}(i);
        elec(i).MNI_x = dat{4}(i);
        elec(i).MNI_y = dat{5}(i);
        elec(i).MNI_z = dat{6}(i);
        elec(i).name_cat = strcat('_', elec(i).name, '_',  num2str(elec(i).plotnr));
        elec(i).color = 0;
        
        for name = string(cfg.LFP.channel)
            if contains(name, elec(i).name_cat)
                elec(i).color = 1;
            end
        end
        
        if elec(i).plotnr == 0
            toclear = [toclear, i];
        end
    end
    elec(toclear) = [];
    
    clear nodes
    for i = 1 : size(elec, 2)
        nodes(i, :) = [elec(i).MNI_x, elec(i).MNI_y, elec(i).MNI_z, elec(i).color, 3, ipatient];
    end

    fname_node  = fullfile(cfg.datasavedir, [cfg.prefix, 'brainnet.node']);
    fname_surf  = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\BrainNetViewer_20191031\Data\SurfTemplate\BrainMesh_Ch2.nv';
    fname_cfg   = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\brainnet_config.mat';
    fname_image = fullfile(cfg.imagesavedir, [cfg.prefix, 'brainnet.jpg']);
    
    dlmwrite(fname_node, nodes, '\t');
    
    % concatinate over patients
    nodes_all   = [nodes_all; nodes];
    
    BrainNet_MapCfg(fname_node, fname_surf, fname_cfg, fname_image); 
end

fname_nodes_all = fullfile(cfg.datasavedir, 'brainnet_all.node');
fname_image_all = fullfile(cfg.imagesavedir, 'brainnet_all.jpg');

% resize nodes that are not analyzed
nodes_all(nodes_all(:, 4) == 0, 5) = 1;

dlmwrite(fname_nodes_all, nodes_all, '\t');
BrainNet_MapCfg(fname_nodes_all, fname_surf, fname_cfg, fname_image_all);



% https://www.nitrc.org/docman/view.php/504/1190/BrainNet%20Viewer%20Manual%201.41.pdf
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there
% are 6 columns: columns 1-3 represent node coordinates, column 4 represents node
% colors, column 5 represents node sizes, and the last column represents node labels.
% Please note, a symbol ‘-‘(no ‘’) in column 6 means no labels. The user may put the
% modular information of the nodes into column 4, like ‘1, 2, 3…’ or other information to
% be shown by color. Column 5 could be set as nodal degree, centrality, T-value, etc. to
% emphasize nodal differences by size. You can generate your nodal file according to the
% requirements.













%% TODO
% Authomatically determine template matching threshold by
% iteration/minimizing/maximizing sensitivity and selectivity compared to
% manual annotation
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