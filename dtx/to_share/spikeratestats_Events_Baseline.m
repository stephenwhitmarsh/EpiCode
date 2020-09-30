function [stats] = spikeratestats_Events_Baseline(cfg,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [stats] = spikeratestats_Events_Baseline(cfg,force)
%
% Compare spike activity during short events versus resting periods.
% Can be used with only baseline data (no events) : set 
% cfg.spike.events_name as 'no'.
% Note : if no events, baseline data must still be cut in trials (eg trials
% of the same length along all the data). Use readSpikeTrials_continuous.m 
% to do that.
% Can compute some stats on LFP and spike waveforms 
%
% ### Necessary inputs to load precomputed stats:
% force                         = whether to redo analyses/plot or read
%                                 previous saved stats file (true/false).
%                                 If force = false and stat file already 
%                                 exists, you just need cfg.prefix and 
%                                 cfg.datasavedir to find and load stat file. 
% cfg.prefix                    = prefix to output files
% cfg.datasavedir               = data directory of results
% 
% ### Necessary inputs to add, to compute stats on baseline :
% cfg.SpikeTrials{parts}.(name) = spike data epoched in FieldTrip trial
%                                 data structure
% cfg.name{names}               = names of the different analysis
% cfg.spike.baseline_name        = name of the analysis which correspond to
%                                 baseline/resting data
% cfg.imagesavedir              = where to save images
% cfg.circus.channel            = channel names which were analyzed by
%                                 Spyking-Circus
% 
% ### Optionnal inputs to add, if want to compute stats on events :
% cfg.spike.events_name          = name(s) of the analysis which will be
%                                 processed as events. Can be 'no' to only 
%                                 use baseline data. Default = 'no'.
% cfg.spike.resamplefs{names}   = defines the frequency of bins in the
%                                 event's spikerate barplot
% cfg.spike.sdflatency.(event_name)  = 'maxperiod' or array of 2 x values
% cfg.spike.sdftimwin.(event_name)   = time window for convolution of spike density.
% cfg.stats.bltoi.(name)        = baseline period for statistics
%                                 measurements of events
% cfg.stats.actoi.(name)        = active period for statistics measurements
%                                 of events
% cfg.stats.alpha               = threshold for selecting significant
%                                 clusters. Default = 0.025
% cfg.stats.numrandomization    = nr of randomizations in the cluster-based
%                                 permutation test. Default = 1000.
% cfg.statstime.timewin             = 10;
% cfg.statstime.slidestep           = 2;
% cfg.statstime.removeempty         = 'yes';
%
% ### Optional cfg fields if want to compute stats on LFP and spike waveforms
% cfg.SpikeWaveforms{parts}{names}= spike wavefoms epoched in FieldTrip trial
%                                 data structure. Default = [].
% cfg.dataLFP{parts}.(name)     = LFP data epoched in FieldTrip trial data
%                                 structure. Default = [].
% cfg.LFP.name                  = name of the LFP analysis (to select which
%                                 one corresponds to cfg.spike.eventname)
% cfg.epoch.toi                 = time period to plot, for LFP data
% 
% ### Optional fields, with default values :
% cfg.circus.postfix            = string postfix appended to spike data 
%                                 results. Default = [].
% cfg.spike.ISIbins             = bins limits to plot ISI histogram, in
%                                 seconds. Default = [0:0.003:0.150].
% cfg.spike.RPV                 = definition of refractory period
%                                 violation, in seconds. Default = 0.003
% cfg.stats.part_list           = list of parts to analyse. Can be an array 
%                                 of integers, or 'all'. Default = 'all'. 
%
% ### OUTPUT
% stats{parts}.(name)       = MATLAB structure with one field per type
%                                 of period analyzed, with infos computed for
%                                 each unit
%    - method                   = 'event' or 'baseline'
%    - label                    = label of each unit
%    - isi                      = output from ft_spike_isi
%    - sdf                      = output from ft_spikedensity
%    - clusterstats             = output from ft_timelockanalysis, with
%                                 additional infos about max cluster and min
%                                 cluster (position, value, percent of
%                                 increase/decrease)
%    - freq_trialavg{i_unit}    = avg freq for each trial, per unit
%    - amplitude_trialavg{i_unit}= avg amplitude for each trial, per unit
%    - template{i_unit}         = measured halfwidth, peaktrough and
%                                 troughpeak for each template, in seconds
%    - LFP.channel /.halfwidth  = measured halfwidth on LFP events, for
%                                 each analyzed channel.
%    - firing{i_unit}        = spike train descriptive stats : meanISI,
%                                 stdISI, meanfreq, %RPV, meanCV2, stdCV2
%    - spikewaveform{i_unit}    = measured halfwidth, peaktrough and
%                                 troughpeak on waveforms for each unit
%    - stats_over_time          = freq and cv2 (with and without bursts) 
%                                 over time for each trial and each unit. 
%                                 For each : cfg, values, time, std, avg.
% 
% ### Dependencies
% - Fieldtrip
% - To create the input structures : 
%       > readSpikeTrials_MuseMarkers.m or readSpikeTrials_continuous
%       > readLFP.m
%       > readSpikeWaveforms.m
% - Used in this script : plot_morpho.m, spikestatsOverTime.m
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load precomputed stats if required
fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikeratestats_events_baseline', cfg.circus.postfix, '.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed spike stats data\n');
    load(fname,'stats');
    return
end
 
% get default cfg fields
cfg.spike                   = ft_getopt(cfg         , 'spike'              	, []);
cfg.stats                   = ft_getopt(cfg         , 'stats'             	, []);
cfg.circus                  = ft_getopt(cfg         , 'circus'           	, []);
cfg.SpikeWaveforms          = ft_getopt(cfg         , 'SpikeWaveforms'    	, []);
cfg.dataLFP                 = ft_getopt(cfg         , 'dataLFP'             , []);
cfg.spike.events_name        = ft_getopt(cfg.spike   , 'events_name'          , 'no');
cfg.spike.ISIbins           = ft_getopt(cfg.spike   , 'ISIbins'             , [0:0.003:0.150]);
cfg.stats.alpha             = ft_getopt(cfg.stats   , 'alpha'               , 0.025);
cfg.stats.numrandomization  = ft_getopt(cfg.stats   , 'numrandomization'    , 1000);
cfg.stats.part_list         = ft_getopt(cfg.stats   , 'part_list'           , 'all');
cfg.circus.postfix          = ft_getopt(cfg.circus  , 'postfix'             , []);

if isempty(cfg.spike.events_name)
    cfg.spike.events_name = 'no';
end

if strcmp(cfg.stats.part_list,'all')
    cfg.stats.part_list = 1:size(cfg.SpikeTrials,2);
end

%% Find event and baseline labels

cfg.spike.events_name    = string(ft_getopt(cfg.spike, 'events_name', []));
if all(strcmp(cfg.spike.events_name,'no'))
    cfg.spike.events_name = [];
end
if isempty(cfg.spike.events_name)
    hasevent = false;
else
    hasevent = true;
end

baseline_name            = string(ft_getopt(cfg.spike, 'baseline_name', []));
if isempty(baseline_name)
    error('Need one baseline label.');
end
if length(baseline_name)>1
    error('Only one baseline is allowed for now\n');
end

%% compute stats over time
cfg.statstime.label_list    = [cfg.spike.events_name, baseline_name];
cfg.statstime.write         = 'yes';%save on disk 
statsovertime               = spikestatsOverTime(cfg, cfg.SpikeTrials,true);
if hasevent
    cfg.statstime.removebursts  = 'yes';
    cfg.statstime.suffix        = '_withoutbursts';
    statsovertime_withoutbursts = spikestatsOverTime(cfg, cfg.SpikeTrials,true);
end

%% go trhough each part and labels
for ipart = cfg.stats.part_list
    
    spike_Fs = cfg.SpikeTrials{ipart}.(baseline_name).hdr.Fs;
    %add a dummy event name if no event, for the for loop
    if isempty(cfg.spike.events_name) && ~hasevent
        cfg.spike.events_name = "dummy";
    end
    
    for event_name = cfg.spike.events_name
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% COMPUTE STATS %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % keep infos about the analysis
        stats{ipart}.(baseline_name).method          = 'baseline';
        if hasevent, stats{ipart}.(event_name).method     = 'event'; end
        stats{ipart}.(baseline_name).label          = cfg.SpikeTrials{ipart}.(baseline_name).label;
        if hasevent, stats{ipart}.(event_name).label     = cfg.SpikeTrials{ipart}.(event_name).label; end
        
        %% Compute ISI event, ISI baseline
        cfgtemp                                             = [];
        cfgtemp.outputunit                                  = 'spikecount';
        cfgtemp.bins                                        = cfg.spike.ISIbins;
        cfgtemp.param                                       = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials                                  = 'yes';
        stats{ipart}.(baseline_name).isi         = ft_spike_isi(cfgtemp,cfg.SpikeTrials{ipart}.(baseline_name));
        if hasevent, stats{ipart}.(event_name).isi    = ft_spike_isi(cfgtemp,cfg.SpikeTrials{ipart}.(event_name)); end
   
        
        %% Compute PSTH for events
        if hasevent
            
            clear sdf_event sdf_baseline
            
            cfgtemp                         = [];
            cfgtemp.binsize                 = cfg.spike.psthbin.(event_name);   
            cfgtemp.keeptrials              = 'yes';
            [psth_event]                     = ft_spike_psth(cfgtemp,cfg.SpikeTrials{ipart}.(event_name));
            stats{ipart}.(event_name).psth = psth_event;
            
            %% compute spike density
            cfgtemp                         = [];
            cfgtemp.fsample                 = cfg.spike.resamplefs.(event_name);   % sample at 1000 hz
            cfgtemp.keeptrials              = 'yes';
            %cfgtemp.timwin                  = cfg.spike.toispikerate.(event_name); %smoothing parameter
            [sdfavg, sdfdata]               = ft_spikedensity(cfgtemp,cfg.SpikeTrials{ipart}.(event_name));
            stats{ipart}.(event_name).sdfavg  = sdfavg;
            stats{ipart}.(event_name).sdfdata = sdfdata;
            
            
            %% Compute cluster stats per unit for events
            
            % prepare dummy data with baseline value per trial for stats
            slim(1)                         = find(psth_event.time > cfg.stats.bltoi.(event_name)(1), 1, 'first');
            slim(2)                         = find(psth_event.time < cfg.stats.bltoi.(event_name)(2), 1, 'last');
            sdf_bl                          = psth_event;
            sdf_bl.trial                    = ones(size(psth_event.trial)) .* nanmean(psth_event.trial(:,:,slim(1):slim(2)),3); % replace with mean
            
            for i_unit = 1 : size(cfg.SpikeTrials{ipart}.(event_name).label, 2)
                
                % statistics
                cfgtemp = [];
                cfgtemp.statistic                       = 'ft_statfun_depsamplesT';
                cfgtemp.alpha                           = cfg.stats.alpha;
                cfgtemp.clusteralpha                    = 0.01;
                cfgtemp.method                          = 'montecarlo';
                cfgtemp.computestat                     = 'yes';
                cfgtemp.correctm                        = 'cluster';
                cfgtemp.latency                         = [cfg.stats.bltoi.(event_name)(2) psth_event.time(end)]; % active period starts after baseline
                cfgtemp.ivar                            = 1;
                cfgtemp.uvar                            = 2;
                cfgtemp.design(1,:)                     = [ones(1,size(psth_event.trial,1)) ones(1,size(sdf_bl.trial,1)) *2];
                cfgtemp.design(2,:)                     = [1 : size(psth_event.trial,1) 1 : size(sdf_bl.trial,1)];
                cfgtemp.numrandomization                = cfg.stats.numrandomization;
                cfgtemp.channel                         = i_unit;
                stats{ipart}.(event_name).clusterstat{i_unit} = ft_timelockstatistics(cfgtemp,psth_event,sdf_bl);
                
                % calculate baseline
                slim(1) = find(psth_event.time > cfg.stats.bltoi.(event_name)(1), 1, 'first');
                slim(2) = find(psth_event.time < cfg.stats.bltoi.(event_name)(2), 1, 'last');
                stats{ipart}.(event_name).clusterstat{i_unit}.bl.avg        = nanmean(psth_event.avg(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(event_name).clusterstat{i_unit}.bl.var        = nanmean(psth_event.var(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(event_name).clusterstat{i_unit}.bl.dof        = nanmean(psth_event.dof(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(event_name).clusterstat{i_unit}.bl.trialavg   = nanmean(psth_event.trial(:,i_unit,slim(1):slim(2)),3);
            end
            
        end %hasevent
        
        %% %%%%%%%%%%%%%%%%%%%%
        %%%%%%%% PLOT %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        for i_unit = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).label, 2)
            
            fig = figure;
            %title for the whole figure
            phy_channr      = cfg.SpikeTrials{ipart}.(baseline_name).template_maxchan(i_unit);
            nlx_channame    = cfg.circus.channel{phy_channr+1}; %+1 because starts at zero
            sgtitle(sprintf('Electrode %d (%s) : %s',phy_channr,nlx_channame, cfg.SpikeTrials{ipart}.(baseline_name).label{i_unit}), 'Fontsize', 22, 'Interpreter', 'none', 'FontWeight', 'bold');
            
            %% Plot firing rate along all data
            if hasevent, subplot(7,4,1:3);hold; else, subplot(4,4,1:3);hold; end
            
            %for baseline
            for itrial = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo, 1)
                x = [statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(1) statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(end)];
                y = nanmean(statsovertime{ipart}.(baseline_name).freq{i_unit}{itrial});
                legend_Bl = plot(x, [y y], 'b', 'LineWidth', 2);
            end
            
            if hasevent
                %for Event
                for itrial = 1:size(cfg.SpikeTrials{ipart}.(event_name).trialinfo, 1)
                    x = statsovertime{ipart}.(event_name).time{i_unit}{itrial}(end);
                    y = nanmean(statsovertime{ipart}.(event_name).freq{i_unit}{itrial});
                    legend_Ev = scatter(x, y, 7, 'o', 'filled', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r');
                end
            end
            
            axis tight
            ax = axis;
            xlim([0 Inf]);
            xticks(0:3600:ax(2));
            xticklabels(xticks/3600); %convert seconds to hours
            xlabel('Time (hours)');
            set(gca, 'YColor', 'k');
            ylabel(sprintf('Mean firing rate \nof each trial (log10)'));
            set(gca, 'YScale', 'log');
            setfig();
            
            if hasevent, legend([legend_Bl, legend_Ev],convertStringsToChars(baseline_name),convertStringsToChars(event_name),'location','eastoutside'); end
            if ~hasevent, legend(legend_Bl,convertStringsToChars(baseline_name),'location','eastoutside'); end
            
            %% Plot amplitude along all data
            if hasevent, subplot(7,4,5:7);hold; else, subplot(4,4,5:7);hold; end
            
            %for baseline
            for itrial = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo, 1)
                x = [statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(1) statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(end)];
                y = nanmean(statsovertime{ipart}.(baseline_name).amplitude{i_unit}{itrial});
                legend_Bl = plot(x, [y y], 'b', 'LineWidth', 2);
            end
            
            if hasevent
                %for Event
                for itrial = 1:size(cfg.SpikeTrials{ipart}.(event_name).trialinfo, 1)
                    x = statsovertime{ipart}.(event_name).time{i_unit}{itrial}(end);
                    y = nanmean(statsovertime{ipart}.(event_name).amplitude{i_unit}{itrial});
                    legend_Ev = scatter(x, y, 7, 'o', 'filled', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r');
                end
            end
            
            axis tight
            xticks(0:3600:ax(2));
            xticklabels(xticks/3600); %convert seconds to hours
            xlim([0 ax(2)]);
            xlabel('Time (hours)');
            set(gca, 'YColor', 'k');
            ylabel(sprintf('Mean amplitude \nof each trial'));
            setfig();
            
            if hasevent, legend([legend_Bl, legend_Ev],convertStringsToChars(baseline_name),convertStringsToChars(event_name),'location','eastoutside'); end
            if ~hasevent, legend(legend_Bl,convertStringsToChars(baseline_name),'location','eastoutside'); end
            

            %% Template
            if hasevent, subplot(7,4,[4 8]); hold; else, subplot(8,4,[8 12]); hold; end
            
            tempsel = cfg.SpikeTrials{ipart}.(baseline_name).template{i_unit}(:,cfg.SpikeTrials{ipart}.(baseline_name).template_maxchan(i_unit)+1,:);%+1 because electrodes nr are zero-based
            temptime = ( (0 : size(cfg.SpikeTrials{ipart}.(baseline_name).template{i_unit},3) - 1) / spike_Fs )';
            %tempsel dimensions : 1:ntemplates 2:maxchan 3:values
            
            % interpolate template
            temptime_int = linspace(temptime(1),temptime(end),10000);
            tempsel_int  = pchip(temptime,tempsel,temptime_int);
            
            %convert to Fieldtrip 'raw' type
            template = [];
            for itrial=1:size(tempsel,1)
                template.time{itrial}   = temptime_int;
                template.trial{itrial}  = tempsel_int(itrial,:,:);
                %permute to remove one of the 2 units dimensions in case of several templates
                if size(template.trial{itrial}, 2) == 1
                    template.trial{itrial} = permute(template.trial{itrial},[1 3 2]);
                end
            end
            template.label{1} = 'template';
            
            cfgtemp                     = [];
            cfgtemp.morpho.channame            = 'template';
            cfgtemp.morpho.measurehalfwidth     = 'yes';
            cfgtemp.morpho.blmethod     = 'min'; 
            cfgtemp.morpho.measurepeaktrough    = 'yes';
            cfgtemp.morpho.toiac               = 'all';
            cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.blmethod = 'min'; 
            [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,template);
            
            title([]);
            ylabel('Template (�V)');
            xlabel('Time (ms)');
            xticklabels(xticks*1000);
            setfig();
            
            stats{ipart}.(baseline_name).template.halfwidth{i_unit}    = halfwidth;
            stats{ipart}.(baseline_name).template.peaktrough{i_unit}   = peaktrough;
            stats{ipart}.(baseline_name).template.troughpeak{i_unit}   = troughpeak;
            stats{ipart}.(baseline_name).template.time{i_unit}         = template.time;
            stats{ipart}.(baseline_name).template.values{i_unit}       = template.trial;
            
            if hasevent
                stats{ipart}.(event_name).template.halfwidth{i_unit}              = halfwidth;
                stats{ipart}.(event_name).template.peaktrough{i_unit}             = peaktrough;
                stats{ipart}.(event_name).template.troughpeak{i_unit}             = troughpeak;
                stats{ipart}.(event_name).template.temptime_int{i_unit}           = template.time;
                stats{ipart}.(event_name).template.tempval_int{i_unit}            = template.trial;
            end
            
            %% Events LFP (same chan as spikes). Ignored if no LFP data
            if hasevent
                
                subplot(7,4,[13 14]);hold;
                
                %find if there is LFP data and correct channel name
                plotLFP = false;
                if any(strcmp(fieldnames(cfg.dataLFP{ipart}),event_name))
                    if ~isempty(cfg.dataLFP{ipart}.(event_name))
                        plotLFP = true;
                        %correct channel name for my data
                        if ~ismember (nlx_channame, cfg.dataLFP{ipart}.(event_name).label)
                            nlx_channame = sprintf('%sLFP',nlx_channame);
                        end
                        if ~ismember (nlx_channame, cfg.dataLFP{ipart}.(event_name).label)
                            warning('Cannot find LFP channel for %s in %s', cfg.SpikeTrials{ipart}.(event_name).label{i_unit}, convertStringsToChars(event_name));
                            plotLFP = false;
                        end
                    end
                end

                if plotLFP
                    %correct baseline
                    cfgtemp = [];
                    cfgtemp.demean = 'yes';
                    cfgtemp.baselinewindow = cfg.stats.bltoi.(event_name);
                    cfg.dataLFP{ipart}.(event_name) = ft_preprocessing(cfgtemp, cfg.dataLFP{ipart}.(event_name));
                    
                    %plot
                    cfgtemp                     = [];
                    cfgtemp.morpho.channame            = nlx_channame;
                    cfgtemp.morpho.plotstd             = 'yes';
                    cfgtemp.morpho.removeoutliers      = 'yes';
                    cfgtemp.morpho.toiplot             = cfg.epoch.toi.(event_name);
                    cfgtemp.morpho.toibl               = cfg.stats.bltoi.(event_name);
                    cfgtemp.morpho.toiac               = cfg.stats.actoi.(event_name);
                    cfgtemp.morpho.measurehalfwidth     = 'yes';
                    cfgtemp.morpho.blmethod     = 'bl';
                    cfgtemp.morpho.name                = event_name;
                    [hw_lfp, ~, ~] = plot_morpho(cfgtemp,cfg.dataLFP{ipart}.(event_name));
                    
                    stats{ipart}.(event_name).LFP.channel{strcmp(cfg.dataLFP{ipart}.(event_name).label,nlx_channame)'}     = nlx_channame;
                    stats{ipart}.(event_name).LFP.halfwidth{strcmp(cfg.dataLFP{ipart}.(event_name).label,nlx_channame)'}   = hw_lfp;
                    
                    xlabel([]);
                    titlepos = title(sprintf('\n%s : %d trials, %d spikes',convertStringsToChars(event_name), size(cfg.SpikeTrials{ipart}.(event_name).trialinfo,1), size(cfg.SpikeTrials{ipart}.(event_name).trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                    titlepos.Position(1) = cfg.epoch.toi.(event_name)(1);
                    ylabel('LFP (�V)');
                    setfig();
                end
                
                if ~plotLFP
                    titlepos = title(sprintf('\n%s : %d trials, %d spikes',convertStringsToChars(event_name), size(cfg.SpikeTrials{ipart}.(event_name).trialinfo,1), size(cfg.SpikeTrials{ipart}.(event_name).trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                    titlepos.Position(1) = 0;
                    axis off;
                end
                
                %% Events raster plot
                subplot(7,4,[17 18]);hold;
                cfgtemp                 = [];
                cfgtemp.spikechannel    = i_unit;
                cfgtemp.latency         = [psth_event.time(1) psth_event.time(end)];
                cfgtemp.trialborders    = 'yes';
                ft_spike_plot_raster(cfgtemp,cfg.SpikeTrials{ipart}.(event_name));
                xlabel([]);
                setfig();
                
                %% Events barplot and stats
                subplot(7,4,[21 22]); hold;
                
                %plot avg and std
                avg = bar(psth_event.time,psth_event.avg(i_unit,:),1,'edgecolor','none','facecolor',[0 0 0]);
                sd  = plot(psth_event.time,sqrt(psth_event.var(i_unit,:))+psth_event.avg(i_unit,:),'-','LineWidth',2,'Color',[0.6 0.6 0.6]);
                
                ylabel('Spikerate (Hz)');
                axis tight
                ax = axis;
                
                baseline = stats{ipart}.(event_name).clusterstat{i_unit}.bl.avg;
                idx_begin_stat = find(psth_event.time > cfg.stats.bltoi.(event_name)(2), 1, 'first');
                
                % plot positive clusters
                if isfield(stats{ipart}.(event_name).clusterstat{i_unit},'posclusters')
                    for ipos = 1 : size(stats{ipart}.(event_name).clusterstat{i_unit}.posclusters,2)
                        if stats{ipart}.(event_name).clusterstat{i_unit}.posclusters(ipos).prob < cfg.stats.alpha
                            sel = [];
                            sel = find(stats{ipart}.(event_name).clusterstat{i_unit}.posclusterslabelmat == ipos);
                            if length(sel) == 1 %correct a bug in the plot if only one sample is positive
                                barwidth = stats{ipart}.(event_name).clusterstat{i_unit}.time(2)-stats{ipart}.(event_name).clusterstat{i_unit}.time(1); 
                            else
                                barwidth = 1;
                            end
                            bar(stats{ipart}.(event_name).clusterstat{i_unit}.time(sel),psth_event.avg(i_unit,sel+idx_begin_stat-1),barwidth,'facecolor','g','edgecolor','g');
                            
                            % compute percentage of max
                            [~,max_idx_temp] = max(psth_event.avg(i_unit,sel+idx_begin_stat-2));
                            max_idx = sel(max_idx_temp) + idx_begin_stat - 2;
                            x_max = psth_event.time(max_idx);
                            y_max = psth_event.avg(i_unit,max_idx);
                            percent_increase = (psth_event.avg(i_unit,max_idx)-baseline) / baseline * 100;
                            %                             y = ax(4)*0.9;
                            %                             text(x_max,y,sprintf('max = %.1fHz \n(+%.1f%%)',y_max,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                            %                             plot(x_max,y_max,'*b');
                            
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.perc{ipos} = percent_increase;
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.x{ipos} = x_max;
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.y{ipos} = y_max;
                        end
                    end
                end
                
                %plot max
                max_freq = max(psth_event.avg(i_unit,:));
                percent_increase = (max_freq-baseline) / baseline * 100;
                x = (cfg.stats.actoi.(event_name)(1) + cfg.stats.actoi.(event_name)(2))/2;
                y = ax(4)*0.9;
                text(x,y,sprintf('max = %.1fHz \n(+%.1f%%)',max_freq,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                stats{ipart}.(event_name).maxfreq.max{i_unit} = max_freq;
                stats{ipart}.(event_name).maxfreq.perc_increase{i_unit} = percent_increase;
                
                % plot negative clusters
                if isfield(stats{ipart}.(event_name).clusterstat{i_unit},'negclusters')
                    for ineg = 1 : size(stats{ipart}.(event_name).clusterstat{i_unit}.negclusters,2)
                        if stats{ipart}.(event_name).clusterstat{i_unit}.negclusters(ineg).prob < cfg.stats.alpha
                            sel = [];
                            sel = find(stats{ipart}.(event_name).clusterstat{i_unit}.negclusterslabelmat == ineg);
                            if length(sel) == 1 %correct a bug in the plot if only one sample is positive
                                barwidth = stats{ipart}.(event_name).clusterstat{i_unit}.time(2)-stats{ipart}.(event_name).clusterstat{i_unit}.time(1);
                            else
                                barwidth = 1;
                            end
                            bar(stats{ipart}.(event_name).clusterstat{i_unit}.time(sel),psth_event.avg(i_unit,sel+idx_begin_stat-1),barwidth,'facecolor','r','edgecolor','r');
                            
                            % compute percentage
                            [~,min_idx_temp] = min(psth_event.avg(i_unit,sel+idx_begin_stat-2));
                            min_idx = sel(min_idx_temp) + idx_begin_stat - 2;
                            x_min = psth_event.time(min_idx);
                            y_min = psth_event.avg(i_unit,min_idx);
                            percent_decrease = (psth_event.avg(i_unit,min_idx)-baseline) / baseline * 100;
                            %                     %text(x_min,y_min+y_min/100,sprintf('%.1fHz (%.1f%%)\n',y_min,percent_decrease),'HorizontalAlignment','center','VerticalAlignment','middle');
                            %                     %plot(x_min,y_min,'*b');
                            
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.perc{ineg} = percent_decrease;
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.x{ineg} = x_min;
                            stats{ipart}.(event_name).clusterstat{i_unit}.maxcluster.y{ineg} = y_min;
                        end
                    end 
                end
                
                %plot sdf
                plot(sdfavg.time, sdfavg.avg(i_unit,:), 'b');
                plot(sdfavg.time, sqrt(sdfavg.var(i_unit,:))+sdfavg.avg(i_unit,:),'b');
                
                % plot baseline patch
                x = [cfg.stats.bltoi.(event_name)(1) cfg.stats.bltoi.(event_name)(2) cfg.stats.bltoi.(event_name)(2) cfg.stats.bltoi.(event_name)(1)];
                y = [ax(3) ax(3) ax(4) ax(4)];
                patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
                
                % plot baseline
                x = [];
                x = [ax(1) ax(2)];
                plot(x,[baseline, baseline],':k');
                
                % plot baseline text
                x = (cfg.stats.bltoi.(event_name)(1) + cfg.stats.bltoi.(event_name)(2))/2;
                y = ax(4)*0.9;
                d = stats{ipart}.(event_name).clusterstat{i_unit}.bl.avg;
                text(x,y,sprintf('baseline = %.1f Hz',d),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                
                legend([avg, sd],'Average','SD');
                xlabel(sprintf('Time from %s (s)', convertStringsToChars(event_name)),'FontSize',10);
                setfig();
                
            end %hasevent
            
            
            %% ISI for event and baseline
            
            if hasevent, plotlist = [event_name, baseline_name]; else, plotlist = baseline_name; end
            
            for iplot = plotlist
                
                if strcmp(iplot, event_name) && hasevent
                    subplot(7,4,25);hold;
                elseif strcmp(iplot,baseline_name)
                    if hasevent, subplot(7,4,27);hold; else, subplot(8,4,[21 22 25 26 29 30]);hold; end
                end
                
                bar(stats{ipart}.(iplot).isi.time,stats{ipart}.(iplot).isi.avg(i_unit,:), 'EdgeColor','none');
                
                RPV         = (length(find(stats{ipart}.(iplot).isi.isi{i_unit} < cfg.spike.RPV)) / length(stats{ipart}.(iplot).isi.isi{i_unit})) * 100;
                meanISI     = nanmean(stats{ipart}.(iplot).isi.isi{i_unit})*1000;
                stdISI      = nanstd(stats{ipart}.(iplot).isi.isi{i_unit})*1000;
                meanfreq    = 1/(meanISI/1000);
                
                %cv2
                arraytemp = stats{ipart}.(iplot).isi.isi{i_unit};
                if length(arraytemp)>2
                    for i = 1:length(arraytemp)-1
                        cv2_data(i) = 2*abs(arraytemp(i)-arraytemp(i+1))/(arraytemp(i)+arraytemp(i+1));
                    end
                    meancv2             = nanmean(cv2_data);
                    stdcv2              = nanstd(cv2_data);
                else
                    meancv2             = NaN;
                    stdcv2              = NaN;
                end
                clear arraytemp
                
                yticklabels(yticks); 
                xticklabels(xticks*1000); %convert in ms
                ylabel('Spike count');
                xlh =xlabel('ISI (ms)');
                text(0,xlh.Position(2),sprintf('RPV = %.1f%% \n<ISI> = %.1f +/- %.1f ms \n<F> = %.1f Hz \nCV2 = %.2f +/- %.2f',RPV, meanISI, stdISI, meanfreq, meancv2, stdcv2),'HorizontalAlignment','left', 'VerticalAlignment','top','FontWeight','bold');
                
                setfig();
                
                stats{ipart}.(iplot).firing.meanISI{i_unit}   = meanISI;
                stats{ipart}.(iplot).firing.stdISI{i_unit}    = stdISI;
                stats{ipart}.(iplot).firing.meanfreq{i_unit}  = meanfreq;
                stats{ipart}.(iplot).firing.RPV{i_unit}       = RPV;
                stats{ipart}.(iplot).firing.meancv2{i_unit}   = meancv2;
                stats{ipart}.(iplot).firing.stdcv2{i_unit}    = stdcv2;
            end
            
            if ~hasevent
                titlepos = title(sprintf('\n%s : %d spikes\n',convertStringsToChars(baseline_name), size(cfg.SpikeTrials{ipart}.(baseline_name).trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                titlepos.Position(1) = 0;
            end
            
            %% Spike Waveforms for events and baseline
            if hasevent, plotlist = [event_name, baseline_name]; else, plotlist = baseline_name; end
            
            for iplot = plotlist
                
                if strcmp(iplot, event_name) && hasevent
                    subplot(7,4,26);hold;
                elseif strcmp(iplot,baseline_name)
                    if hasevent, subplot(7,4,28);hold; else, subplot(8,4,[23 24 27 28 31 32]);hold; end
                end
                
                if ~isempty(cfg.SpikeWaveforms)
                    if ~isempty(cfg.SpikeWaveforms{ipart}.(iplot){i_unit})
                        
                        cfgtemp                     = [];
                        cfgtemp.morpho.channame            = cfg.SpikeWaveforms{ipart}.(iplot){i_unit}.label{1};
                        cfgtemp.morpho.plotstd             = 'yes';
                        cfgtemp.morpho.removeoutliers      = 'yes'; %if big noise, impair seeing real data. Still present in avg and std.
                        cfgtemp.morpho.measurehalfwidth     = 'yes';
                        cfgtemp.morpho.blmethod            = 'min';
                        cfgtemp.morpho.measurepeaktrough    = 'yes';
                        cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.blmethod     = 'min';
                        cfgtemp.morpho.toiac               = 'all';
                        cfgtemp.morpho.name                = iplot;
                        [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,cfg.SpikeWaveforms{ipart}.(iplot){i_unit});
                        
                        xlabel('Time (ms)');
                        xticklabels(xticks*1000); %convert in ms
                        ylabel('Spike waveform (�V)');
                        title([]);
                        setfig();
                        
                        stats{ipart}.(iplot).spikewaveform.halfwidth{i_unit}      = halfwidth;
                        stats{ipart}.(iplot).spikewaveform.peaktrough{i_unit}     = peaktrough;
                        stats{ipart}.(iplot).spikewaveform.troughpeak{i_unit}     = troughpeak;
                    else
                        axis off;
                    end
                end
            end
            
            %% baseline freq over time
            if hasevent
                subplot(7,4,[15 16]); hold;
                
                %Color darker for the last trials
                c_greys = 0.9 : -0.9 / size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo,1) : 0; 
                
                set_y_lim = []; %ONLY FOR PAUL : resize freq ylim with the max value at end
                for itrial = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo, 1)
                    color = [c_greys(itrial), c_greys(itrial), c_greys(itrial)];
                    x = statsovertime{ipart}.(baseline_name).time{i_unit}{itrial} - statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(end); %timelock to the end
                    plot(x, statsovertime{ipart}.(baseline_name).freq{i_unit}{itrial}, 'Color', color)
                    set_y_lim(itrial) = statsovertime{ipart}.(baseline_name).freq{i_unit}{itrial}(end);
                end
                
                ylabel('Spikerate (Hz)');
                setfig();
                ax = axis;
                try ylim([0, max(set_y_lim)]); end
                titlepos = title(sprintf('\n%s : %d trials, %d spikes',convertStringsToChars(baseline_name), size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo,1), size(cfg.SpikeTrials{ipart}.(baseline_name).trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                titlepos.Position(1) = ax(1);
                
                %% baseline CV2 over time
                subplot(7,4,[19 20]); hold;
                
                for itrial = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo, 1)
                    color = [c_greys(itrial), c_greys(itrial), c_greys(itrial)];
                    x = statsovertime{ipart}.(baseline_name).time{i_unit}{itrial} - statsovertime{ipart}.(baseline_name).time{i_unit}{itrial}(end); %timelock to the end
                    plot(x, statsovertime{ipart}.(baseline_name).cv2{i_unit}{itrial}, 'Color', color)
                end
                title([]);
                ylabel('Spikerate CV2');
                setfig();
                
                %% baseline CV2 over time, without bursts
                subplot(7,4,[23 24]); hold;
                
                for itrial = 1:size(cfg.SpikeTrials{ipart}.(baseline_name).trialinfo, 1)
                    color = [c_greys(itrial), c_greys(itrial), c_greys(itrial)];
                    x = statsovertime_withoutbursts{ipart}.(baseline_name).time{i_unit}{itrial} - statsovertime_withoutbursts{ipart}.(baseline_name).time{i_unit}{itrial}(end); %timelock to the end
                    plot(x, statsovertime_withoutbursts{ipart}.(baseline_name).cv2{i_unit}{itrial}, 'Color', color)
                end
                
                ylabel(sprintf('CV2 without \nbursts'));
                xlabel(sprintf('Time from end of %s',convertStringsToChars(baseline_name)));
                setfig();
              
            end %hasevent
            
            %% saveplot
            if ~(exist(cfg.imagesavedir)==7)
                mkdir(cfg.imagesavedir);
                ffigure
                printf('Create forlder %s',cfg.imagesavedir);
            end
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            fig.PaperType = 'A2';
            
            if hasevent
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}.(event_name).label{i_unit},'-spikestats_',convertStringsToChars(event_name),'_',convertStringsToChars(baseline_name),'.pdf']),'-r600');
                print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}.(event_name).label{i_unit},'-spikestats_',convertStringsToChars(event_name),'_',convertStringsToChars(baseline_name),'.png']),'-r600');
            else
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}.(baseline_name).label{i_unit},'-spikestats_',convertStringsToChars(baseline_name),'.pdf']),'-r600');
                print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}.(baseline_name).label{i_unit},'-spikestats_',convertStringsToChars(baseline_name),'.png']),'-r600');
            end
            
            close all
            
        end %i_unit
    end %event_name
end %ipart

%save stats output
save(fname, 'stats', '-v7.3');

end

% setfig = @() set(gca,'FontWeight','bold','TickDir','out');
function setfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set common features to all the subplots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
end
