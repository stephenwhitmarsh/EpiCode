function [stats] = spikeratestats_Events_Baseline(cfg,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [stats] = spikeratestats_Events_Baseline(cfg,force)
%
% Compare spike activity during short events versus resting periods.
% Can be used with only baseline data (no events) : set 
% cfg.spike.eventsname as 'no'.
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
% cfg.SpikeTrials{parts}{names} = spike data epoched in FieldTrip trial
%                                 data structure
% cfg.name{names}               = names of the different analysis
% cfg.spike.baselinename        = name of the analysis which correspond to
%                                 baseline/resting data
% cfg.imagesavedir              = where to save images
% cfg.circus.channel            = channel names which were analyzed by
%                                 Spyking-Circus
% 
% ### Optionnal inputs to add, if want to compute stats on events :
% cfg.spike.eventsname          = name(s) of the analysis which will be
%                                 processed as events. Can be 'no' to only 
%                                 use baseline data. Default = 'no'.
% cfg.spike.resamplefs{names}   = defines the frequency of bins in the
%                                 event's spikerate barplot
% cfg.spike.sdflatency{ievent}  = 'maxperiod' or array of 2 x values
% cfg.spike.sdftimwin{ievent}   = time window for convolution of spike density.
% cfg.stats.bltoi{names}        = baseline period for statistics
%                                 measurements of events
% cfg.stats.actoi{names}        = active period for statistics measurements
%                                 of events
% cfg.stats.alpha               = threshold for selecting significant
%                                 clusters. Default = 0.025
% cfg.stats.numrandomization    = nr of randomizations in the cluster-based
%                                 permutation test. Default = 1000.
%
% ### Optional cfg fields if want to compute stats on LFP and spike waveforms
% cfg.SpikeWaveforms{parts}{names}= spike wavefoms epoched in FieldTrip trial
%                                 data structure. Default = [].
% cfg.dataLFP{parts}{names}     = LFP data epoched in FieldTrip trial data
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
% stats{parts}.(cfg.name)       = MATLAB structure with one field per type
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
%    - discharge{i_unit}        = spike train descriptive stats : meanISI,
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
% Paul Baudin (paul.baudin@live.fr)
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
cfg.spike.eventsname        = ft_getopt(cfg.spike   , 'eventsname'          , 'no');
cfg.spike.ISIbins           = ft_getopt(cfg.spike   , 'ISIbins'             , [0:0.003:0.150]);
cfg.spike.RPV               = ft_getopt(cfg.spike   , 'RPV'                 , 0.003);
cfg.stats.alpha             = ft_getopt(cfg.stats   , 'alpha'               , 0.025);
cfg.stats.numrandomization  = ft_getopt(cfg.stats   , 'numrandomization'    , 1000);
cfg.stats.part_list         = ft_getopt(cfg.stats   , 'part_list'           , 'all');
cfg.circus.postfix          = ft_getopt(cfg.circus  , 'postfix'             , []);

if isempty(cfg.spike.eventsname)
    cfg.spike.eventsname = 'no';
end

if strcmp(cfg.stats.part_list,'all')
    cfg.stats.part_list = 1:size(cfg.SpikeTrials,2);
end

%% Find which labels in cfg.name correspond to events, baseline or LFP
event_indexes           = [];
LFP_indexes             = [];
baseline_index        = [];

%event indexes. if no indexes, hasevent = false.
hasevent = false;
if ~any(strcmp(cfg.spike.eventsname, 'no'))
    for i = 1 : size(cfg.spike.eventsname,2)
        for j = 1:size(cfg.name,2)
            if strcmp(cfg.spike.eventsname{i}, cfg.name{j})
                event_indexes = [event_indexes, j];
                hasevent = true;
            end
        end
    end
end

if ~hasevent, event_indexes = 1; end

% LFP indexes
if ~isempty(cfg.dataLFP)
    for i = 1 : size(event_indexes,2)
        hasLFP = false;
        for j = 1:size(cfg.LFP.name,2)
            if strcmp(cfg.name{event_indexes(i)}, cfg.LFP.name{j})
                hasLFP = true;
                LFP_indexes{i} = j;
            end
        end
        if ~hasLFP
            LFP_indexes{i} = [];
        end
    end
end

%baseline index
for i = 1 : size(cfg.name,2)
    if strcmp(cfg.name{i},cfg.spike.baselinename)
        baseline_index = i;
    end
end

%FIXME see if it is usefull to adapt this function to plot only events without baseline
if isempty(baseline_index)   , error('No baseline label found.\n')    ; end


for ipart = cfg.stats.part_list
    
    spike_Fs = cfg.SpikeTrials{ipart}{baseline_index}.hdr.Fs;
    
    for ievent = event_indexes
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% COMPUTE STATS %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % keep infos about the analysis
        stats{ipart}.(cfg.name{baseline_index}).method          = 'baseline';
        if hasevent, stats{ipart}.(cfg.name{ievent}).method     = 'event'; end
        stats{ipart}.(cfg.name{baseline_index}).label          = cfg.SpikeTrials{ipart}{baseline_index}.label;
        if hasevent, stats{ipart}.(cfg.name{ievent}).label     = cfg.SpikeTrials{ipart}{ievent}.label; end
        
        %% Compute ISI event, ISI baseline
        cfgtemp                                             = [];
        cfgtemp.outputunit                                  = 'spikecount';
        cfgtemp.bins                                        = cfg.spike.ISIbins;
        cfgtemp.param                                       = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials                                  = 'yes';
        stats{ipart}.(cfg.name{baseline_index}).isi         = ft_spike_isi(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});
        if hasevent, stats{ipart}.(cfg.name{ievent}).isi    = ft_spike_isi(cfgtemp,cfg.SpikeTrials{ipart}{ievent}); end
        
        
        %% Compute spike density for events
        if hasevent
            
            clear sdf_event sdf_baseline
            
            cfgtemp                         = [];
            cfgtemp.fsample                 = cfg.spike.resamplefs{ievent};   % sample at 1000 hz
            cfgtemp.keeptrials              = 'yes';
            cfgtemp.latency                 = cfg.spike.sdflatency{ievent};
            cfgtemp.timwin                  = cfg.spike.sdftimewin{ievent};
            [sdf_event]                     = ft_spikedensity(cfgtemp,cfg.SpikeTrials{ipart}{ievent});
            stats{ipart}.(cfg.name{ievent}).sdf = sdf_event;
            
            %% Compute cluster stats per unit for events
            
            % prepare dummy data with baseline value per trial for stats
            slim(1)                         = find(sdf_event.time > cfg.stats.bltoi{ievent}(1), 1, 'first');
            slim(2)                         = find(sdf_event.time < cfg.stats.bltoi{ievent}(2), 1, 'last');
            sdf_bl                          = sdf_event;
            sdf_bl.trial                    = ones(size(sdf_event.trial)) .* nanmean(sdf_event.trial(:,:,slim(1):slim(2)),3); % replace with mean
            
            for i_unit = 1 : size(cfg.SpikeTrials{ipart}{ievent}.label, 2)
                
                % statistics
                cfgtemp = [];
                cfgtemp.statistic                       = 'ft_statfun_depsamplesT';
                cfgtemp.alpha                           = cfg.stats.alpha;
                cfgtemp.clusteralpha                    = 0.01;
                cfgtemp.method                          = 'montecarlo';
                cfgtemp.computestat                     = 'yes';
                cfgtemp.correctm                        = 'cluster';
                cfgtemp.latency                         = [cfg.stats.bltoi{ievent}(2) sdf_event.time(end)]; % active period starts after baseline
                cfgtemp.ivar                            = 1;
                cfgtemp.uvar                            = 2;
                cfgtemp.design(1,:)                     = [ones(1,size(sdf_event.trial,1)) ones(1,size(sdf_bl.trial,1)) *2];
                cfgtemp.design(2,:)                     = [1 : size(sdf_event.trial,1) 1 : size(sdf_bl.trial,1)];
                cfgtemp.numrandomization                = cfg.stats.numrandomization;
                cfgtemp.channel                         = i_unit;
                stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit} = ft_timelockstatistics(cfgtemp,sdf_event,sdf_bl);
                
                % calculate baseline
                slim(1) = find(sdf_event.time > cfg.stats.bltoi{ievent}(1), 1, 'first');
                slim(2) = find(sdf_event.time < cfg.stats.bltoi{ievent}(2), 1, 'last');
                stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.avg        = nanmean(sdf_event.avg(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.var        = nanmean(sdf_event.var(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.dof        = nanmean(sdf_event.dof(i_unit,slim(1):slim(2)),2);
                stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.trialavg   = nanmean(sdf_event.trial(:,i_unit,slim(1):slim(2)),3);
            end
            
        end %hasevent
        
        %% %%%%%%%%%%%%%%%%%%%%
        %%%%%%%% PLOT %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        for i_unit = 1:size(cfg.SpikeTrials{ipart}{baseline_index}.label, 2)
            
            fig = figure;
            %title for the whole figure
            phy_channr      = cfg.SpikeTrials{ipart}{baseline_index}.template_maxchan(i_unit);
            nlx_channame    = cfg.circus.channel{phy_channr+1}; %+1 because starts at zero
            sgtitle(sprintf('Electrode %d (%s) : %s',phy_channr,nlx_channame, cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit}), 'Fontsize', 22, 'Interpreter', 'none', 'FontWeight', 'bold');
            
            %% Plot firing rate along all data
            if hasevent, subplot(7,4,1:3);hold; else, subplot(4,4,1:3);hold; end
            
            %for baseline
            cfgtemp                 = [];
            cfgtemp.spikechannel    = cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit};
            cfgtemp.cutlength       = 0; %no cut length with trialavg
            cfgtemp.method          = 'freq';
            cfgtemp.timelock        = 'no';
            cfgtemp.removeempty     = 'yes';
            cfgtemp.plot            = 'trialavg';
            cfgtemp.color           = 'b';              
            [stats{ipart}.(cfg.name{baseline_index}).freq_trialavg{i_unit}, legend_Bl] = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});

            if hasevent
                %for Event
                cfgtemp.plot            = 'scatter';
                cfgtemp.color           = 'r';
                [stats{ipart}.(cfg.name{ievent}).freq_trialavg{i_unit}, legend_Ev] = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{ievent});
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
            
            if hasevent, legend([legend_Bl, legend_Ev],cfg.name{baseline_index},cfg.name{ievent},'location','eastoutside'); end
            if ~hasevent, legend(legend_Bl,cfg.name{baseline_index},'location','eastoutside'); end
            
            %% Plot amplitude along all data
            if hasevent, subplot(7,4,5:7);hold; else, subplot(4,4,5:7);hold; end
            
            %for baseline
            cfgtemp.cutlength       = 0; %no cut length with trialavg
            cfgtemp.method          = 'amplitude';
            cfgtemp.plot            = 'trialavg';
            cfgtemp.color           = 'b';
            
            [stats{ipart}.(cfg.name{baseline_index}).amplitude_trialavg{i_unit}, legend_Bl] = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});

            if hasevent
                %for events
                cfgtemp.plot            = 'scatter';
                cfgtemp.color           = 'r';
                [stats{ipart}.(cfg.name{ievent}).amplitude_trialavg{i_unit}, legend_Ev] = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{ievent});
            end
            
            axis tight
            xticks(0:3600:ax(2));
            xticklabels(xticks/3600); %convert seconds to hours
            xlim([0 ax(2)]);
            xlabel('Time (hours)');
            set(gca, 'YColor', 'k');
            ylabel(sprintf('Mean amplitude \nof each trial'));
            setfig();
            
            if hasevent, legend([legend_Bl, legend_Ev],cfg.name{baseline_index},cfg.name{ievent},'location','eastoutside'); end
            if ~hasevent, legend(legend_Bl,cfg.name{baseline_index},'location','eastoutside'); end
            

            %% Template
            if hasevent, subplot(7,4,[4 8]); hold; else, subplot(8,4,[8 12]); hold; end
            
            tempsel = cfg.SpikeTrials{ipart}{baseline_index}.template{i_unit}(:,cfg.SpikeTrials{ipart}{baseline_index}.template_maxchan(i_unit)+1,:);%+1 because electrodes nr are zero-based
            temptime = ( (0 : size(cfg.SpikeTrials{ipart}{baseline_index}.template{i_unit},3) - 1) / spike_Fs )';
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
            cfgtemp.channame            = 'template';
            cfgtemp.mesurehalfwidth     = 'yes';
            cfgtemp.halfwidthmethod     = 'min'; 
            cfgtemp.mesurepeaktrough    = 'yes';
            cfgtemp.toiac               = 'all';
            cfgtemp.toibl               = []; %no need of bl if cfgtemp.halfwidthmethod = 'min'; 
            [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,template);
            
            title([]);
            ylabel('Template (µV)');
            xlabel('Time (ms)');
            xticklabels(xticks*1000);
            setfig();
            
            stats{ipart}.(cfg.name{baseline_index}).template.halfwidth{i_unit}    = halfwidth;
            stats{ipart}.(cfg.name{baseline_index}).template.peaktrough{i_unit}   = peaktrough;
            stats{ipart}.(cfg.name{baseline_index}).template.troughpeak{i_unit}   = troughpeak;
            stats{ipart}.(cfg.name{baseline_index}).template.time{i_unit}         = template.time;
            stats{ipart}.(cfg.name{baseline_index}).template.values{i_unit}       = template.trial;
            
            if hasevent
                stats{ipart}.(cfg.name{ievent}).template.halfwidth{i_unit}              = halfwidth;
                stats{ipart}.(cfg.name{ievent}).template.peaktrough{i_unit}             = peaktrough;
                stats{ipart}.(cfg.name{ievent}).template.troughpeak{i_unit}             = troughpeak;
                stats{ipart}.(cfg.name{ievent}).template.temptime_int{i_unit}           = template.time;
                stats{ipart}.(cfg.name{ievent}).template.tempval_int{i_unit}            = template.trial;
            end
            
            %% Events LFP (same chan as spikes). Ignored if no LFP data
            if hasevent
                
                subplot(7,4,[13 14]);hold;
                
                %check if has LFP data to plot
                plotLFP = false;
                if ~isempty(LFP_indexes)
                    iLFP = LFP_indexes{ievent};
                    
                    if ~isempty(iLFP)
                        if ~isempty(cfg.dataLFP{ipart}{iLFP})
                            plotLFP = true;
                            
                            %correct channel name for my data
                            if ~ismember (nlx_channame, cfg.dataLFP{ipart}{iLFP}.label)
                                nlx_channame = sprintf('%sLFP',nlx_channame);
                            end
                            if ~ismember (nlx_channame, cfg.dataLFP{ipart}{iLFP}.label)
                                warning('Cannot find LFP channel for %s in %s', cfg.SpikeTrials{ipart}{ievent}.label{i_unit}, cfg.name{ievent});
                                plotLFP = false;
                            end
                        end
                    end
                end
                
                if plotLFP
                    %correct baseline
                    cfgtemp = [];
                    cfgtemp.demean = 'yes';
                    cfgtemp.baselinewindow = cfg.stats.bltoi{ievent};
                    cfg.dataLFP{ipart}{iLFP} = ft_preprocessing(cfgtemp, cfg.dataLFP{ipart}{iLFP});
                    
                    %plot
                    cfgtemp                     = [];
                    cfgtemp.channame            = nlx_channame;
                    cfgtemp.plotstd             = 'yes';
                    cfgtemp.removeoutliers      = 'yes';
                    cfgtemp.toiplot             = cfg.epoch.toi{ievent};
                    cfgtemp.toibl               = cfg.stats.bltoi{ievent};
                    cfgtemp.toiac               = cfg.stats.actoi{ievent};
                    cfgtemp.mesurehalfwidth     = 'yes';
                    cfgtemp.halfwidthmethod     = 'bl';
                    cfgtemp.name                = cfg.LFP.name{ievent};
                    [hw_lfp, ~, ~] = plot_morpho(cfgtemp,cfg.dataLFP{ipart}{iLFP});
                    
                    stats{ipart}.(cfg.name{ievent}).LFP.channel{strcmp(cfg.dataLFP{ipart}{iLFP}.label,nlx_channame)'}     = nlx_channame;
                    stats{ipart}.(cfg.name{ievent}).LFP.halfwidth{strcmp(cfg.dataLFP{ipart}{iLFP}.label,nlx_channame)'}   = hw_lfp;
                    
                    xlabel([]);
                    titlepos = title(sprintf('\n%s : %d trials, %d spikes',cfg.name{ievent}, size(cfg.SpikeTrials{ipart}{ievent}.trialinfo,1), size(cfg.SpikeTrials{ipart}{ievent}.trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                    titlepos.Position(1) = cfg.epoch.toi{ievent}(1);
                    ylabel('LFP (µV)');
                    setfig();
                end
                
                if ~plotLFP
                    titlepos = title(sprintf('\n%s : %d trials, %d spikes',cfg.name{ievent}, size(cfg.SpikeTrials{ipart}{ievent}.trialinfo,1), size(cfg.SpikeTrials{ipart}{ievent}.trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                    titlepos.Position(1) = 0;
                    axis off;
                end
                
                %% Events raster plot
                subplot(7,4,[17 18]);hold;
                cfgtemp                 = [];
                cfgtemp.spikechannel    = i_unit;
                cfgtemp.latency         = [sdf_event.time(1) sdf_event.time(end)];
                cfgtemp.trialborders    = 'yes';
                ft_spike_plot_raster(cfgtemp,cfg.SpikeTrials{ipart}{ievent});
                xlabel([]);
                setfig();
                
                %% Events barplot and stats
                subplot(7,4,[21 22]); hold;
                
                %plot avg and std
                avg = bar(sdf_event.time,sdf_event.avg(i_unit,:),1,'edgecolor','none','facecolor',[0 0 0]);
                sd  = plot(sdf_event.time,sqrt(sdf_event.var(i_unit,:))+sdf_event.avg(i_unit,:),'-','LineWidth',2,'Color',[0.6 0.6 0.6]);
                
                ylabel('Spikerate (Hz)');
                axis tight
                ax = axis;
                
                baseline = stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.avg;
                idx_begin_stat = find(sdf_event.time > cfg.stats.bltoi{ievent}(2), 1, 'first');
                
                % plot positive clusters
                if isfield(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit},'posclusters')
                    for ipos = 1 : size(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.posclusters,2)
                        if stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.posclusters(ipos).prob < cfg.stats.alpha
                            sel = [];
                            sel = find(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.posclusterslabelmat == ipos);
                            if length(sel) == 1 %correct a bug in the plot if only one sample is positive
                                barwidth = stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(2)-stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(1); 
                            else
                                barwidth = 1;
                            end
                            bar(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),barwidth,'facecolor','g','edgecolor','g');
                            
                            % compute percentage of max
                            [~,max_idx_temp] = max(sdf_event.avg(i_unit,sel+idx_begin_stat-2));
                            max_idx = sel(max_idx_temp) + idx_begin_stat - 2;
                            x_max = sdf_event.time(max_idx);
                            y_max = sdf_event.avg(i_unit,max_idx);
                            percent_increase = (sdf_event.avg(i_unit,max_idx)-baseline) / baseline * 100;
                            %                             y = ax(4)*0.9;
                            %                             text(x_max,y,sprintf('max = %.1fHz \n(+%.1f%%)',y_max,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                            %                             plot(x_max,y_max,'*b');
                            
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.perc{ipos} = percent_increase;
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.x{ipos} = x_max;
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.y{ipos} = y_max;
                        end
                    end
                end
                
                %plot max
                max_freq = max(sdf_event.avg(i_unit,:));
                percent_increase = (max_freq-baseline) / baseline * 100;
                x = (cfg.stats.actoi{ievent}(1) + cfg.stats.actoi{ievent}(2))/2;
                y = ax(4)*0.9;
                text(x,y,sprintf('max = %.1fHz \n(+%.1f%%)',max_freq,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                stats{ipart}.(cfg.name{ievent}).maxfreq.max{i_unit} = max_freq;
                stats{ipart}.(cfg.name{ievent}).maxfreq.perc_increase{i_unit} = percent_increase;
                
                % plot negative clusters
                if isfield(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit},'negclusters')
                    for ineg = 1 : size(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusters,2)
                        if stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusters(ineg).prob < cfg.stats.alpha
                            sel = [];
                            sel = find(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusterslabelmat == ineg);
                            if length(sel) == 1 %correct a bug in the plot if only one sample is positive
                                barwidth = stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(2)-stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(1);
                            else
                                barwidth = 1;
                            end
                            bar(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),barwidth,'facecolor','r','edgecolor','r');
                            
                            % compute percentage
                            [~,min_idx_temp] = min(sdf_event.avg(i_unit,sel+idx_begin_stat-2));
                            min_idx = sel(min_idx_temp) + idx_begin_stat - 2;
                            x_min = sdf_event.time(min_idx);
                            y_min = sdf_event.avg(i_unit,min_idx);
                            percent_decrease = (sdf_event.avg(i_unit,min_idx)-baseline) / baseline * 100;
                            %                     %text(x_min,y_min+y_min/100,sprintf('%.1fHz (%.1f%%)\n',y_min,percent_decrease),'HorizontalAlignment','center','VerticalAlignment','middle');
                            %                     %plot(x_min,y_min,'*b');
                            
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.perc{ineg} = percent_decrease;
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.x{ineg} = x_min;
                            stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.y{ineg} = y_min;
                        end
                    end
                end
                
                % plot baseline patch
                x = [cfg.stats.bltoi{ievent}(1) cfg.stats.bltoi{ievent}(2) cfg.stats.bltoi{ievent}(2) cfg.stats.bltoi{ievent}(1)];
                y = [ax(3) ax(3) ax(4) ax(4)];
                patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
                
                % plot baseline
                x = [];
                x = [ax(1) ax(2)];
                plot(x,[baseline, baseline],':k');
                
                % plot baseline text
                x = (cfg.stats.bltoi{ievent}(1) + cfg.stats.bltoi{ievent}(2))/2;
                y = ax(4)*0.9;
                d = stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.bl.avg;
                text(x,y,sprintf('baseline = %.1f Hz',d),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                
                legend([avg, sd],'Average','SD');
                xlabel(sprintf('Time from %s (s)', cfg.name{ievent}),'FontSize',10);
                setfig();
                
            end %hasevent
            
            
            %% ISI for event and baseline
            
            if hasevent, plotlist = [ievent, baseline_index]; else, plotlist = baseline_index; end
            
            for iplot = plotlist
                
                if iplot == ievent && hasevent
                    subplot(7,4,25);hold;
                elseif iplot == baseline_index
                    if hasevent, subplot(7,4,27);hold; else, subplot(8,4,[21 22 25 26 29 30]);hold; end
                end
                
                bar(stats{ipart}.(cfg.name{iplot}).isi.time,stats{ipart}.(cfg.name{iplot}).isi.avg(i_unit,:), 'EdgeColor','none');
                
                RPV         = (length(find(stats{ipart}.(cfg.name{iplot}).isi.isi{i_unit} < cfg.spike.RPV)) / length(stats{ipart}.(cfg.name{iplot}).isi.isi{i_unit})) * 100;
                meanISI     = nanmean(stats{ipart}.(cfg.name{iplot}).isi.isi{i_unit})*1000;
                stdISI      = nanstd(stats{ipart}.(cfg.name{iplot}).isi.isi{i_unit})*1000;
                meanfreq    = 1/(meanISI/1000);
                
                %cv2
                arraytemp = stats{ipart}.(cfg.name{iplot}).isi.isi{i_unit};
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
                
                stats{ipart}.(cfg.name{iplot}).discharge.meanISI{i_unit}   = meanISI;
                stats{ipart}.(cfg.name{iplot}).discharge.stdISI{i_unit}    = stdISI;
                stats{ipart}.(cfg.name{iplot}).discharge.meanfreq{i_unit}  = meanfreq;
                stats{ipart}.(cfg.name{iplot}).discharge.RPV{i_unit}       = RPV;
                stats{ipart}.(cfg.name{iplot}).discharge.meancv2{i_unit}   = meancv2;
                stats{ipart}.(cfg.name{iplot}).discharge.stdcv2{i_unit}    = stdcv2;
            end
            
            if ~hasevent
                titlepos = title(sprintf('\n%s : %d spikes\n',cfg.name{baseline_index}, size(cfg.SpikeTrials{ipart}{baseline_index}.trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                titlepos.Position(1) = 0;
            end
            
            %% Spike Waveforms for events and baseline
            for iplot = plotlist
                
                if iplot == ievent && hasevent
                    subplot(7,4,26);hold;
                elseif iplot == baseline_index
                    if hasevent, subplot(7,4,28);hold; else, subplot(8,4,[23 24 27 28 31 32]);hold; end
                end
                
                if ~isempty(cfg.SpikeWaveforms)
                    if ~isempty(cfg.SpikeWaveforms{ipart}{iplot}{i_unit})
                        
                        cfgtemp                     = [];
                        cfgtemp.channame            = cfg.SpikeWaveforms{ipart}{iplot}{i_unit}.label{1};
                        cfgtemp.plotstd             = 'yes';
                        cfgtemp.removeoutliers      = 'yes'; %if big noise, impair seeing real data. Still present in avg and std.
                        cfgtemp.mesurehalfwidth     = 'yes';
                        cfgtemp.halfwidthmethod     = 'min';
                        cfgtemp.mesurepeaktrough    = 'yes';
                        cfgtemp.toibl               = []; %no need of bl if cfgtemp.halfwidthmethod     = 'min';
                        cfgtemp.toiac               = 'all';
                        cfgtemp.name                = cfg.name{iplot};
                        [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,cfg.SpikeWaveforms{ipart}{iplot}{i_unit});
                        
                        xlabel('Time (ms)');
                        xticklabels(xticks*1000); %convert in ms
                        ylabel('Spike waveform (µV)');
                        title([]);
                        setfig();
                        
                        stats{ipart}.(cfg.name{iplot}).spikewaveform.halfwidth{i_unit}      = halfwidth;
                        stats{ipart}.(cfg.name{iplot}).spikewaveform.peaktrough{i_unit}     = peaktrough;
                        stats{ipart}.(cfg.name{iplot}).spikewaveform.troughpeak{i_unit}     = troughpeak;
                    else
                        axis off;
                    end
                end
            end
            
            %% baseline freq over time
            if hasevent
                
                subplot(7,4,[15 16]); hold;
                
                cfgtemp                 = [];
                cfgtemp.spikechannel    = cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit};
                cfgtemp.cutlength       = 10; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
                cfgtemp.method          = 'freq';
                cfgtemp.timelock        = 'end';
                cfgtemp.removeempty     = 'yes';
                cfgtemp.removeoutlier   = 'yes';
                cfgtemp.plot            = 'movmean';               
                stats{ipart}.(cfg.name{baseline_index}).stats_over_time.freq{i_unit} = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});
                
                ylabel('Spikerate (Hz)');
                setfig();
                ax = axis;
                titlepos = title(sprintf('\n%s : %d trials, %d spikes',cfg.name{baseline_index}, size(cfg.SpikeTrials{ipart}{baseline_index}.trialinfo,1), size(cfg.SpikeTrials{ipart}{baseline_index}.trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                titlepos.Position(1) = ax(1);
                
                %% baseline CV2 over time
                subplot(7,4,[19 20]); hold;
                
                cfgtemp.method          = 'cv2';
                stats{ipart}.(cfg.name{baseline_index}).stats_over_time.cv2{i_unit} = spikestatsOverTime(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});
                
                title([]);
                ylabel('Spikerate CV2');
                setfig();
                
                %% baseline CV2 over time, without bursts
                subplot(7,4,[23 24]); hold;
                
                cfgtemp.method          = 'cv2';
                cfgtemp.removebursts    = 'yes';
                
                stats{ipart}.(cfg.name{baseline_index}).stats_over_time.cv2_withoutbursts{i_unit} = spikestatsOverTime(cfgtemp, cfg.SpikeTrials{ipart}{baseline_index});
                
                ylabel(sprintf('CV2 without \nbursts'));
                xlabel(sprintf('Time from end of %s',cfg.name{baseline_index}));
                setfig();
              
            end %hasevent
            
            %% saveplot
            if ~(exist(cfg.imagesavedir)==7)
                mkdir(cfg.imagesavedir);
                fprintf('Create forlder %s',cfg.imagesavedir);
            end
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            fig.PaperType = 'A2';
            
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}{ievent}.label{i_unit},'-spikestats_',cfg.name{ievent},'_',cfg.name{baseline_index},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}{ievent}.label{i_unit},'-spikestats_',cfg.name{ievent},'_',cfg.name{baseline_index},'.png']),'-r600');
            close all
            
        end %i_unit
    end %ievent
end %ipart

%save stats output
save(fname, 'stats');

end

% setfig = @() set(gca,'FontWeight','bold','TickDir','out');
function setfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set common features to all the subplots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
end
