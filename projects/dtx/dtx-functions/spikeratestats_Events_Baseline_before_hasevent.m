function [stats] = spikeratestats_Events_Baseline(cfg,SpikeTrials,SpikeWaveforms,dat_LFP,parts_to_read,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [stats] = spikeratestats_Events_Baseline(cfg,SpikeTrials,SpikeWaveforms,dat_LFP,parts_to_read,force)
% 
% Compare spike activity during short events versus resting periods.
% Can be used with only baseline data (no events) : leave
% cfg.spike.eventsname empty, or 'no'.
% 
% ### INPUT
% 
% # Necessary fields
% SpikeTrials{parts}{names}     = spike data epoched in FieldTrip trial 
%                                 data structure
% cfg.name{names}               = names of all the different analysis
% cfg.spike.eventsname          = name(s) of the analysis which will be
%                                 processed as events. Can be empty or 'no'
%                                 to only use baseline data
% cfg.spike.baselinename        = name of the analysis which correspond to
%                                 baseline/resting data
% cfg.spike.ISIbins             = bins limits to plot ISI histogram, in seconds
% cfg.spike.RPV                 = definition of refractory period
%                                 violation, in seconds
% cfg.spike.resamplefs{names}   = defines the frequency of bins in the
%                                 event's spikerate barplot
% cfg.stats.bltoi{names}        = baseline period for statistics measurements
% cfg.stats.actoi{names}        = active period for statistics measurements
% cfg.stats.alpha               = threshold for selecting significant clusters
% cfg.stats.numrandomization    = nr of randomizations in the cluster-based
%                                 permutation test
% cfg.prefix                    = prefix to output files
% cfg.datasavedir               = data directory of results
% cfg.imagesavedir              = where to save images
% cfg.circus.channel            = channel names which were analyzed by 
%                                 Spyking-Circus
% force                         = whether to redo analyses/plot or read 
%                                 previous saved stats file (true/false)
% 
% # Optional fields (corresponding plots/stats are empty if those data miss)
% SpikeWaveforms{parts}{names}  = spike wavefoms epoched in FieldTrip trial 
%                                 data structure
% dat_LFP{parts}{names}         = LFP data epoched in FieldTrip trial data 
%                                 structure
% cfg.LFP.name                  = name of the LFP analysis (to select which
%                                 one corresponds to cfg.spike.eventname)
% cfg.epoch.toi                 = time period to plot, for LFP data
%
% ### OUTPUT 
% stats{parts}.(cfg.name)       = MATLAB structure with one field per type 
%                                 of period analyzed, with infos computed for
%                                 each unit
%    - method                   = 'event' or 'baseline'
%    - isi                      = output from ft_spike_isi
%    - sdf                      = output from ft_spikedensity
%    - clusterstats             = output from ft_timelockanalysis, with
%                                 additional infos about max cluster and min 
%                                 cluster (position, value, percent of 
%                                 increase/decrease)
%    - freq_trialavg{i_unit}    = avg freq for each trial, per unit
%    - template{i_unit}         = measured halfwidth, peaktrough and
%                                 troughpeak for each template, in seconds
%    - LFP.channel /.halfwidth  = measured halfwidth on LFP events, for
%                                 each analyzed channel.
%    - discharge{i_unit}        = spike train descriptive stats : meanISI, 
%                                 stdISI, meanfreq, %RPV, meanCV2, stdCV2
%    - spikewaveform{i_unit}    = measured halfwidth, peaktrough and
%                                 troughpeak on waveforms for each unit
%    - stats_over_time          = for cv2 and freq, with and without time
%                                 normalization. For each : values, time, 
%                                 std, avg, per trial
% 
% ### Dependencies
% - Fieldtrip
% - To create the input structures : readSpykingCircus.m, readLFP.m, readSpikeWaveforms.m
% - Used in this script : plot_morpho.m, spikestatsOverTime.m
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikeratestats_events_baseline', cfg.circus.postfix, '.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed spike stats data\n');
    load(fname,'stats');
    return
end

%% Find which labels in cfg.name correspond to events, baseline or LFP
event_indexes           = [];
LFP_indexes             = [];
baseline_index        = [];

%event indexes. if no indexes, hasevent = false.
hasevent = false;
if isfield(cfg.spike, 'eventsname')
    if ~any(strcmp(cfg.spike.eventsname, 'no')) || ~isempty(cfg.spike.eventsname)
        for i = 1 : size(cfg.spike.eventsname,2)
            for j = 1:size(cfg.name,2)
                if strcmp(cfg.spike.eventsname{i}, cfg.name{j})
                    event_indexes = [event_indexes, j];
                    hasevent = true;
                end
            end
        end
    end
end
if ~hasevent, event_indexes = 1; end

% LFP indexes
if ~isempty(dat_LFP)
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

%FIXME see if it is usefull to adapt this function to plot only event without baseline
if isempty(event_indexes)   , error('No events labels found.\n')    ; end

if strcmp(parts_to_read,'all'), parts_to_read = 1:size(cfg.directorylist,2); end

for ipart = parts_to_read
    
    spike_Fs = SpikeTrials{ipart}{baseline_index}.hdr.Fs;

    for ievent = event_indexes
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% COMPUTE STATS %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        stats{ipart}.(cfg.name{ievent}).method              = 'event';
        stats{ipart}.(cfg.name{baseline_index}).method    = 'baseline';    
        
        %% Compute ISI event, ISI baseline
        cfgtemp                                             = [];
        cfgtemp.outputunit                                  = 'proportion';
        cfgtemp.bins                                        = cfg.spike.ISIbins;   
        cfgtemp.param                                       = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials                                  = 'yes';
        stats{ipart}.(cfg.name{ievent}).isi                 = ft_spike_isi(cfgtemp,SpikeTrials{ipart}{ievent});
        stats{ipart}.(cfg.name{baseline_index}).isi       = ft_spike_isi(cfgtemp,SpikeTrials{ipart}{baseline_index});
        
       
         
        %% Compute spike density for events 
        %trialinfo: column 3 StartSample 4 EndSample
        %voir si ajuster timwin
        clear sdf_event sdf_baseline
        
        cfgtemp                         = [];
        cfgtemp.fsample                 = cfg.spike.resamplefs{ievent};   % sample at 1000 hz
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.latency                 = 'maxperiod';%[cfg.epoch.toi{ievent}(1), cfg.epoch.toi{ievent}(2)];
        cfgtemp.timwin                  = [-1/cfg.spike.resamplefs{ievent}/2 1/cfg.spike.resamplefs{ievent}/2];%timewin same length as period
        [sdf_event]                     = ft_spikedensity(cfgtemp,SpikeTrials{ipart}{ievent});
        stats{ipart}.(cfg.name{ievent}).sdf = sdf_event;
        
        %% Compute cluster stats per unit for events

        % prepare dummy data with baseline value per trial for stats
        slim(1)                         = find(sdf_event.time > cfg.stats.bltoi{ievent}(1), 1, 'first');
        slim(2)                         = find(sdf_event.time < cfg.stats.bltoi{ievent}(2), 1, 'last');
        sdf_bl                          = sdf_event;
        sdf_bl.trial                    = ones(size(sdf_event.trial)) .* nanmean(sdf_event.trial(:,:,slim(1):slim(2)),3); % replace with mean
        
         for i_unit = 1 : size(SpikeTrials{ipart}{ievent}.label, 2) 
             
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
         

         %% %%%%%%%%%%%%%%%%%%%%
         %%%%%%%% PLOT %%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%
         
         for i_unit = 1:size(SpikeTrials{ipart}{ievent}.label, 2)
         
             
         fig = figure;
         
         phy_channr = SpikeTrials{ipart}{ievent}.template_maxchan(i_unit);
         nlx_channame = cfg.circus.channel{phy_channr+1}; %+1 because starts at zero
         sgtitle(sprintf('Electrode %d (%s) : %s',phy_channr,nlx_channame, SpikeTrials{ipart}{ievent}.label{i_unit}), 'Fontsize', 22, 'Interpreter', 'none', 'FontWeight', 'bold');
         
        %% Plot firing rate along all data
        subplot(7,4,1:3);hold;
        
        %for baseline
        for itrial = 1:size(SpikeTrials{ipart}{baseline_index}.trialinfo,1)
            yyaxis left
            clear spike_indx
            spike_indx              = SpikeTrials{ipart}{baseline_index}.trial{i_unit} == SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
            trial_begin             = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
            trial_freq_avg(itrial)  = 1 / nanmean(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit}(spike_indx));
            baseline              = plot([trial_begin trial_end], log10([trial_freq_avg(itrial) trial_freq_avg(itrial)]),'-b','LineWidth',2);
            
            yyaxis right
            if log10(trial_freq_avg(itrial)) ~= Inf && log10(trial_freq_avg(itrial)) ~= -Inf && ~isnan(trial_freq_avg(itrial))
                hasdata = plot([trial_begin trial_end],[95 95],'-k','LineWidth', 2);
            end
        end
        stats{ipart}.(cfg.name{baseline_index}).freq_trialavg{i_unit} = trial_freq_avg;
        
        clear trial_freq_avg
        %for event
        for itrial = 1:size(SpikeTrials{ipart}{ievent}.trialinfo,1)
            yyaxis left
            clear spike_indx
            spike_indx              = SpikeTrials{ipart}{ievent}.trial{i_unit} == SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
            trial_begin             = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 4) / spike_Fs;
            trial_freq_avg(itrial)  = 1 / nanmean(stats{ipart}.(cfg.name{ievent}).isi.isi{i_unit}(spike_indx));
            event                   = scatter(trial_begin, log10(trial_freq_avg(itrial)),7,'or','filled');%, 'filled');
            yyaxis right
            if log10(trial_freq_avg(itrial)) ~= Inf && log10(trial_freq_avg(itrial)) ~= -Inf && ~isnan(trial_freq_avg(itrial))
                hasdata = plot([trial_begin trial_end],[95 95],'-k','LineWidth', 2);
            end
        end
        stats{ipart}.(cfg.name{ievent}).freq_trialavg{i_unit} = trial_freq_avg;
        
        yyaxis left 
        axis tight
        ax = axis;
        xticks(0:3600:ax(2));
        xticklabels(xticks/3600); %convert seconds to hours
        xlim([0 Inf]);
        xlabel('Time (hours)');
        set(gca, 'YColor', 'k');
        ylabel(sprintf('Mean firing rate \nof each trial (log10)'));
        yticklabels(10.^yticks);
        setfig();
        
        yyaxis right
        ylim([0 100]);
        set(gca, 'YColor', 'none');
        
        if exist('hasdata', 'var')
            legend([baseline, event, hasdata],cfg.name{baseline_index},cfg.name{ievent},'Has data','location','eastoutside');
        end
        
          %% Plot amplitude along all data %FIXME : do it after solving the issue of selecting amplitudes in trials
        subplot(7,4,5:7);hold;
        
        clear baseline event hasdata
        
        %for baseline
        for itrial = 1:size(SpikeTrials{ipart}{baseline_index}.trialinfo,1)
            clear spike_indx
            spike_indx              = SpikeTrials{ipart}{baseline_index}.trial{i_unit} == SpikeTrials{ipart}{baseline_index}.trialinfo(itrial,7);
            trial_begin             = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
            trial_amp_avg(itrial)   = nanmean(SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
            trial_amp_std(itrial)   = nanstd(SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
            baseline              = plot([trial_begin trial_end], [trial_amp_avg(itrial) trial_amp_avg(itrial)],'-b','LineWidth',2);
        end
        stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.avg{i_unit} = trial_amp_avg;
        stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.std{i_unit} = trial_amp_std;

        clear trial_amp_avg trial_amp_avg
        %for event
        for itrial = 1:size(SpikeTrials{ipart}{ievent}.trialinfo,1)
            clear spike_indx
            spike_indx              = SpikeTrials{ipart}{ievent}.trial{i_unit} == SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
            trial_begin             = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = SpikeTrials{ipart}{ievent}.trialinfo(itrial, 4) / spike_Fs;
            trial_amp_avg(itrial)   = nanmean(SpikeTrials{ipart}{ievent}.amplitude{i_unit}(spike_indx));
            trial_amp_std(itrial)   = nanstd(SpikeTrials{ipart}{ievent}.amplitude{i_unit}(spike_indx));
            event                   = scatter(trial_begin, trial_amp_avg(itrial),7,'or','filled');

        end
        stats{ipart}.(cfg.name{ievent}).amp_trialavg.avg{i_unit} = trial_amp_avg;
        stats{ipart}.(cfg.name{ievent}).amp_trialavg.std{i_unit} = trial_amp_std;

        axis tight
        ax = axis;
        plot([0 ax(2)], [1 1], '--k'); %amplitude 1 = template
        xticks(0:3600:ax(2));
        xticklabels(xticks/3600); %convert seconds to hours
        xlim([0 Inf]);
        ylim([0.5 1.5]);
        xlabel('Time (hours)');
        set(gca, 'YColor', 'k');
        ylabel(sprintf('Mean spike amplitude \nof each trial'));
        legend([baseline, event],cfg.name{baseline_index},cfg.name{ievent},'location','eastoutside');
        setfig();
        
        yyaxis right %to set the same size of plot as firing rata
        set(gca, 'YColor', 'none');

        
        %% Template 
        subplot(7,4,[4 8]); hold;
        
        tempsel = SpikeTrials{ipart}{ievent}.template{i_unit}(:,SpikeTrials{ipart}{ievent}.template_maxchan(i_unit)+1,:);%+1 because electrodes nr are zero-based
        temptime = ( (0 : size(SpikeTrials{ipart}{ievent}.template{i_unit},3) - 1) / spike_Fs )';
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
        cfgtemp.plotstd             = 'no';
        cfgtemp.removeoutliers      = 'no';
        cfgtemp.mesurehalfwidth     = 'yes';
        cfgtemp.halfwidthmethod     = 'min'; %or bl
        cfgtemp.mesurepeaktrough    = 'yes';
        cfgtemp.toiplot             = 'all';
        cfgtemp.toiac               = 'all';
        cfgtemp.toibl               = [];
        cfgtemp.name                = [];
        cfgtemp.saveplot            = 'no';
        cfgtemp.imagesavedir        = [];
        cfgtemp.prefix              = [];
        [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,template);
        
        title([]);
        ylabel('Template');
        xlabel('Time (ms)');
        xticklabels(xticks*1000);
        setfig();
        
        stats{ipart}.(cfg.name{ievent}).template.halfwidth{i_unit}              = halfwidth;
        stats{ipart}.(cfg.name{ievent}).template.peaktrough{i_unit}             = peaktrough;
        stats{ipart}.(cfg.name{ievent}).template.troughpeak{i_unit}             = troughpeak;        
        stats{ipart}.(cfg.name{baseline_index}).template.halfwidth{i_unit}    = halfwidth;
        stats{ipart}.(cfg.name{baseline_index}).template.peaktrough{i_unit}   = peaktrough;
        stats{ipart}.(cfg.name{baseline_index}).template.troughpeak{i_unit}   = troughpeak;
        
        %% Events LFP (same chan as spikes). Ignored if no LFP data
        subplot(7,4,[13 14]);hold;
        plotLFP = false;
        
        if ~isempty(LFP_indexes)
            iLFP = LFP_indexes{ievent};
            
            if ~isempty(iLFP)
                if ~isempty(dat_LFP{ipart}{iLFP})
                    plotLFP = true;
                    
                    %correct channel name for my data
                    if ~ismember (nlx_channame, dat_LFP{ipart}{iLFP}.label)
                        nlx_channame = sprintf('%sLFP',nlx_channame);
                    end
                    if ~ismember (nlx_channame, dat_LFP{ipart}{iLFP}.label)
                        warning('Cannot find LFP channel for %s in %s', SpikeTrials{ipart}{ievent}.label{i_unit}, cfg.name{ievent});
                        plotLFP = false;
                    end
                    
                    if plotLFP
                        %correct baseline
                        cfgtemp = [];
                        cfgtemp.demean = 'yes';
                        cfgtemp.baselinewindow = cfg.stats.bltoi{ievent};
                        dat_LFP{ipart}{iLFP} = ft_preprocessing(cfgtemp, dat_LFP{ipart}{iLFP});
                        
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
                        cfgtemp.mesurepeaktrough    = 'no';
                        cfgtemp.name                = cfg.LFP.name{ievent};
                        cfgtemp.saveplot            = 'no';
                        cfgtemp.imagesavedir        = cfg.imagesavedir;
                        cfgtemp.prefix              = cfg.prefix;
                        [hw_lfp, ~, ~] = plot_morpho(cfgtemp,dat_LFP{ipart}{iLFP});
                        
                        stats{ipart}.(cfg.name{ievent}).LFP.channel{strcmp(dat_LFP{ipart}{iLFP}.label,nlx_channame)'}     = nlx_channame;
                        stats{ipart}.(cfg.name{ievent}).LFP.halfwidth{strcmp(dat_LFP{ipart}{iLFP}.label,nlx_channame)'}   = hw_lfp;
                        
                        xlabel([]);
                        titlepos = title(sprintf('\n%s : %d trials',cfg.name{ievent}, size(SpikeTrials{ipart}{ievent}.trialinfo,1)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
                        titlepos.Position(1) = cfg.epoch.toi{ievent}(1);
                        ylabel('LFP (µV)');
                        setfig();
                    end
                end
            end
        end
        
        if ~plotLFP 
            titlepos = title(sprintf('\n%s : %d trials',cfg.name{ievent}, size(SpikeTrials{ipart}{ievent}.trialinfo,1)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
            titlepos.Position(1) = 0;
            axis off;
        end
            
        
         %% Events raster plot
        subplot(7,4,[17 18]);hold;
        cfgtemp                 = [];
        cfgtemp.spikechannel    = i_unit;
        cfgtemp.latency         = [sdf_event.time(1) sdf_event.time(end)]; 
        cfgtemp.trialborders    = 'yes';
        ft_spike_plot_raster(cfgtemp,SpikeTrials{ipart}{ievent});
        xlabel([]);
        setfig();
        
        %% Events barplot and stats
        subplot(7,4,[21 22]); hold;
        
        bar(sdf_event.time,sdf_event.avg(i_unit,:),1,'edgecolor','none','facecolor',[0 0 0]);
        plot(sdf_event.time,sqrt(sdf_event.var(i_unit,:))+sdf_event.avg(i_unit,:),'-','LineWidth',2,'Color',[0.6 0.6 0.6]);
        
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
                    bar(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),1,'facecolor','g','edgecolor','g');
                    
                    % plot percentage 
                    [~,max_idx_temp] = max(sdf_event.avg(i_unit,sel+idx_begin_stat-2));
                    max_idx = sel(max_idx_temp) + idx_begin_stat - 2;
                    x_max = sdf_event.time(max_idx);
                    y_max = sdf_event.avg(i_unit,max_idx);
                    percent_increase = (sdf_event.avg(i_unit,max_idx)-baseline) / baseline * 100;
                    y = ax(4)*0.9;
                    text(x_max,y,sprintf('max = %.1fHz \n(+%.1f%%)',y_max,percent_increase),'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight', 'bold', 'FontSize',10);
                    plot(x_max,y_max,'*b');
                    
                    stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.perc{ipos} = percent_increase;
                    stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.x{ipos} = x_max;
                    stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.maxcluster.y{ipos} = y_max;
                end
            end
        end
                
        % plot negative clusters
        if isfield(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit},'negclusters')
            for ineg = 1 : size(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusters,2)
                if stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusters(ineg).prob < cfg.stats.alpha
                    sel = [];
                    sel = find(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.negclusterslabelmat == ineg);
                    bar(stats{ipart}.(cfg.name{ievent}).clusterstat{i_unit}.time(sel),sdf_event.avg(i_unit,sel+idx_begin_stat-2),1,'facecolor','r','edgecolor','r');
                    
                    % plot percentage 
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
        
        legend('Average','SD');
        xlabel(sprintf('Time from %s (s)', cfg.name{ievent}),'FontSize',10);
        setfig();

       
        
        %% ISI for event and baseline
        for iplot = [ievent, baseline_index]
            
            if iplot == ievent
                subplot(7,4,25);hold;
            elseif iplot == baseline_index
                subplot(7,4,27);hold;
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
            
            yticklabels(yticks*100); %convert in percent
            xticklabels(xticks*1000); %convert in ms
            ylabel('Proportion (%)');
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
    
        %% Spike Waveforms for events and baseline
        for iplot = [ievent, baseline_index]
            
            if iplot == ievent
                subplot(7,4,26);hold;
            elseif iplot == baseline_index
                subplot(7,4,28);hold;
            end
            
            if ~isempty(SpikeWaveforms)
                if ~isempty(SpikeWaveforms{ipart}{iplot}{i_unit})
                    
                    cfgtemp                     = [];
                    cfgtemp.channame            = SpikeWaveforms{ipart}{iplot}{i_unit}.label{1};
                    cfgtemp.plotstd             = 'yes';
                    cfgtemp.toiplot             = 'all';
                    cfgtemp.removeoutliers      = 'yes';
                    cfgtemp.mesurehalfwidth     = 'yes';
                    cfgtemp.halfwidthmethod     = 'min';
                    cfgtemp.mesurepeaktrough    = 'yes';
                    cfgtemp.toibl               = [];
                    cfgtemp.toiac               = 'all';
                    cfgtemp.name                = cfg.name{iplot};
                    cfgtemp.saveplot            = 'no';
                    cfgtemp.imagesavedir        = [];
                    cfgtemp.prefix              = [];
                    [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,SpikeWaveforms{ipart}{iplot}{i_unit});
                    
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
        
        %% baseline freq over time, normalized
        subplot(7,4,15); hold;
        
        cfgtemp                 = [];
        cfgtemp.spikechannel    = SpikeTrials{ipart}{baseline_index}.label{i_unit};
        cfgtemp.normtime        = 'yes';
        cfgtemp.normvalues      = 'no';
        cfgtemp.cutlength       = 0.1; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
        cfgtemp.method          = 'freq';
        cfgtemp.removeempty     = 'yes';
        cfgtemp.removeoutlier   = 'yes';
        cfgtemp.saveplot        = 'no';
        cfgtemp.name            = [];
        cfgtemp.prefix          = [];
        cfgtemp.imagesavedir    = [];
        
        stats{ipart}.(cfg.name{baseline_index}).stats_over_time.freq_timenorm{i_unit} = spikestatsOverTime(cfgtemp, SpikeTrials{ipart}{baseline_index});
        
        titlepos = title(sprintf('\n%s : %d trials',cfg.name{baseline_index}, size(SpikeTrials{ipart}{baseline_index}.trialinfo,1)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
        titlepos.Position(1) = cfgtemp.cutlength;
        ylabel('Spikerate, time norm');
        setfig();

        %% baseline CV2 over time, normalized
        subplot(7,4,16); hold;
        
        cfgtemp                 = [];
        cfgtemp.spikechannel    = SpikeTrials{ipart}{baseline_index}.label{i_unit};
        cfgtemp.normtime        = 'yes';
        cfgtemp.normvalues      = 'no';
        cfgtemp.cutlength       = 0.1; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
        cfgtemp.method          = 'cv2';
        cfgtemp.removeempty     = 'yes';
        cfgtemp.removeoutlier   = 'yes';
        cfgtemp.saveplot        = 'no';
        cfgtemp.name            = [];
        cfgtemp.prefix          = [];
        cfgtemp.imagesavedir    = [];

        stats{ipart}.(cfg.name{baseline_index}).stats_over_time.cv2_timenorm{i_unit} = spikestatsOverTime(cfgtemp,SpikeTrials{ipart}{baseline_index});
        
        title([]);
        ylabel('CV2, time norm');
        setfig();
                        
        %% baseline freq over time, non norm
        subplot(7,4,[19 20]); hold;
        
        cfgtemp                 = [];
        cfgtemp.spikechannel    = SpikeTrials{ipart}{baseline_index}.label{i_unit};
        cfgtemp.normtime        = 'no';
        cfgtemp.normvalues      = 'no';
        cfgtemp.cutlength       = 5; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
        cfgtemp.method          = 'freq';
        cfgtemp.removeempty     = 'yes';
        cfgtemp.removeoutlier   = 'yes';
        cfgtemp.saveplot        = 'no';
        cfgtemp.name            = [];
        cfgtemp.prefix          = [];
        cfgtemp.imagesavedir    = [];
        
        stats{ipart}.(cfg.name{baseline_index}).stats_over_time.freq{i_unit} = spikestatsOverTime(cfgtemp,SpikeTrials{ipart}{baseline_index});
        
        ylabel('Spikerate (Hz)');
        setfig();        
        
        %% baseline CV2 over time, non norm
        subplot(7,4,[23 24]); hold;
        
        cfgtemp                 = [];
        cfgtemp.spikechannel    = SpikeTrials{ipart}{baseline_index}.label{i_unit};
        cfgtemp.normtime        = 'no';
        cfgtemp.normvalues      = 'no';
        cfgtemp.cutlength       = 5; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
        cfgtemp.method          = 'cv2';
        cfgtemp.removeempty     = 'yes';
        cfgtemp.removeoutlier   = 'yes';
        cfgtemp.saveplot        = 'no';
        cfgtemp.name            = [];
        cfgtemp.prefix          = [];
        cfgtemp.imagesavedir    = [];
        
        stats{ipart}.(cfg.name{baseline_index}).stats_over_time.cv2{i_unit} = spikestatsOverTime(cfgtemp,SpikeTrials{ipart}{baseline_index});
        
        title([]);
        ylabel('Spikerate CV2');
        xlabel('Time from begining of trial (s)');
        setfig();
        
        %% saveplot
        if ~(exist(cfg.imagesavedir)==7)
            mkdir(cfg.imagesavedir);
            fprintf('Create forlder %s',cfg.imagesavedir);
        end
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        fig.PaperType = 'A2';
        
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',SpikeTrials{ipart}{ievent}.label{i_unit},'-spikestats_',cfg.name{ievent},'_',cfg.name{baseline_index},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',SpikeTrials{ipart}{ievent}.label{i_unit},'-spikestats_',cfg.name{ievent},'_',cfg.name{baseline_index},'.png']),'-r600');
        close all
        

        
         end %i_unit
    end %ievent
end %ipart

%save stats output
save(fname, 'stats');

end

function setfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set common features to all the subplots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
end
