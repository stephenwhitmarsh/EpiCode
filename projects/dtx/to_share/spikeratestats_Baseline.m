function [stats] = spikeratestats_Baseline(cfg,force)

% [stats] = spikeratestats_Baseline(cfg,force)
%
% Compute general stas on spike activity along all some data period. 
% Input spike data must be cut into trials to compute avg freq and
% amplitude over trials. 
% Use readSpikeTrials_continuous, and removetrialsMuseMarkers to
% respectively cut the data and remove trials which intersect BAD markers.
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
% ### Optional cfg fields if want to compute stats on spike waveforms
% cfg.SpikeWaveforms{parts}{names}= spike wavefoms epoched in FieldTrip trial
%                                 data structure. Default = [].
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
%    - label                    = label of each unit
%    - isi                      = output from ft_spike_isi
%    - freq_trialavg{i_unit}    = avg freq for each trial, per unit
%    - amplitude_trialavg{i_unit}= avg amplitude for each trial, per unit
%    - template{i_unit}         = measured halfwidth, peaktrough and
%                                 troughpeak for each template, in seconds
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


% load precomputed stats if required
fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikeratestats_baseline', cfg.circus.postfix, '.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed spike stats data\n');
    load(fname,'stats');
    return
end

% get default cfg fields
cfg.spike                   = ft_getopt(cfg         , 'spike'              	, []);
cfg.circus                  = ft_getopt(cfg         , 'circus'           	, []);
cfg.SpikeWaveforms          = ft_getopt(cfg         , 'SpikeWaveforms'    	, []);
cfg.spike.ISIbins           = ft_getopt(cfg.spike   , 'ISIbins'             , [0:0.003:0.150]);
cfg.stats.part_list         = ft_getopt(cfg.stats   , 'part_list'           , 'all');
cfg.circus.postfix          = ft_getopt(cfg.circus  , 'postfix'             , []);

if strcmp(cfg.stats.part_list,'all')
    cfg.stats.part_list = 1:size(cfg.SpikeTrials,2);
end

%Find which label in cfg.name correspond baseline
baseline_index        = [];
for i = 1 : size(cfg.name,2)
    if strcmp(cfg.name{i},cfg.spike.baselinename)
        baseline_index = i;
    end
end
if isempty(baseline_index)   , error('No baseline label found.\n')    ; end

%% compute stats over time
cfg.statstime.label_list    = baseline_index;
statsovertime               = spikestatsOverTime(cfg, cfg.SpikeTrials,true);


%% go trhough each part and labels
for ipart = cfg.stats.part_list
    
    spike_Fs = cfg.SpikeTrials{ipart}{baseline_index}.hdr.Fs;
    
    
    stats{ipart}.(cfg.name{baseline_index}).label          = cfg.SpikeTrials{ipart}{baseline_index}.label;
    
    %% Compute ISI baseline
    cfgtemp                                             = [];
    cfgtemp.outputunit                                  = 'spikecount';
    cfgtemp.bins                                        = cfg.spike.ISIbins;
    cfgtemp.param                                       = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
    cfgtemp.keeptrials                                  = 'yes';
    stats{ipart}.(cfg.name{baseline_index}).isi         = ft_spike_isi(cfgtemp,cfg.SpikeTrials{ipart}{baseline_index});
    
    %plot data for each unit
    for i_unit = 1:size(cfg.SpikeTrials{ipart}{baseline_index}.label, 2)
        
        fig = figure;
        %title for the whole figure
        phy_channr      = cfg.SpikeTrials{ipart}{baseline_index}.template_maxchan(i_unit);
        nlx_channame    = cfg.circus.channel{phy_channr+1}; %+1 because starts at zero
        sgtitle(sprintf('Electrode %d (%s) : %s',phy_channr,nlx_channame, cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit}), 'Fontsize', 22, 'Interpreter', 'none', 'FontWeight', 'bold');
        
        %% Plot firing rate along all data
        subplot(4,4,1:3);hold;
        
        for itrial = 1:size(cfg.SpikeTrials{ipart}{baseline_index}.trialinfo,1)
            yyaxis left
            clear spike_indx
            spike_indx              = cfg.SpikeTrials{ipart}{baseline_index}.trial{i_unit} == cfg.SpikeTrials{ipart}{ievent}.trialinfo(itrial,7);
            trial_begin             = cfg.SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = cfg.SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
            trial_freq_avg(itrial)  = 1 / nanmean(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit}(spike_indx));
            legend_Bl               = plot([trial_begin trial_end], log10([trial_freq_avg(itrial) trial_freq_avg(itrial)]),'-b','LineWidth',2);
        end
        stats{ipart}.(cfg.name{baseline_index}).freq_trialavg{i_unit} = trial_freq_avg;
        
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
        
        legend(legend_Bl,cfg.name{baseline_index},'location','eastoutside');
        
        %% Plot amplitude along all data
        subplot(4,4,5:7);hold;
        
        %for baseline
        for itrial = 1:size(cfg.SpikeTrials{ipart}{baseline_index}.trialinfo,1)
            clear spike_indx
            spike_indx              = cfg.SpikeTrials{ipart}{baseline_index}.trial{i_unit} == cfg.SpikeTrials{ipart}{baseline_index}.trialinfo(itrial,7);
            trial_begin             = cfg.SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 3) / spike_Fs;
            trial_end               = cfg.SpikeTrials{ipart}{baseline_index}.trialinfo(itrial, 4) / spike_Fs;
            trial_amp_avg(itrial)   = nanmean(cfg.SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
            trial_amp_std(itrial)   = nanstd(cfg.SpikeTrials{ipart}{baseline_index}.amplitude{i_unit}(spike_indx));
            legend_Bl               = plot([trial_begin trial_end], [trial_amp_avg(itrial) trial_amp_avg(itrial)],'-b','LineWidth',2);
        end
        stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.avg{i_unit} = trial_amp_avg;
        stats{ipart}.(cfg.name{baseline_index}).amp_trialavg.std{i_unit} = trial_amp_std;
        
        axis tight
        xticks(0:3600:ax(2));
        xticklabels(xticks/3600); %convert seconds to hours
        xlim([0 ax(2)]);
        xlabel('Time (hours)');
        set(gca, 'YColor', 'k');
        ylabel(sprintf('Mean amplitude \nof each trial'));
        setfig();
        
        legend(legend_Bl,cfg.name{baseline_index},'location','eastoutside');
        
        
        %% Template
        subplot(8,4,[8 12]); hold;
        
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
        cfgtemp.morpho.channame            = 'template';
        cfgtemp.morpho.mesurehalfwidth     = 'yes';
        cfgtemp.morpho.blmethod     = 'min';
        cfgtemp.morpho.mesurepeaktrough    = 'yes';
        cfgtemp.morpho.toiac               = 'all';
        cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.blmethod = 'min';
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
        
        
        %% ISI
        
        
        subplot(8,4,[21 22 25 26 29 30]);hold;
        
        bar(stats{ipart}.(cfg.name{baseline_index}).isi.time,stats{ipart}.(cfg.name{baseline_index}).isi.avg(i_unit,:), 'EdgeColor','none');
        
        RPV         = (length(find(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit} < cfg.spike.RPV)) / length(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit})) * 100;
        meanISI     = nanmean(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit})*1000;
        stdISI      = nanstd(stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit})*1000;
        meanfreq    = 1/(meanISI/1000);
        
        %cv2
        arraytemp = stats{ipart}.(cfg.name{baseline_index}).isi.isi{i_unit};
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
        
        stats{ipart}.(cfg.name{baseline_index}).firing.meanISI{i_unit}   = meanISI;
        stats{ipart}.(cfg.name{baseline_index}).firing.stdISI{i_unit}    = stdISI;
        stats{ipart}.(cfg.name{baseline_index}).firing.meanfreq{i_unit}  = meanfreq;
        stats{ipart}.(cfg.name{baseline_index}).firing.RPV{i_unit}       = RPV;
        stats{ipart}.(cfg.name{baseline_index}).firing.meancv2{i_unit}   = meancv2;
        stats{ipart}.(cfg.name{baseline_index}).firing.stdcv2{i_unit}    = stdcv2;
        
        titlepos = title(sprintf('\n%s : %d spikes\n',cfg.name{baseline_index}, size(cfg.SpikeTrials{ipart}{baseline_index}.trial{i_unit},2)),'Fontsize',22,'Interpreter','none','HorizontalAlignment','left');
        titlepos.Position(1) = 0;
        
        %% Spike Waveforms for baseline
        subplot(8,4,[23 24 27 28 31 32]);hold;
        
        if ~isempty(cfg.SpikeWaveforms)
            if ~isempty(cfg.SpikeWaveforms{ipart}{baseline_index}{i_unit})
                
                cfgtemp                     = [];
                cfgtemp.morpho.channame            = cfg.SpikeWaveforms{ipart}{baseline_index}{i_unit}.label{1};
                cfgtemp.morpho.plotstd             = 'yes';
                cfgtemp.morpho.removeoutliers      = 'yes'; %if big noise, impair seeing real data. Still present in avg and std.
                cfgtemp.morpho.mesurehalfwidth     = 'yes';
                cfgtemp.morpho.blmethod            = 'min';
                cfgtemp.morpho.mesurepeaktrough    = 'yes';
                cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.blmethod     = 'min';
                cfgtemp.morpho.toiac               = 'all';
                cfgtemp.morpho.name                = cfg.name{baseline_index};
                [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,cfg.SpikeWaveforms{ipart}{baseline_index}{i_unit});
                
                xlabel('Time (ms)');
                xticklabels(xticks*1000); %convert in ms
                ylabel('Spike waveform (µV)');
                title([]);
                setfig();
                
                stats{ipart}.(cfg.name{baseline_index}).spikewaveform.halfwidth{i_unit}      = halfwidth;
                stats{ipart}.(cfg.name{baseline_index}).spikewaveform.peaktrough{i_unit}     = peaktrough;
                stats{ipart}.(cfg.name{baseline_index}).spikewaveform.troughpeak{i_unit}     = troughpeak;
            else
                axis off;
            end
        end
        
        
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
        
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit},'-spikestats_',cfg.name{baseline_index},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',cfg.SpikeTrials{ipart}{baseline_index}.label{i_unit},'-spikestats_',cfg.name{baseline_index},'.png']),'-r600');
        
        close all
        
    end %i_unit
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
