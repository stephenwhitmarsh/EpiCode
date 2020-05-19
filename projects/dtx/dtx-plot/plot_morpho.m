function 	[halfwidth, peaktrough, troughpeak] = plot_morpho(cfg,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [halfwidth, peaktrough, troughpeak] = plot_morpho(cfg,data)
% Plot trials and average of trials. If asked, do measurements on it.
%
% Note :
% - If there are more than 1000 trials, only 1000 random-selected trials
% are plotted. The average remains the average of all the trials.
% - Measurements may not work well if data events are unclear or artefacted
%
% ### INPUT
% data                    = raw Fieldtrip data structure epoched in trials
% cfg.channame            = name of the channel to analyse in data.label
% cfg.removeoutliers      = whether to plot outlier trials (>10*std to the
%                           mean). Average is not modified.
% cfg.mesurehalfwidth     = 'yes' or 'no', whether to compute halfwidth
% cfg.halfwidthmethod     = 'min' or 'bl' : reference for peak-amplitude
%                           measurement if mesurehalfwidth = 'yes'
% cfg.mesurepeaktrough    = 'yes' or 'no', whether to compute peak-trough
%                           and trough-peak
% cfg.toiac               = active period for measurements. Can be 'all'
% cfg.toibl               = baseline period if cfg.halfwidthmethod = 'bl'
% cfg.toiplot             = x limits of plot. Can be 'all'
% cfg.saveplot            = 'yes' or 'no', whether to save and close the
%                           plot or to output it (respectively)
% cfg.name                = name of the analysis (for title of plot and of
%                           saved file)
% cfg.imagesavedir        = if saveplot = 'yes', where to save the plot.
% cfg.prefix              = if saveplot = 'yes', prefix attached to the
%                           name of the saved-image
%
% ### OUTPUT
% halfwidth               = mesured halfwidth value, in seconds, or [].
% peaktrough              = mesured peaktrough value, in seconds, or [].
% troughpeak              = mesured troughpeak value, in seconds, or [].
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg, 'removeoutliers') %TEMPORARY, to avoid error with running script
    cfg.removeoutliers = 'no';
end

%convert some string inputs to usable values
if strcmp(cfg.toiplot,'all')
    cfg.toiplot = [-Inf Inf];
end
if strcmp(cfg.toiac,'all')
    cfg.toiac = [-Inf Inf];
end

%select channel to analyse
cfgtemp         = [];
cfgtemp.channel = cfg.channame;
data            = ft_selectdata(cfgtemp, data);
if isempty(data), error('Channel %s was not found in data', channame); end

%prepare figure
if strcmp(cfg.saveplot, 'yes')
    fig = figure;
end
hold on;

%% plot overdraw and avg

%compute avg 
data_avg        = ft_timelockanalysis([],data);

%select trials to plot, if more than 1000
if size(data.trial,2) > 1000
    trial_list = randperm(size(data.trial,2), 1000);
else
    trial_list = 1 : size(data.trial,2);
end
cfgtemp = [];
cfgtemp.trials = trial_list;
data = ft_selectdata(cfgtemp, data);

%search for outlier trials if asked
isOutlierTrial = false(1, size(data.trial,2));
if strcmp(cfg.removeoutliers, 'yes')
    for itrial = 1:size(data.trial,2)
        if any(data.trial{itrial} > 10 * sqrt(data_avg.var) + abs(data_avg.avg)) || any(data.trial{itrial} < -10 * sqrt(data_avg.var) + abs(data_avg.avg))
            isOutlierTrial(itrial) = true;
        end
    end
end

%plot trials. Reject outlier trials if asked
data.time = data.time(~isOutlierTrial);
data.trial = data.trial(~isOutlierTrial);
for itrial = 1:size(data.time, 2)
    plot(data.time{itrial},data.trial{itrial},'color',[0.6 0.6 0.6]);
end

%plot avg of all trials
plot(data_avg.time,data_avg.avg,'k','LineWidth', 2);

%select active period
cfgtemp         = [];
cfgtemp.latency = cfg.toiac;
data_avg_ac     = ft_selectdata(cfgtemp, data_avg);

%% Mesure and plot half width
halfwidth = [];

if strcmp(cfg.mesurehalfwidth, 'yes')
    
    %measure peak and half hamp
    if strcmp(cfg.halfwidthmethod, 'bl')
        cfgtemp         = [];
        cfgtemp.latency = cfg.toibl;
        data_avg_bl     = ft_selectdata(cfgtemp, data_avg);
        bl          = mean(data_avg_bl.avg);
    elseif strcmp(cfg.halfwidthmethod, 'min')
        bl          = min(data_avg_ac.avg);
    else
        error('%s is not a method for mesuring half width. Set ''bl'' or ''min''',cfg.halfwidthmethod);
    end
    peak        = max(data_avg_ac.avg);
    halfamp     = double(bl+(peak-bl)/2);
    
    %find halfamp indexes in data active period
    zci  = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %find time cross zero
    indx = zci(data_avg_ac.avg - halfamp);
    
    %If found at least 2 indexes, use the 2 first. Otherwise do not compute/plot halfwidth
    if length(indx)>=2
        plot(data_avg_ac.time(indx(1:2)),ones(1,length(indx(1:2)))*(halfamp),'-x','Color','g','MarkerFaceColor','g','MarkerEdgeColor','g');
        
        halfwidth = (data_avg_ac.time(indx(2))-data_avg_ac.time(indx(1)));
        [unit, halfwidth_corr] = setunit(halfwidth); %adapt halfwidth unit
        
        text(data_avg_ac.time(indx(2)),halfamp,sprintf('   %.1f %s',halfwidth_corr,unit),'Color','k','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
    end
end

%% Measure and plot pt and tp
peaktrough = [];
troughpeak = [];

if strcmp(cfg.mesurepeaktrough, 'yes')
    % Find the higher positive peak :
    [~,Xpos] = findpeaks(data_avg_ac.avg, data_avg_ac.time,'NPeaks',1,'SortStr','descend','WidthReference','Halfheight'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
    % Find throughs
    [Yneg,Xneg_temp] = findpeaks(-data_avg_ac.avg,data_avg_ac.time);
    
    if length(Xneg_temp) >= 2
        % Search first through before and after XPos
        [Xneg(1),x_idx(1)] = max(Xneg_temp(Xneg_temp-Xpos < 0));
        [Xneg(2), x_idx(2)] = min(Xneg_temp(Xneg_temp-Xpos > 0));
        Yneg = Yneg(x_idx);
        
        %compute values
        peaktrough = abs(Xpos-Xneg(1));
        troughpeak = abs(Xpos-Xneg(2));
        
        %plot horizontal lines from peak to trough
        plot([Xpos,Xneg(1)], [Yneg(1)*0.3, Yneg(1)*0.3],'-x','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-x','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
        
        %add text for each value
        x = double((Xpos + Xneg(1))/2);
        y = double(Yneg(1)*0.3);
        [unit, peaktrough_corr] = setunit(peaktrough);
        text(x,y,sprintf('%.1f%s   ',peaktrough_corr,unit),'Color','k','HorizontalAlignment','right','VerticalAlignment','middle','FontWeight', 'bold','FontSize',10);
        
        x = double((Xpos + Xneg(2))/2);
        y = double(-Yneg(2)*0.3);
        [unit, troughpeak_corr] = setunit(troughpeak);
        text(x,y,sprintf('%.1f%s',troughpeak_corr,unit),'Color','k','HorizontalAlignment','center','VerticalAlignment','top','FontWeight', 'bold','FontSize',10);
    end
end

axis tight;

set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlabel('Time (s)');
ylabel('µV');
xlim(cfg.toiplot);

if strcmp(cfg.removeoutliers, 'yes')
    title(sprintf('%s chan %s: %d trials (%d outliers removed)', cfg.name, data.label{1},sum(~isOutlierTrial), sum(isOutlierTrial)), 'Fontsize',18, 'Interpreter','none');
else
    title(sprintf('%s chan %s: %d trials', cfg.name, data.label{1},length(data.trial)), 'Fontsize',18, 'Interpreter','none');
end


%% sava data
if strcmp(cfg.saveplot, 'yes')
    
    set(gca,'Fontsize',15);
    
    if ~(exist(cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprintf('Create forlder %s',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',channame,'_morphology_scale',strrep(num2str(cfg.toiplot),'  ','_'),'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',channame,'_morphology_scale',strrep(num2str(cfg.toiplot),'  ','_'),'.png']),'-r600');
    close all
end


end


function [unit, value] = setunit(value)
%convert unit from second to the most relevant unit
unit = 's';
if value < 1
    value = value * 1000;
    unit = 'ms';
end
if value < 1
    value = value * 1000;
    unit = 'us';
end
if value < 1
    value = value * 1000;
    unit = 'ns';
end
end
