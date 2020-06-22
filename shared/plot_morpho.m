function 	[halfwidth, peaktrough, troughpeak] = plot_morpho(cfg,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [halfwidth, peaktrough, troughpeak] = plot_morpho(cfg,data)
% Plot trials and average of trials. If asked, do measurements on it.
%
% Note :
% - If there are more than 1000 trials, only 1000 random-selected trials
% are plotted. The average (and std) remain calculated on all the trials.
% - Measurements may not work well if data events are unclear or artefacted
%
% ### Necessary input
% data                    = raw Fieldtrip data structure epoched in trials
% cfg.channame            = name of the channel to analyse (in data.label)
% 
% ### Optionnal cfg fields :
% cfg.plotstd             = 'yes' or 'no', whether to plot std or not.
%                           Default = 'no'.
% cfg.toiplot             = x limits of plot. Can be 'all'. Default = 'all'.
% cfg.saveplot            = 'yes' or 'no', whether to save and close the
%                           plot or to output it (respectively). Default =
%                           'no'.
% cfg.removeoutliers      = whether to plot outlier trials (>10*std to the
%                           mean). Average is not modified. Default = 'no'.
% cfg.mesurehalfwidth     = 'yes' or 'no', whether to compute halfwidth.
%                           Default = 'no'.
% cfg.halfwidthmethod     = 'min' or 'bl' : reference for peak-amplitude
%                           measurement if mesurehalfwidth = 'yes'. 
%                           Default = 'bl' (baseline)
% cfg.mesurepeaktrough    = 'yes' or 'no', whether to compute peak-trough
%                           and trough-peak. Default = 'no'.
% 
% ### Necessary cfg fields if cfg.mesurehalfwidth = 'yes' or cfg.mesurepeaktrough = 'yes'
% cfg.toiac               = active period for measurements. Can be 'all'
% cfg.toibl               = baseline period if cfg.halfwidthmethod = 'bl'
% 
% ### Necessary cfg fields if cfg.saveplot = 'yes'
% cfg.name                = name of the analysis (title of plot and of
%                           output image file)
% cfg.imagesavedir        = where to save the image.
% cfg.prefix              = prefix attached to the name of the saved-image
%
% ### OUTPUT
% halfwidth               = mesured halfwidth value, in seconds, or [].
% peaktrough              = mesured peaktrough value, in seconds, or [].
% troughpeak              = mesured troughpeak value, in seconds, or [].
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Get default cfg parameters
cfg.plotstd                     = ft_getopt(cfg, 'plotstd'          , 'no');
cfg.removeoutliers              = ft_getopt(cfg, 'removeoutliers'  	, 'no');
cfg.mesurehalfwidth             = ft_getopt(cfg, 'mesurehalfwidth'	, 'no');
cfg.halfwidthmethod             = ft_getopt(cfg, 'halfwidthmethod' 	, 'bl');
cfg.mesurepeaktrough            = ft_getopt(cfg, 'mesurepeaktrough'	, 'no');
cfg.toiplot                     = ft_getopt(cfg, 'toiplot'          , 'all');
cfg.saveplot                    = ft_getopt(cfg, 'saveplot'         , 'no');

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

%% plot overdraw, avg, std

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


%plot std if required
if strcmp(cfg.plotstd, 'yes')
    %plot(data_avg.time,data_avg.avg + sqrt(data_avg.var),'k','LineWidth', 2);
    %plot(data_avg.time,data_avg.avg - sqrt(data_avg.var),'k','LineWidth', 2);
    std = sqrt(data_avg.var);
    x = data_avg.time;
    y = [data_avg.avg - std; std; std]';
    filled_SD = area(x,y);
    filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
    filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
    filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
    filled_SD(1).ShowBaseLine = 'off';
end

%plot avg of all trials
plot(data_avg.time, data_avg.avg, 'k', 'LineWidth', 2);


%select active period
cfgtemp         = [];
cfgtemp.latency = cfg.toiac;
data_avg_ac     = ft_selectdata(cfgtemp, data_avg);

%% Mesure and plot half width
halfwidth = NaN;

if strcmp(cfg.mesurehalfwidth, 'yes')
    try %TEMPORARY - CHECK DTX4 unit just after 53
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
            
            %increase time precision *1000 by linear interpolation between 2 samples
            sel_1           = data_avg_ac.avg(indx(1)) : (data_avg_ac.avg(indx(1)+1)-data_avg_ac.avg(indx(1)))/1000 : data_avg_ac.avg(indx(1)+1);
            sel_2           = data_avg_ac.avg(indx(2)) : (data_avg_ac.avg(indx(2)+1)-data_avg_ac.avg(indx(2)))/1000 : data_avg_ac.avg(indx(2)+1);
            indx_temp(1)    = find(sel_1 > halfamp, 1, 'first');
            indx_temp(2)    = find(sel_2 < halfamp, 1, 'first');
            samp_interv     = data_avg_ac.time(2) - data_avg_ac.time(1);
            x_precise(1)    = data_avg_ac.time(indx(1)) + indx_temp(1)*samp_interv/1000;
            x_precise(2)    = data_avg_ac.time(indx(2)) + indx_temp(2)*samp_interv/1000;
            
            plot(x_precise,ones(1,length(indx(1:2)))*(halfamp),'-x','Color','g','MarkerFaceColor','g','MarkerEdgeColor','g');
            
            halfwidth = x_precise(2)-x_precise(1);
            [unit, halfwidth_corr] = setunit(halfwidth); %adapt halfwidth unit
            
            text(x_precise(2),halfamp,sprintf('   %.1f %s',halfwidth_corr,unit),'Color','k','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
        end
    end %REMOVEME
end

%% Measure and plot pt and tp
peaktrough = NaN;
troughpeak = NaN;

if strcmp(cfg.mesurepeaktrough, 'yes')
    try %TEMPORARY - CHECK DTX7 unit just after 26
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
        x = Xneg(1);
        y = double(Yneg(1)*0.3);
        [unit, peaktrough_corr] = setunit(peaktrough);
        text(x,y,sprintf('%.1f%s   ',peaktrough_corr,unit),'Color','k','HorizontalAlignment','right','VerticalAlignment','middle','FontWeight', 'bold','FontSize',10);
        
        x = Xneg(2);
        y = double(-Yneg(2)*0.3);
        [unit, troughpeak_corr] = setunit(troughpeak);
        text(x,y,sprintf('   %.1f%s',troughpeak_corr,unit),'Color','k','HorizontalAlignment','left','VerticalAlignment','middle','FontWeight', 'bold','FontSize',10);
    end
     end %REMOVEME
end

axis tight;

set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
xlabel('Time (s)');
ylabel('µV');
xlim(cfg.toiplot);

if strcmp(cfg.removeoutliers, 'yes')
    title(sprintf('Chan %s: %d trials (%d outliers removed)', data.label{1},sum(~isOutlierTrial), sum(isOutlierTrial)), 'Fontsize',18, 'Interpreter','none');
else
    title(sprintf('Chan %s: %d trials', data.label{1},length(data.trial)), 'Fontsize',18, 'Interpreter','none');
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
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',cfg.channame,'_morpho_scale',strrep(num2str(cfg.toiplot),'  ','_'),'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',cfg.channame,'_morpho_scale',strrep(num2str(cfg.toiplot),'  ','_'),'.png']),'-r600');
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
