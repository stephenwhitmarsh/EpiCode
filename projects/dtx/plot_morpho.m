function 	[halfwidth, peaktrough, troughpeak, amplitude] = plot_morpho(cfg,data)

% Plot trials and average of trials + compute some parameters, see the
% following cfg inputs : 

cfg.morpho                             = ft_getopt(cfg, 'morpho', []);
cfg.morpho.channame                    = ft_getopt(cfg.morpho, 'channame', 1);
cfg.morpho.negpeak                     = ft_getopt(cfg.morpho, 'negpeak'            , 'no');
cfg.morpho.plotstd                     = ft_getopt(cfg.morpho, 'plotstd'            , 'no');
cfg.morpho.plotavg                     = ft_getopt(cfg.morpho, 'plotavg'            , 'yes');
cfg.morpho.plotraw                     = ft_getopt(cfg.morpho, 'plotraw'            , 'yes');
cfg.morpho.raw_facealpha               = ft_getopt(cfg.morpho, 'raw_facealpha'      , 1);
cfg.morpho.measurehalfwidth            = ft_getopt(cfg.morpho, 'measurehalfwidth'	, 'no');
cfg.morpho.measureamplitude            = ft_getopt(cfg.morpho, 'measureamplitude'	, 'no');
cfg.morpho.blmethod                    = ft_getopt(cfg.morpho, 'blmethod' 	        , 'bl');
cfg.morpho.toiplot                     = ft_getopt(cfg.morpho, 'toiplot'            , 'all');
cfg.morpho.toiac                       = ft_getopt(cfg.morpho, 'toiac'              , 'all');
cfg.morpho.saveplot                    = ft_getopt(cfg.morpho, 'saveplot'           , 'no');

if strcmp(cfg.morpho.toiplot,'all')
    cfg.morpho.toiplot = [-Inf Inf];
end
if strcmp(cfg.morpho.toiac,'all')
    cfg.morpho.toiac = [-Inf Inf];
end

%select channel to analyse
cfgtemp         = [];
cfgtemp.channel = cfg.morpho.channame;
data            = ft_selectdata(cfgtemp, data);
if isempty(data), error('Channel %s was not found in data', channame); end

%prepare figure
if strcmp(cfg.morpho.saveplot, 'yes')
    fig = figure;
end
hold on;

%% plot overdraw, avg, std
%flip data if required
if istrue(cfg.morpho.negpeak)
    flip = -1;
else
    flip = 1;
end

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

%plot trials
if strcmp(cfg.morpho.plotraw, 'yes')
    for itrial = 1:size(data.time, 2)
        p = plot(data.time{itrial},data.trial{itrial},'k');
        p.Color(4) = cfg.morpho.raw_facealpha;
    end
end

%plot std
if strcmp(cfg.morpho.plotstd, 'yes')
    std_data = sqrt(data_avg.var);
    x = data_avg.time;
    y = [data_avg.avg - std_data; std_data; std_data]';
    filled_SD = area(x,y);
    filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
    filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
    filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
    filled_SD(1).ShowBaseLine = 'off';
end

%plot avg of all trials
if strcmp(cfg.morpho.plotavg, 'yes')
    if size(data.time, 2) ==1
        plot(data_avg.time, data_avg.avg, 'k', 'LineWidth', 1);
    else
        plot(data_avg.time, data_avg.avg, 'k', 'LineWidth', 2);
    end
end

%select active period
cfgtemp         = [];
cfgtemp.latency = cfg.morpho.toiac;
data_avg_ac     = ft_selectdata(cfgtemp, data_avg);
data_avg_ac.avg = data_avg_ac.avg.*flip;%flip data if required

%% measure and plot half width
halfwidth = NaN;
if strcmp(cfg.morpho.measurehalfwidth, 'yes')
    %measure peak and half hamp
    if strcmp(cfg.morpho.blmethod, 'bl')
        cfgtemp         = [];
        cfgtemp.latency = cfg.morpho.toibl;
        data_avg_bl     = ft_selectdata(cfgtemp, data_avg);
        data_avg_bl.avg = data_avg_bl.avg.*flip;%flip data if required
        bl          = mean(data_avg_bl.avg);
    elseif strcmp(cfg.morpho.blmethod, 'min')
        bl          = min(data_avg_ac.avg).*flip;
    else
        error('%s is not a method for mesuring half width. Set ''bl'' or ''min''',cfg.morpho.blmethod);
    end
    peak        = max(data_avg_ac.avg);
    halfamp     = double(bl+(peak-bl)/2);
    x1 = data_avg_ac.time;
    y1 = data_avg_ac.avg;
    x2 = x1;
    y2 = ones(size(x1)) .* halfamp;
    [x_intersect, y_intersect] = intersections(x1,y1,x2,y2,true);
    if length(x_intersect) < 2
        halfamp = nan;
    elseif all(x_intersect <=0) || all(x_intersect >=0)
        halfamp = nan;
    else
        indx(1) = find(x_intersect <0, 1, 'last');
        indx(2) = find(x_intersect >0, 1, 'first');
        plot(x_intersect(indx),y_intersect(indx).*flip,'-x','Color','g','MarkerFaceColor','g','MarkerEdgeColor','g');
        
        halfwidth = diff(x_intersect(indx));
        [unit, halfwidth_corr] = setunit(halfwidth); %adapt halfwidth unit
        text(x_intersect(indx(2)),halfamp.*flip,sprintf('   %.1f %s',halfwidth_corr,unit),'Color','k','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
    end
end

%% measure amplitude
amplitude = NaN;
if strcmp(cfg.morpho.measureamplitude, 'yes')
    if strcmp(cfg.morpho.blmethod, 'bl')
        cfgtemp         = [];
        cfgtemp.latency = cfg.morpho.toibl;
        data_avg_bl     = ft_selectdata(cfgtemp, data_avg);
        data_avg_bl.avg = data_avg_bl.avg.*flip;%flip data if required
        bl          = mean(data_avg_bl.avg);
    elseif strcmp(cfg.morpho.blmethod, 'min')
        bl          = min(data_avg_ac.avg).*flip;
    else
        error('%s is not a method for mesuring half width. Set ''bl'' or ''min''',cfg.morpho.blmethod);
    end
    [peak, peak_idx]        = max(data_avg_ac.avg);
    amplitude               = (peak-bl).*flip;
    amplitude_loc           = data_avg_ac.time(peak_idx);
end

%% save data
if strcmp(cfg.morpho.saveplot, 'yes')
    
    set(gca,'Fontsize',15);
    
    if ~(exist(cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprintf('Create forlder %s',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',cfg.morpho.channame,'_morpho_scale',strrep(num2str(cfg.morpho.toiplot),'  ','_'),'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name,'_',cfg.morpho.channame,'_morpho_scale',strrep(num2str(cfg.morpho.toiplot),'  ','_'),'.png']),'-r600');
    close all
    
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