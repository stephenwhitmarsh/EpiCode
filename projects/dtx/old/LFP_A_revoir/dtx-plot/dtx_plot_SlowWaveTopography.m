function dtx_plot_SlowWaveTopography(cfg,data)
%plot topography of event related to one or two markers, according to 1020
%human EEG layout.
%cfg.labels.emg.
%ATTENTION : marker_list must have only 1 or 2 markers
%cfg parameters to save figure : 
% cfg.name
% cfg.prefix
% cfg.imagesavedir

% Other config parameters : 
suffix          = ft_getopt(cfg.topoplot, 'suffix'          , []);
part_list       = ft_getopt(cfg.topoplot, 'part_list'       , 'all');
marker_list     = ft_getopt(cfg.topoplot, 'marker_list'     , 'all');
toi_topoplot 	= ft_getopt(cfg.topoplot, 'toi_topoplot'    , 'all'); % [-2 2]
toi_multiplot 	= ft_getopt(cfg.topoplot, 'toi_multiplot'   , 'all'); % [-2 2]

if strcmp(part_list, 'all')
    part_list = 1:size(data,2);
end

for ipart = part_list
    
    if strcmp(marker_list, 'all')
        marker_list = 1:size(data{ipart},2);
    end
    
    if length(marker_list) > 2
        error('only 2 markers maximum are allowed for this function');
    end
    
    % compute avg
    for imarker = marker_list
        cfgtemp                 = [];
        cfgtemp.channel         = {'EEG'}; %ignore EMG channels
        dat_EEG_avg{imarker}    = ft_timelockanalysis(cfgtemp,data{ipart}{imarker});
    end
    
    fig = figure;
    
    %Test and see if no error without this : 
%     set(gcf, 'Renderer', 'Painters');
%     set(0, 'defaultFigureRenderer', 'painters')
%     set(groot, 'defaultFigureRenderer', 'painters')
    
    %multiplot
    subplot(2,2,[1,3])
    hold;
    
    cfgtemp = [];
    cfgtemp.layout        = 'EEG1020';
    cfgtemp.xlim          = toi_multiplot;
    cfgtemp.axes          = 'yes';
    cfgtemp.box           = 'no';
    cfgtemp.showlabels    = 'yes' ;
    cfgtemp.showoutline   = 'yes' ;
    cfgtemp.showscale     = 'no' ;
    cfgtemp.fontsize      = 15;
    cfgtemp.interactive   = 'no';
    cfgtemp.renderer      = 'painters';
    cfgtemp.linecolor     = 'br';
    cfgtemp.linewidth     = 1;
    cfgtemp.comment       = '\n';
    ft_multiplotER(cfgtemp,dat_EEG_avg{:});
    
    title(sprintf('Slow deflexion topography :'),'Interpreter','none','Fontsize',18);
    
    %topoplot
    for imarker = marker_list
        subplot(2,2,imarker*2)
        
        cfgtemp = [];
        cfgtemp.layout        = 'EEG1020';
        cfgtemp.colorbar      = 'yes';
        cfgtemp.zlim          = 'maxabs';
        cfgtemp.xlim          = toi_topoplot;
        cfgtemp.comment       = 'xlim';
        cfgtemp.fontsize      = 15;
        cfgtemp.renderer      = 'painters';
        ft_topoplotER(cfgtemp,dat_EEG_avg{imarker});
        
        if imarker ==1
            title(sprintf('%s (%d trials) :',cfg.LFP.name{imarker},size(data{ipart}{imarker}.trial,2)),'Interpreter','none','Fontsize',18,'Color','b');
        elseif imarker == 2
            title(sprintf('%s (%d trials) :',cfg.LFP.name{imarker},size(data{ipart}{imarker}.trial,2)),'Interpreter','none','Fontsize',18,'Color','r');
        end
        
    end
    
    %if only one marker, write it in the figure
    if  length(data{ipart}) == 1
        subplot(2,2,4)
        set(gca,'TickLength',[0 0]);
        yticklabels([]);
        xticklabels([]);
        set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
        text(0.2,0.5,sprintf('No controlateral seizure for \n%s', cfg.prefix(1:end-1)),'Interpreter','none','Fontsize',15);
    end
    
    %% print to file
    
    if ~(exist(cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('create dir %s\n',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-topoplot',suffix]),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-topoplot',suffix]),'-r600');
    close all
    
end
end

