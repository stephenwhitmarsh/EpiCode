function dtx_plot_SlowWaveTopography(cfg,data,separate_R_L,saveplot)
%plot topography of event related to one or two markers, according to 1020
%human EEG layout.
%cfg.labels.emg.

%ATTENTION : if there are more than 2 markers which define trials, this
%script is not yet adapted.

%if separate_R_L = true, subplot right and left and indicate absence of one
%side of seizure if so

abscisse_limits = 2; %s



for imarker = 1:length(data)
    % avg
    cfgtemp = [];
    cfgtemp.channel = {'EEG'};
    dat_temp = ft_timelockanalysis(cfgtemp,data{imarker});
    dat_EEG_avg{imarker} = dat_temp;
    
%     % baseline correction
%     cfgtemp = [];
%     cfgtemp.baseline = [-10 -5];
%     dat_temp = ft_timelockbaseline(cfgtemp,dat_temp);
%     dat_EEG_avg_bl{imarker} = dat_temp;
end

fig = figure;

set(gcf, 'Renderer', 'Painters');
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

%% topography Right and Left

%EEG avg
subplot(2,2,[1,3])
hold;

cfgtemp = [];
cfgtemp.layout        = 'EEG1020';
cfgtemp.xlim          = [-abscisse_limits abscisse_limits];
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

if  separate_R_L == true && length(data)>1
    ft_multiplotER(cfgtemp,dat_EEG_avg{1},dat_EEG_avg{2}); %script adapted only for 2 markers maximum
    
else
    ft_multiplotER(cfgtemp,dat_EEG_avg{1});
    
end

title(sprintf('Slow deflexion topography :'),'Interpreter','none','Fontsize',18);


% Topography power

for imarker = 1:length(data)
    
    
    subplot(2,2,imarker*2)
    
    cfgtemp = [];
    cfgtemp.layout        = 'EEG1020';
    cfgtemp.colorbar      = 'yes';
    cfgtemp.zlim          = 'maxabs';%[-100 180];%'maxabs';
    cfgtemp.xlim          = [-0.125 0.125];
    cfgtemp.comment       = 'xlim';
    cfgtemp.fontsize      = 15;
    cfgtemp.renderer      = 'painters';
    
    ft_topoplotER(cfgtemp,dat_EEG_avg{imarker});
    
    if imarker ==1
        title(sprintf('%d %s :',size(data{imarker}.trial,2),cfg.LFP.name{imarker}),'Interpreter','none','Fontsize',18,'Color','b');
    elseif imarker == 2
        title(sprintf('%d %s :',size(data{imarker}.trial,2),cfg.LFP.name{imarker}),'Interpreter','none','Fontsize',18,'Color','r');
    end


end

if  separate_R_L == true && length(data) == 1
    subplot(2,2,4)
    set(gca,'TickLength',[0 0]);
    yticklabels([]);
    xticklabels([]);
    set(gca,'YColor', [1 1 1],'XColor', [1 1 1]);
    text(0.2,0.5,sprintf('No controlateral seizure for \n%s', cfg.prefix(1:end-1)),'Interpreter','none','Fontsize',15);
end




%% print to file
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'topography']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'topography']),'-r600');
    close all
end

end

