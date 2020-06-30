function dtx_plot_SlowWaveTopographyTimecouse(cfg,data)
%plot topography of event related to one or two markers, according to 1020
%human EEG layout.
%cfg.labels.emg.
%One figure per marker

% %rename prefix in case of "merge" data
% if isfield(cfg, 'merge')
%     if cfg.merge == true
%         if ipart > 1 && ipart == length(cfg.directorylist) %last part = merge (except if only one part, nothing to merge)
%             cfg.prefix = [cfg.prefix, 'MERGED-'];
%         else
%             cfg.prefix = [cfg.prefix, cfg.directorylist{ipart}{:}, '-'];
%         end
%     end
% end

suffix = ft_getopt(cfg.topoplot, 'suffix', []);
part_list       = ft_getopt(cfg.topoplot, 'part_list'       , 'all');
marker_list     = ft_getopt(cfg.topoplot, 'marker_list'     , 'all');

if strcmp(part_list, 'all')
    part_list = 1:size(data,2);
end

for ipart = part_list 

data = data{ipart};

abscisse_limits = 1; %s

    if strcmp(marker_list, 'all')
        marker_list = 1:size(data{ipart},2);
    end
    
    for imarker = marker_list

% avg
cfgtemp = [];
cfgtemp.channel = {'EEG'};
dat_EEG_avg = ft_timelockanalysis(cfgtemp,data{imarker});

%     % baseline correction
%     cfgtemp = [];
%     cfgtemp.baseline = [-10 -5];
%     dat_EEG_avg_bl = ft_timelockbaseline(cfgtemp,dat_EEG_avg);



%Find zlim according to values at 0

figtemp = figure;

cfgtemp = [];
cfgtemp.layout        = 'EEG1020';
cfgtemp.colorbar      = 'yes';
cfgtemp.zlim          = 'maxabs';
cfgtemp.xlim          = [-0.125 0.125];
cfgtemp.comment       = 'xlim';
cfgtemp.fontsize      = 15;
cfgtemp.renderer      = 'painters';

ft_topoplotER(cfgtemp,dat_EEG_avg);

zaxis_t0 = get(gca,'CLim');
close(figtemp);

%% Topography timecourse

fig1 = figure;

cfgtemp               = [];
cfgtemp.layout        = 'EEG1020';
cfgtemp.colorbar      = 'no';
cfgtemp.zlim          = zaxis_t0;
cfgtemp.xlim          = [-1:0.125:1];
cfgtemp.fontsize      = 15;
cfgtemp.interactive   = 'no';
cfgtemp.renderer      = 'painters';
cfgtemp.comment       = 'xlim';

ft_topoplotER(cfgtemp,dat_EEG_avg);


%% movie
fig2 = figure;
cfgtemp = [];
cfgtemp.layout          = 'EEG1020';
cfgtemp.framespersec    = 1;
cfgtemp.samperframe     = 256;
cfgtemp.colorbar        = 'yes';
cfgtemp.zlim            = zaxis_t0;
cfgtemp.xlim            = [-3 10];
%cfgtemp.samperframe     = 100;


ft_movieplotER(cfgtemp,dat_EEG_avg);





%% print to file
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        warning('%s did not exist for saving images, create now',cfg.imagesavedir);
    end

    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition', [0 0 1 1]);
    set(fig1,'Renderer','Painters');
    print(fig1, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_topography_timecourse',suffix]),'-r600');
    print(fig1, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_topography_timecourse',suffix]),'-r600');
    
    savefig(fig2,fullfile(cfg.imagesavedir,[cfg.prefix,'topography_movie_',cfg.LFP.name{imarker},suffix,'.fig']));
    
    close all
    end
end
end


