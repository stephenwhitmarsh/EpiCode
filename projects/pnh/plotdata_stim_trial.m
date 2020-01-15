function plotdata_stim_trial(cfg,dat_micro, dat_macro)

% select data for plotting

if ~strcmp(cfg.frequency,'all')
    selfreq = cfg.markers.Frequency == cfg.frequency;
else
    selfreq = ones(height(cfg.markers),1);
end

if ~strcmp(cfg.current,'all')
    selamp = cfg.markers.Current == cfg.current;
else
    selamp = ones(height(cfg.markers),1);
end

if ~strcmp(cfg.contact1,'all')
    selcontact1 = [];
    for i = 1 : height(cfg.markers)
        if strcmp(cfg.markers.contact1(i,:),cfg.contact1)
            selcontact1(i,1) = 1;
        else
            selcontact1(i,1) = 0;
        end
    end
else
    selcontact1 = ones(height(cfg.markers),1);
end

if ~strcmp(cfg.contact2,'all')
    selcontact2 = [];
    for i = 1 : height(cfg.markers)
        if strcmp(cfg.markers.contact2(i,:),cfg.contact2)
            selcontact2(i,1) = 1;
        else
            selcontact2(i,1) = 0;
        end
    end
else
    selcontact2 = ones(height(cfg.markers),1);
end


dat_micro = ft_appenddata([],dat_micro{selfreq&selamp&selcontact1&selcontact2});
% 
% cfgtemp = [];
% cfgtemp.trials = find(selfreq&selamp&selcontact1&selcontact2);
% dat_micro = ft_selectdata(cfgtemp,markerdat_micro);
% dat_macro = ft_appenddata([],markerdat_macro{selfreq&selamp&selcontact1&selcontact2});
cfg.markers = cfg.markers(selfreq&selamp&selcontact1&selcontact2,:);

% clear markerdat_micro markerdat_macro

% select channels and filter
cfgtemp                 = [];
cfgtemp.hpfilter        = cfg.hpfilter;
cfgtemp.hpfreq          = cfg.hpfreq;
cfgtemp.channel         = cfg.channel_micro;
dat_micro               = ft_preprocessing(cfgtemp,dat_micro);

% cfgtemp.channel         = cfg.channel_macro;
% dat_macro               = ft_preprocessing(cfgtemp,dat_macro);

cfgtemp = [];
cfgtemp.latency         = cfg.latency;
dat_micro_temp          = ft_selectdata(cfgtemp,dat_micro);
% dat_macro_temp          = ft_selectdata(cfgtemp,dat_macro);

% 
% cfgtemp = [];
% cfgtemp.viewmode = 'vertical';
% ft_databrowser(cfgtemp,dat_micro)

set(0, 'DefaultTextInterpreter', 'none')
for istim = 1 : size(dat_micro_temp.trial,2)
    
    fig = figure;
    %     subplot(2,1,1);
    %     hold;
    %     for ichan = 1 : size(dat_macro_temp.label,1)
    %         plot(dat_macro_temp.time{istim},dat_macro_temp.trial{istim}(ichan,:));
    %     end
    %     title(sprintf('MACRO %s %dHz %0.1fmA',cfg.markers.Label{istim},cfg.markers.Frequency(istim,:),cfg.markers.Current(istim,:)));
    
%     subplot(2,1,2);
    hold;
    for ichan = 1 : size(dat_micro_temp.label,1)
        plot(dat_micro_temp.time{istim},dat_micro_temp.trial{istim}(ichan,:));
    end
    titletext = sprintf('MICRO_%s_%dHz_%0.1fmA',cfg.markers.Label{istim},cfg.markers.Frequency(istim,:),cfg.markers.Current(istim,:));
    title(titletext);
    
    % write to PDF
    fig.Renderer='Painters';
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[titletext,'.pdf']),'-r600');
    %
end
