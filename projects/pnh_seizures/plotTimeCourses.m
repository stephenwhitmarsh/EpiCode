function plotTimeCourses(cfg,dat_micro,dat_macro)


for imarker = 1 : size(cfg.name,2)
    
    fig = figure;
    
    for ifile = 1 : size(cfg.plot.fname{imarker},2)

        cfgtemp                 = [];
        cfgtemp.dataset         = cfg.plot.fname{imarker}{ifile};
        hdr                     = ft_read_header(cfgtemp.dataset);
        cfgtemp.trl             = cfg.plot.toi{imarker}*hdr.Fs;
        cfgtemp.hpfilter        = cfg.plot.hpfilter{imarker}{ifile};
        cfgtemp.hpfreq          = cfg.plot.hpfreq{imarker}{ifile};
        cfgtemp.lpfilter        = cfg.plot.lpfilter{imarker}{ifile};
        cfgtemp.lpfreq          = cfg.plot.lpfreq{imarker}{ifile};
        
        dat{imarker}{ifile}     = ft_preprocessing(cfgtemp);
        
        if ~isempty(cfg.plot.refname{imarker}{ifile})
            cfgtemp.dataset                 = cfg.plot.refname{imarker}{ifile};
            ref                             = ft_preprocessing(cfgtemp);
            dat{imarker}{ifile}.trial{1}    = dat{imarker}{ifile}.trial{1} - ref.trial{1};
        end
        
        
%         subplot(size(cfg.plot.fname{imarker},2),1,ifile);
        
        subaxis(size(cfg.plot.fname{imarker},2),1,ifile, 'SpacingVert',0.02);

        
        ph = plot(dat{imarker}{ifile}.time{1},dat{imarker}{ifile}.trial{1},'k');
        ph.LineWidth = 0.2;
        axis tight
        ylim(cfg.plot.scale{imarker}{ifile});
        
        ax = gca;
        ax.Clipping = 'off';
        
        box off
        set(gca,'TickDir','out');
        if ifile ~= size(cfg.plot.fname{imarker},2)
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
        end
        
    end
    
    % print to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap   
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'example_timecourses_',cfg.name{imarker},'.pdf']),'-r600');
    close all

end

h = size(cfg.plot.fname{imarker},2);
w = size(cfg.name,2);
w = 3;
p_marker = repmat([1:w],h,1);
p_file   = repmat([1:h]',1,w);

fig = figure;
for imarker = 1 : size(cfg.name,2)
    
    for ifile = 1 : size(cfg.plot.fname{imarker},2)
        
        
%         subplot(h,w,find((imarker == p_marker & ifile == p_file)'));
        subaxis(h,w,find((imarker == p_marker & ifile == p_file)'), 'SpacingVert',0.02);
        ph = plot(dat{imarker}{ifile}.time{1},dat{imarker}{ifile}.trial{1},'k');
        ph.LineWidth = 0.2;
        axis tight
        ylim(cfg.plot.scale{imarker}{ifile});
        box off
        set(gca,'TickDir','out');  
        ax = gca;
        ax.Clipping = 'off'; 
        if ifile ~= size(cfg.plot.fname{imarker},2)
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
        end
    end
end

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'example_timecourses_all.pdf']),'-r600');
close all
