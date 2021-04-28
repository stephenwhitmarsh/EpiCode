function plotTimeCoursesExamples(cfg)


fig = figure('visible', true);
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
% set(fig, 'PaperOrientation', 'landscape');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0.2 0.2 0.8 0.8]);

ncols   = 4;
icol    = 1;

for markername = string(fields(cfg.plot.toi))'
    
    nrows   = size(cfg.plot.fname.(markername), 2);

    for ifile = 1 : size(cfg.plot.fname.(markername), 2)
        
        cfgtemp                 = [];
        cfgtemp.dataset         = cfg.plot.fname.(markername){ifile};
        hdr                     = ft_read_header(cfgtemp.dataset);
        cfgtemp.trl             = cfg.plot.toi.(markername) * hdr.Fs;
        cfgtemp.hpfilter        = cfg.plot.hpfilter.(markername){ifile};
        cfgtemp.hpfreq          = cfg.plot.hpfreq.(markername){ifile};
        cfgtemp.lpfilter        = cfg.plot.lpfilter.(markername){ifile};
        cfgtemp.lpfreq          = cfg.plot.lpfreq.(markername){ifile};
        dat.(markername){ifile} = ft_preprocessing(cfgtemp);
        
        if ~isempty(cfg.plot.refname.(markername){ifile})
            cfgtemp.dataset                     = cfg.plot.refname.(markername){ifile};
            ref                                 = ft_preprocessing(cfgtemp);
            dat.(markername){ifile}.trial{1}    = dat.(markername){ifile}.trial{1} - ref.trial{1};
        end
                
        subaxis(nrows, ncols, (ifile-1) * ncols + icol, 'SpacingVert', 0.02);

        disp((ifile-1) * ncols + icol)
        ph = plot(dat.(markername){ifile}.time{1}, dat.(markername){ifile}.trial{1}, 'k');
        
        ph.LineWidth = 0.2;
        axis tight
        ylim(cfg.plot.scale.(markername){ifile});
        
        ax = gca;
        ax.Clipping = 'off';
        
        box off
        set(gca,'TickDir','out');
        if ifile ~= size(cfg.plot.fname.(markername), 2)
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'XColor','none')            
        else
            xlabel('Time');
            ylabel('uV');
        end
        set(gca,'FontSize',6);
    end
    icol = icol + 1;
end

% print to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
fname_fig = fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'example_timecourses'));
exportgraphics(fig, [fname_fig, '.jpg'], 'resolution', 600);
exportgraphics(fig, [fname_fig, '.tiff'], 'resolution', 600);
exportgraphics(fig, [fname_fig, '.pdf']);

close all


% 
% h = size(cfg.plot.fname{imarker},2);
% w = size(cfg.name,2);
% w = 3;
% p_marker = repmat([1:w],h,1);
% p_file   = repmat([1:h]',1,w);
% 
% fig = figure;
% 
% for markername = string(fields(cfg.plot.toi))'
%     
%     for ifile = 1 : size(cfg.plot.fname{imarker},2)
%         
%         
% %         subplot(h,w,find((imarker == p_marker & ifile == p_file)'));
%         subaxis(h, w, find((imarker == p_marker & ifile == p_file)'), 'SpacingVert', 0.02);
%         ph = plot(dat{imarker}{ifile}.time{1},dat{imarker}{ifile}.trial{1},'k');
%         ph.LineWidth = 0.2;
%         axis tight
%         ylim(cfg.plot.scale{imarker}{ifile});
%         box off
%         set(gca,'TickDir','out');  
%         ax = gca;
%         ax.Clipping = 'off'; 
%         if ifile ~= size(cfg.plot.fname{imarker},2)
%             set(gca,'xtick',[]);
%             set(gca,'xticklabel',[]);
%         end
%     end
% end
% 
% % print to file
% fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% set(fig,'PaperPosition', [0 0 1 1]);
% print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'example_timecourses_all.pdf']),'-r600');
% close all
