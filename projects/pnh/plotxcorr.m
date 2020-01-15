function plotxcorr(cfg,stat)

fig = figure;
set(fig, 'units','normalized','position', [0 0 1 0.5]);
i = 1;
for ix = 1 : size(stat.xcorr,1)
    for iy = 1 : size(stat.xcorr,2)
        
        if ix > iy
            c = [0 0 0];
        end
        
        if ix < iy
            c = [0 0 0];
        end
        
        if ix == iy
            c = [0 0 0.8];
        end
        
        h = subplot(size(stat.xcorr,1),size(stat.xcorr,2),i);
        bar(stat.time,squeeze(stat.xcorr(ix,iy,:)),1,'facecolor',c,'edgecolor',c);
        axis tight
%         set(h,'yticklabel',{[]});
        t = sprintf('%dx%d',ix,iy);
        title(t);
        pbaspect([1 1 1])
        i = i + 1;
        grid on
    end
end

% print to file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir,['xcorr_',cfg.prefix,'.pdf']),'-r600');

