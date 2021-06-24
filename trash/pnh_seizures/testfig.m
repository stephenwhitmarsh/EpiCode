function testfig(stat_x)
set(gcf,'Renderer','zbuffer')
f = figure;
set(f, 'units','normalized','position', [0 0 1 0.5]);
i = 1;
for ix = 1 : size(stat_x.xcorr,1)
    for iy = 1 : size(stat_x.xcorr,2)
        
        if ix > iy
            c = [0 0.8 0];
        end
        
        if ix < iy
            c = [0 0.8 0];
        end
        
        if ix == iy
            c = [0 0 0.8];
        end
        
        h = subplot(size(stat_x.xcorr,1),size(stat_x.xcorr,2),i);
        
        bar(stat_x.time,squeeze(stat_x.xcorr(iy,ix,:)),'facecolor',c);
        axis tight
        %         axis square
%         set(h,'xticklabel',{[]});
        set(h,'yticklabel',{[]});
        t = sprintf('%dx%d',ix,iy);
        title(t);
        %         p = [(1/size(stat_x.xcorr,1))*(ix), (1/size(stat_x.xcorr,2))*(iy-1),(1/size(stat_x.xcorr,1))*0.5,(1/size(stat_x.xcorr,2))*0.5]
        
        
        %         set(h,'pos', p);
        pbaspect([1 1 1])
        i = i + 1;
        
    end
end

