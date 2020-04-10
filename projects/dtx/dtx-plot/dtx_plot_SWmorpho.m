function halfwidth = dtx_plot_SWmorpho(cfg,data,ipart,imarker,channame, toi, donormalize, saveplot)

%fonction écrite rapidement, à nettoyer
%verif taille de channame

cfgtemp = [];
cfgtemp.channel = channame;
data = ft_selectdata(cfgtemp, data{ipart}{imarker});



% rename prefix in case of "merge" data
if isfield(cfg, 'merge')
    if cfg.merge == true
        if ipart > 1 && ipart == length(cfg.directorylist) %last part = merge (except if only one part, nothing to merge)
            cfg.prefix = [cfg.prefix, 'MERGED-'];
        else
            cfg.prefix = [cfg.prefix, cfg.directorylist{ipart}{:}, '-'];
        end
    end
end



fig = figure;
hold;

%% plot overdraw 
                

for itrial = 1 : size(data.trial,2)
    plot(data.time{itrial},data.trial{itrial}(1,:),'color',[0.6 0.6 0.6]); %first on top
end


set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
axis tight;
xlim(toi);
ylim([-10000 10000]); %TEMPORAIRE, A RETIRER
%ylim([y_min y_max]);
% set(gca,'ycolor',C(imarker,:));
set(gca,'Fontsize',15);

title(sprintf('%s : chan %s (µV)',cfg.LFP.name{imarker},data.label{1}),'Fontsize',10,'Interpreter','none');

%Average chan1

cfgtemp                 = [];
cfgtemp.vartrllength    = 2;
data_EEG_rptavg        = ft_timelockanalysis(cfgtemp,data);

plot(data_EEG_rptavg.time,data_EEG_rptavg.avg(1,:),'r','LineWidth', 2);

t           = data_EEG_rptavg.time;
bl          = mean(data_EEG_rptavg.avg(1,(t>-2&t<1)));
peak        = data_EEG_rptavg.avg(1,(t==0));
halfamp     = (peak-bl)/2;

zci  = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %find time cross zero
indx_temp = zci(data_EEG_rptavg.avg(1,:) - halfamp);
j=0;
for i = 1:length(indx_temp)
    if ~(t(indx_temp(i))<-2 || t(indx_temp(i))>2)
        j = j+1;
        indx(j) = indx_temp(i);
    end
end
indx = indx(1:2);%TEMPORAIRE A RETIRER
% indx = indx([ceil(length(indx)/2-0.5), ceil(length(indx)/2+0.5)]);
plot(data_EEG_rptavg.time(indx),ones(1,length(indx))*halfamp,'-o','Color','b','MarkerFaceColor','b','MarkerEdgeColor','b');

x = sum(data_EEG_rptavg.time(indx))/length(indx);
y = halfamp-halfamp*0.5;
halfwidth = (data_EEG_rptavg.time(indx(2))-data_EEG_rptavg.time(indx(1)))*1000;
text(x,y,sprintf('%.0fms',halfwidth),'HorizontalAlignment','center');

               
               

%% sava data
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprinf('Create forlder %s',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_',channame,'_Morphology.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.LFP.name{imarker},'_',channame,'_Morphology.png']),'-r600');
    close all
end


end

