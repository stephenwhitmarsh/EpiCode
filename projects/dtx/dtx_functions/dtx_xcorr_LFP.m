function dat_xcorr = test_xcorr_LFP(cfg, data)

%pour calculer la propagation trial par trial

cfg.lfp.xcorr.suffix = '_csd';
cfg.lfp.xcorr.xchan  = 'C4';
cfg.lfp.xcorr.toi    = [-0.5 0.5];
cfg.lfp.xcorr.plotdata = 'yes';

if strcmp(cfg.lfp.xcorr.plotdata, 'yes')
    fig = figure;hold;
    sgtitle(sprintf('x = %s, y = lfp chan indicated',cfg.lfp.xcorr.xchan),'Interpreter','none','Color','r');
    set(fig, 'units','normalized','position', [0 0 1 0.5]);
end

cfgtemp         = [];
cfgtemp.channel = cfg.lfp.xcorr.xchan;
cfgtemp.latency = cfg.lfp.xcorr.toi;
data_x          = ft_selectdata(cfgtemp, data);

i=1;
for iy = 1:size(data.label, 1)
    
    if strcmp(cfg.lfp.xcorr.plotdata, 'yes')
        subplot(size(data.label, 1),4,[i+1 i+2]);hold;%subplot(size(data.label, 1),size(data.label, 1),i);hold;
    end
    
    cfgtemp         = [];
    cfgtemp.channel = iy;
    cfgtemp.latency = cfg.lfp.xcorr.toi;
    data_y           = ft_selectdata(cfgtemp, data);   
    
    for itrial = 1:size(data_x.trial,2)
        xcorr_trials(itrial,:) = xcorr(data_x.trial{itrial}, data_y.trial{itrial});
    end
    
    if strcmp(cfg.lfp.xcorr.xchan, data_y.label{1})
        c = 'b';
    else
        c = 'k';
    end
    
    maxcorr(iy) = max(max(xcorr_trials));
    xcorr_trials_norm = normalize(xcorr_trials, 2, 'range');
    x = linspace(data_x.time{1}(1),data_x.time{1}(end), size(xcorr_trials,2));
    [~, lag_idx] = max(nanmean(xcorr_trials_norm));
    lag(iy) = x(lag_idx);
    
    if strcmp(cfg.lfp.xcorr.plotdata, 'yes')
        plot(x,xcorr_trials_norm, 'Color', [0.6 0.6 0.6]);
        plot(x,nanmean(xcorr_trials_norm), c, 'LineWidth',2);
        ax = axis;
        plot([0 0], [ax(3), ax(4)],'r');
        plot([lag(iy) lag(iy)], [ax(3), ax(4)],'b');
        ylabel(sprintf('%s', data_y.label{1}),'FontSize',5);
        
        title(sprintf('lag = %g, maxcorr = %g', lag(iy),maxcorr(iy)), 'FontSize',5);
        xticks([])
        yticks([]);
        
        i=i+4;
    end
    
    dat_xcorr.avg{iy}      = nanmean(xcorr_trials);
    dat_xcorr.avg_norm{iy} = nanmean(xcorr_trials_norm);
    dat_xcorr.xlabel{iy}   = data_x.label{1};
    dat_xcorr.ylabel{iy}   = data_y.label{1};
    dat_xcorr.time{iy}     = x;
    
end

%print to file, if a figure was created
if strcmp(cfg.lfp.xcorr.plotdata, 'yes')
    fprintf('Print to file and save data\n');
    
    if ~isfolder(fullfile(cfg.imagesavedir,'xcorrs_lfp'))
        mkdir(fullfile(cfg.imagesavedir,'xcorrs_lfp'));
    end
    
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'xcorrs_lfp',[cfg.prefix,'p',num2str(ipart),'-xcorr_lfp_',cfg.name{ilabel},cfg.lfp.xcorr.suffix,'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,'xcorrs_lfp',[cfg.prefix,'p',num2str(ipart),'-xcorr_lfp_',cfg.name{ilabel},cfg.lfp.xcorr.suffix,'.png']),'-r600');
    
    close all
end

return
%% Test steps, not used in the function

%plot all non norm xcorr
figure;hold;
for ichan = 1:size(dat_xcorr.avg,2)
    plot(dat_xcorr.time{1},dat_xcorr.avg{ichan});
end
legend(dat_xcorr.ylabel{:});
ax = axis;
plot(cfg.lfp.xcorr.toi, [0 0]);
plot([0 0], [ax(3) ax(4)]);

%find chans with positive avg correlogram
chan_sel = cellfun(@(v) max(v)>500000,dat_xcorr.avg);
chan_sel = cellfun(@(v) max(v)>0.01,dat_xcorr.avg);

%normalized xcorrs of selected channels
figure;hold;
for ichan = 1:size(dat_xcorr.avg,2)
    if chan_sel(ichan)
        plot(dat_xcorr.time{1},normalize(dat_xcorr.avg_norm{ichan},2,'range'), 'LineWidth',2);
    end
end
ax = axis;
plot([0 0], [ax(3), ax(4)]);
legend(dat_xcorr.ylabel{chan_sel});

%normalized xcorrs of selected channels : normalisation afterwards
figure;hold;
for ichan = 1:size(dat_xcorr.avg,2)
    if chan_sel(ichan)
        plot(dat_xcorr.time{1},normalize(dat_xcorr.avg{ichan},2,'range'), 'LineWidth',2);
    end
end
ax = axis;
plot([0 0], [ax(3), ax(4)]);
legend(dat_xcorr.ylabel{chan_sel});

%plot of real LFP of selected channels
LFP_avg = ft_timelockanalysis([], data);
figure; hold;
plot(LFP_avg.time, LFP_avg.avg((chan_sel)',:)*-1, 'LineWidth',2);
xlim([-0.5 0.5]);
legend(LFP_avg.label{chan_sel});

%test xcorr
a = 1:1000;
b = [100:1000, 1:99];
c = xcorr(a,b);

figure;hold;plot(-999:999,c);
% a faire pour chaque couple d'électrode, identifier s'il y a un pic, et si
% oui calculer son lag.