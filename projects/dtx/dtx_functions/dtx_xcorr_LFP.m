function dat_xcorr = dtx_xcorr_LFP(cfg, data)

%pour calculer la propagation trial par trial

% cfg.LFP.xcorr.suffix = '_csd';
% cfg.LFP.xcorr.xchan  = 'C4';
% cfg.LFP.xcorr.toi    = [-0.5 0.5];%can be 'all'
% cfg.LFP.xcorr.plotdata = 'yes';

colors = linspecer(size(data.label,1));

cfgtemp         = [];
cfgtemp.channel = cfg.LFP.xcorr.xchan;
cfgtemp.latency = cfg.LFP.xcorr.toi;
data_x_alltrials= ft_selectdata(cfgtemp, data);


if strcmp(cfg.LFP.xcorr.plotdata, 'yes')
    fig = figure('visible','off');hold;
    set(fig, 'units','normalized','position', [0 0 1 0.5]);
end

cfgtemp         = [];
%cfgtemp.trials  = itrial;
data_x          = ft_selectdata(cfgtemp, data_x_alltrials);

i=1;
for iy = 1:size(data.label, 1)
    
    %         if strcmp(cfg.LFP.xcorr.plotdata, 'yes')
    %             subplot(size(data.label, 1),4,[i+1 i+2]);hold;%subplot(size(data.label, 1),size(data.label, 1),i);hold;
    %         end
    
    cfgtemp         = [];
    cfgtemp.channel = iy;
    cfgtemp.latency = cfg.LFP.xcorr.toi;
    data_y          = ft_selectdata(cfgtemp, data);
    
    %[xcorr_trial, x] = xcorr(data_x.trial{1}, data_y.trial{1});
    [xcorr_trial, x] = xcorr(data_x.avg, data_y.avg);
    
    %convert x from sample to seconds
    %Fs = 1/diff(data_y.time{1}(1:2));
    Fs = 1/diff(data_y.time(1:2));
    x = x./Fs;
    
    if strcmp(cfg.LFP.xcorr.xchan, data_y.label{1})
        c = 'b';
    else
        c = 'k';
    end
    
    [~, shift] = findpeaks(xcorr_trial);
    %scatter(x(shift), xcorr_trial(shift), 'rx');
    if isempty(shift)
        lag(iy) = nan;
    else
        shift = shift(x(shift)~=0); %remove peak at zero
        if isempty(shift)
            lag(iy) = 0;
        else
            [~, sel] = max(xcorr_trial(shift));
            shift = shift(sel); %select maximum peak
            if xcorr_trial(shift) > xcorr_trial(x==0)/2
                lag(iy) = x(shift);
            else
                lag(iy) = 0;
            end
        end
    end
    
    if strcmp(cfg.LFP.xcorr.plotdata, 'yes')
        %plot(x,xcorr_trials_norm, c);
        leg{iy} = plot(x,xcorr_trial, 'color', colors(iy,:));
        ax = axis;
        plot([0 0], [ax(3), ax(4)],'r');
        plot([lag(iy) lag(iy)], [ax(3), ax(4)],'b');
        %             xticks([])
        %             yticks([]);
        
        i=i+4;
    end
    
    dat_xcorr.xcorr{iy}         = xcorr_trial;
    dat_xcorr.xlabel{iy}        = data_x.label{1};
    dat_xcorr.ylabel{iy}        = data_y.label{1};
    dat_xcorr.time{iy}          = x;
    dat_xcorr.lag(iy)           = lag(iy);
    
end

legend([leg{:}],data.label');

%print to file, if a figure was created
if strcmp(cfg.LFP.xcorr.plotdata, 'yes')
    fprintf('Print to file\n');
    
    if ~isfolder(fullfile(cfg.imagesavedir,'xcorrs_lfp'))
        mkdir(fullfile(cfg.imagesavedir,'xcorrs_lfp'));
    end
    
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
%     print(fig, '-dpdf', fullfile(cfg.imagesavedir,'xcorrs_lfp',[cfg.prefix,'xcorr_lfp_',cfg.LFP.xcorr.suffix,'.pdf']),'-r600');
%     print(fig, '-dpng', fullfile(cfg.imagesavedir,'xcorrs_lfp',[cfg.prefix,'xcorr_lfp_',cfg.LFP.xcorr.suffix,'.png']),'-r600');
    
%     close all
end

%end %itrial

return
%% Test steps, not used in the function

%plot all non norm xcorr
figure;hold;
for ichan = 1:size(dat_xcorr.avg,2)
    plot(dat_xcorr.time{1},dat_xcorr.avg{ichan});
end
legend(dat_xcorr.ylabel{:});
ax = axis;
plot(cfg.LFP.xcorr.toi, [0 0]);
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