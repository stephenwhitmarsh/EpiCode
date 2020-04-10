function dtx_plot_comparison_severaleeg(cfg,datalist,ipart,channames, toi, donormalize, saveplot)
% datalist : cell array with the data to compare {data1, data2, data3}
%channames : cell array with the channel names to use for each data (one channel per data) {name1, name2, name3}
%toi : time of interest [low, high]

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


if length(datalist) ~= length(channames)
    error('List of data must have the same size as list of channel names');
end
fig = figure;
C = linspecer(length(datalist));
n_subplots = round((length(datalist))/2)+1;


for idata = 1:length(datalist)
    if ~isempty(datalist{idata})
        if ~isempty(datalist{idata}.trial)
            if ~strcmp(channames{idata},'no') && ~isempty(channames{idata})
                subplot(n_subplots,2,idata);
                hold;
                
                %     if idata == 1
                %         title(sprintf('%s - EEG comparison : %d trials \nChannel %s',cfg.prefix(1:end-1),length(datalist{idata}.trial),channames{idata}),'Fontsize',18,'Interpreter','none');
                %     end
                
                %select data
                cfgtemp = [];
                cfgtemp.channel = channames{idata};
                data_EEG{idata} = ft_selectdata(cfgtemp,datalist{idata});
                
                %find good Y lim to avoid flattening by too much noise
                for itrial = 1 : length(data_EEG{idata}.trial)
                    t_0 = -(cfg.epoch.toi{1}(1)-cfg.epoch.pad{1}(1))*data_EEG{idata}.fsample; % offset for which t = 0;
                    y_max_temp(itrial) = max(data_EEG{idata}.trial{itrial}(1,round(-0.5*data_EEG{idata}.fsample)+t_0: round(0.5*data_EEG{idata}.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
                    y_min_temp(itrial) = min(data_EEG{idata}.trial{itrial}(1,round(-0.5*data_EEG{idata}.fsample)+t_0: round(0.5*data_EEG{idata}.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
                end
                y_max = max(y_max_temp);
                y_min = min(y_min_temp);
                
                
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                for itrial = 1 : size(data_EEG{idata}.trial,2)
                    plot(data_EEG{idata}.time{itrial},data_EEG{idata}.trial{itrial}(1,:),'color',[0.6 0.6 0.6]); %first on top
                end
                
                
                set(gca,'FontWeight','bold' );
                set(gca,'TickDir','out');
                axis tight;
                xlim(toi);
                %ylim([y_min y_max]);
                set(gca,'ycolor',C(idata,:));
                set(gca,'Fontsize',15);
                xticks([]);
                if strcmp(cfg.LFP.name{idata},data_EEG{idata}.label{1})
                    title(sprintf('%s (µV)',cfg.LFP.name{idata}),'Fontsize',10,'Interpreter','none');
                else
                    title(sprintf('%s : chan %s (µV)',cfg.LFP.name{idata},data_EEG{idata}.label{1}),'Fontsize',10,'Interpreter','none');
                end
                
                %Average chan1
                
                cfgtemp                 = [];
                cfgtemp.vartrllength    = 2;
                data_EEG_rptavg{idata}        = ft_timelockanalysis(cfgtemp,data_EEG{idata});
                
                plot(data_EEG_rptavg{idata}.time,data_EEG_rptavg{idata}.avg(1,:),'color',C(idata,:),'LineWidth', 2);
            end
        end
        
    end
end

%% Overdraw chans avg

subplot(n_subplots, 2, [n_subplots*2-1, n_subplots*2])
hold;

for idata = 1:length(datalist)
    if ~isempty(datalist{idata})
        if ~isempty(datalist{idata}.trial)
            if ~strcmp(channames{idata},'no') && ~isempty(channames{idata})
                
                if donormalize
                    data_EEG_rptavg{idata}.avg(1,:) = data_EEG_rptavg{idata}.avg(1,:) / max(data_EEG_rptavg{idata}.avg(1,:));
                end
                
                plot(data_EEG_rptavg{idata}.time,data_EEG_rptavg{idata}.avg(1,:),'color',C(idata,:),'LineWidth', 2);
                
            end
        end
    end
end

%ylabel(sprintf('EEG average (µV) \n %s %s %s %s %s %s %s %s %s %s',channames{:}), 'Fontsize',15);

axis tight
xlim(toi);

set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
set(gca,'Fontsize',15);

xlabel('Time (s)','Fontsize',18);


%% sava data
if saveplot
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprinf('Create forlder %s',cfg.imagesavedir);
    end
    
    eegname = [];
    for idata = 1:length(datalist)
        if ~strcmp(channames{idata},'no') && ~isempty(channames{idata})
            if idata == 1
                eegname = [eegname, channames{idata}];
            else
                eegname = [eegname, '_', channames{idata}];
            end
        end
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'ComparisonSeveralEEG',eegname,'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'ComparisonSeveralEEG','_comparisonseveraleeg_',eegname,'.png']),'-r600');
    close all
end


end

