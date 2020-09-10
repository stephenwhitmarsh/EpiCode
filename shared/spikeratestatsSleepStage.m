
function [stats] = spikeratestatsSleepStage(cfg, SpikeRaw, SpikeTrials, Hypnogram, force)

% SPIKERATESTATSEVENTS calculates and plots spike statistics
%
% use as
%   spikeratestatsEvents(cfg,SpikeRaw,SpikeTrials,force,varargin)

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikestatsSleepStages.mat']);
if exist(fname,'file') && force == false
    load(fname,'stats');
else

    for ipart = 1:length(SpikeRaw)

        % fix this for consistency
        temp                = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,['p',num2str(ipart)],'-multifile-',cfg.circus.channel{1},'.ncs']));
        hdr_fname           = fullfile(temp(1).folder,temp(1).name);
        hdr                 = ft_read_header(hdr_fname); % take the first file to extract the header of the data

        % ISI on all data
        fig = figure;
        for itemp = 1 : length(SpikeRaw{ipart}.label)
            stats{ipart}.isi_all_data{itemp} = diff(double(SpikeRaw{ipart}.samples{itemp})) / hdr.Fs;
            subplot(round(length(SpikeRaw{ipart}.label)/2+0.25),2,itemp);
            histogram(stats{ipart}.isi_all_data{itemp}*1000,'BinWidth',0.5,'BinLimits',[0, 50],'LineStyle','None','FaceColor',[0 0 0]);
            [y,indx] = max(stats{ipart}.isi_all_data{itemp});
            title(sprintf('Unit: %d, Max ISI: %.1fms',itemp,stats{ipart}.isi_all_data{itemp}(indx)));
            yticks([]);
            axis tight
            set(gca,'fontsize',6);
        end

        % print ISI to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'ISI_all_data_sleepstage.pdf']),'-r600');
        close all

        % ISI per stage
        cfgtemp                                 = [];
        cfgtemp.trials                          = 'all';
        cfgtemp.outputunit                      = 'proportion';
        cfgtemp.bins                            = [0 : 0.0005 : 0.100]; %cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
        cfgtemp.param                           = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials                      = 'yes';
        stats{ipart}.isi_all_trials             = ft_spike_isi(cfgtemp,SpikeTrials{ipart});

        % ISI descriptives per sleep stage
        for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'

            % ISI
            cfgtemp                                 = [];
            cfgtemp.trials                          = find(SpikeTrials{ipart}.trialinfo.stage == istage);
            cfgtemp.outputunit                      = 'proportion';
            cfgtemp.bins                            = [0 : 0.0005 : 0.100]; %cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
            cfgtemp.param                           = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
            stats{ipart}.isi_sleepstage{istage+2}   = ft_spike_isi(cfgtemp,SpikeTrials{ipart});

            %     RPV = (length(find(ISI < 2)) / length(ISI)) * 100

            % plot ISI for each cluster
            fig = figure; hold;
            for itemp = 1 : length(SpikeRaw{ipart}.label)
                subplot(round(length(SpikeRaw{ipart}.label)/2+0.5),2,itemp);
                bar(stats{ipart}.isi_sleepstage{istage+2}.time*1000,stats{ipart}.isi_sleepstage{istage+2}.avg(itemp,:),1);
                [y,indx] = max(stats{ipart}.isi_sleepstage{istage+2}.avg(itemp,:));
                title(sprintf('Unit: %d, Max ISI: %.1fms',itemp,stats{ipart}.isi_sleepstage{istage+2}.time(indx)*1000));
                xlabel('ms');
                axis tight
                set(gca,'fontsize',6);
            end

            % print to file
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-ISI_FT_sleepstage',num2str(istage),'.pdf']),'-r600');
            close all

            % create stats for each cluster x stage
            for itemp = 1 : length(SpikeTrials{ipart}.label)

                trials = find(SpikeTrials{ipart}.trialinfo.stage == istage);
                i = 1;
                clear trialavg_isi
                for itrial = trials'
                    indx    = SpikeTrials{ipart}.trial{itemp} == itrial;
                    isi     = stats{ipart}.isi_all_trials.isi{itemp}(indx);
                    trialavg_isi(i) = nanmean(isi);
                    short(i) = sum(isi < 0.005);
                    long(i)  = sum(isi < 0.100);
                    i = i + 1;
                end
                stats{ipart}.mean_freq{itemp}(istage+2)     = 1 / nanmean(trialavg_isi);
                stats{ipart}.stdev_freq{itemp}(istage+2)    = 1 / nanstd(trialavg_isi);
                stats{ipart}.median_freq{itemp}(istage+2)   = 1 / nanmedian(trialavg_isi);
                stats{ipart}.mean_isi{itemp}(istage+2)      = nanmean(trialavg_isi);
                stats{ipart}.stdev_isi{itemp}(istage+2)     = nanstd(trialavg_isi);
                stats{ipart}.var_isi{itemp}(istage+2)       = nanstd(trialavg_isi)^2;
                stats{ipart}.fano_isi{itemp}(istage+2)      = nanstd(trialavg_isi)^2  / nanmean(trialavg_isi);
                stats{ipart}.CV_isi{itemp}(istage+2)        = nanstd(trialavg_isi) / nanmean(trialavg_isi);
                stats{ipart}.burstindex{itemp}(istage+2)    = sum(short) / sum(long);

            end
        end


        % ISI descriptives combined over sleep stages
        fig = figure; hold;
        C = linspecer(5); % colormap for nr. of sleep stages

        % plot ISI for each cluster
        for itemp = 1 : length(SpikeRaw{ipart}.label)
            subplot(round(length(SpikeRaw{ipart}.label)/2+0.5),2,itemp); hold;

            for istage = 0 : 4 %unique(SpikeTrials{ipart}.trialinfo.stage)' % to remove -1, i.e. all non-scored

                bar(stats{ipart}.isi_sleepstage{istage+2}.time*1000,stats{ipart}.isi_sleepstage{istage+2}.avg(itemp,:),1,'FaceColor',C(istage+1,:),'FaceAlpha',0.5);
                [y,indx] = max(stats{ipart}.isi_sleepstage{istage+2}.avg(itemp,:));
                title(sprintf('Unit: %d, Max ISI: %.1fms',itemp,stats{ipart}.isi_sleepstage{istage+2}.time(indx)*1000));
                xlabel('ms');
                axis tight
                set(gca,'fontsize',6);
            end
            legend({'W','Stage 1','Stage 2','Stage 3','REM'},'location','eastoutside');
        end

        % print to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-ISI_FT_sleepstages_combined.pdf']),'-r600');
        %         print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-ISI_FT_sleepstages_combined.png']),'-r600');

        close all




        % plot hypnogram with spikerate for each cluster
        Hypnogram_part = Hypnogram(Hypnogram.part == ipart,:);
        for itemp = 1 : length(SpikeTrials{ipart}.label)

            fig = figure;
            subplot(3,1,1); hold;

            X = [];
            Y = [];
            for im = 1 : height(Hypnogram_part)
                if ~isempty(X)
                    % if there's a gap, 'fill' with 0
                    if Hypnogram_part.starttime(im) ~= X(end)
                        X = [X, X(end) Hypnogram_part.starttime(im)];
                        Y = [Y, 0,  0];
                    end
                end
                X = [X, Hypnogram_part.starttime(im), Hypnogram_part.endtime(im)];

                % height in hypnogram is based on order of config.hyp.contains
                switch cell2mat(Hypnogram_part.stagelabel(im))
                    case 'AWAKE'
                        y = 5;
                    case 'REM'
                        y = 4;
                    case 'STAGE 1'
                        y = 3;
                    case 'STAGE 2'
                        y = 2;
                    case 'STAGE 3'
                        y = 1;
                end
                Y = [Y, y, y];
            end

            for i = 1 : length(X)-1
                if Y(i) ~= 0 && Y(i+1) ~= 0
                    if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
                        plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
                    else
                        plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
                    end
                end
            end
            set(gca,'Layer','top');
            set(gca,'Ytick', 1 : 5,'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'},'TickDir','out');
            title('Hypnogram');
            axis tight;
            xl = xlim;
            ylim([0 6]);

            % create stats for each trail (no sepation on stages)
            for itrial = 1 : height(SpikeTrials{ipart}.trialinfo)
                indx    = SpikeTrials{ipart}.trial{itemp} == itrial;
                isi     = stats{ipart}.isi_all_trials.isi{itemp}(indx);
                stats{ipart}.trialavg_isi{itemp}(itrial) = nanmean(isi);
                stats{ipart}.trialavg_freq{itemp}(itrial) = 1 / nanmean(isi);
            end

            % plot spikerates
            subplot(3,1,2);
            x = (SpikeTrials{ipart}.trialinfo.starttime + (SpikeTrials{ipart}.trialinfo.endtime - SpikeTrials{ipart}.trialinfo.starttime)/2);
            plot(x,log(stats{ipart}.trialavg_freq{itemp}));
            xlim(xl); ylabel('Log( Firingrate(Hz) )');
            title('Firingrate');

            % plot stats
            subplot(3,1,3);
            stageindx = [zeros(1,length(stats{ipart}.isi_sleepstage{1}.isi{1})) ...
                ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{2})), ...
                ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{3}))*2, ...
                ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{4}))*3, ...
                ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{5}))*4, ...
                ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{6}))*5];
            stageISI = [stats{ipart}.isi_sleepstage{1}.isi{1}, ...
                stats{ipart}.isi_sleepstage{1}.isi{2}, ...
                stats{ipart}.isi_sleepstage{1}.isi{3}, ...
                stats{ipart}.isi_sleepstage{1}.isi{4}, ...
                stats{ipart}.isi_sleepstage{1}.isi{5}, ...
                stats{ipart}.isi_sleepstage{1}.isi{6}];

            subplot(3,4,9);
            hold on
            errorbar(1:6,stats{ipart}.mean_freq{itemp},zeros(size(stats{ipart}.stdev_freq{itemp})),stats{ipart}.stdev_freq{itemp},'k','LineStyle','none');
            bar(1:6,stats{ipart}.mean_freq{itemp});
            set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
            axis tight
            box off
            xlim([1.5,6.5]); % remove x, aka 'rest of data'
            title('Freq+SD');

            subplot(3,4,10);
            bar(1:6,stats{ipart}.fano_isi{itemp});
            set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
            axis tight
            box off
            xlim([1.5,6.5]); % remove x, aka 'rest of data'
            title('FANO');

            subplot(3,4,11);
            bar(1:6,stats{ipart}.CV_isi{itemp});
            set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
            axis tight
            box off
            xlim([1.5,6.5]); % remove x, aka 'rest of data'
            title('CV1');

            subplot(3,4,12);
            bar(1:6,stats{ipart}.burstindex{itemp});
            set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
            axis tight
            box off
            xlim([1.5,6.5]); % remove x, aka 'rest of data'
            title('BI');

            % print to file
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-cluster',num2str(itemp),'_hypnospikestats.pdf']),'-r600');
            close all

        end


        %% plot hypnogram with spikerate for all clusters at once, per night (part)
        Hypnogram_part = Hypnogram(Hypnogram.part == ipart,:);

        fig = figure;
        subplot(3,1,1); hold;

        X = [];
        Y = [];
        for im = 1 : height(Hypnogram_part)
            if ~isempty(X)
                % if there's a gap, 'fill' with 0
                if Hypnogram_part.starttime(im) ~= X(end)
                    X = [X, X(end) Hypnogram_part.starttime(im)];
                    Y = [Y, 0,  0];
                end
            end
            X = [X, Hypnogram_part.starttime(im), Hypnogram_part.endtime(im)];

            % height in hypnogram is based on order of config.hyp.contains
            switch cell2mat(Hypnogram_part.stagelabel(im))
                case 'AWAKE'
                    y = 5;
                case 'REM'
                    y = 4;
                case 'STAGE 1'
                    y = 3;
                case 'STAGE 2'
                    y = 2;
                case 'STAGE 3'
                    y = 1;
            end
            Y = [Y, y, y];
        end

        for i = 1 : length(X)-1
            if Y(i) ~= 0 && Y(i+1) ~= 0
                if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
                    plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
                else
                    plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
                end
            end
        end
        set(gca,'Layer','top');
        set(gca,'Ytick', 1 : 5,'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'},'TickDir','out');
        title('Hypnogram');
        axis tight;
        xl = xlim;
        ylim([0 6]);
        legend('stages','location','eastoutside');


        % plot spikerates
        subplot(3,1,2); hold;
        C = linspecer(length(SpikeTrials{ipart}.label));

        x = (SpikeTrials{ipart}.trialinfo.starttime + (SpikeTrials{ipart}.trialinfo.endtime - SpikeTrials{ipart}.trialinfo.starttime)/2);
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(x,log(stats{ipart}.trialavg_freq{itemp}),'color',C(itemp,:));
        end
        xlim(xl); ylabel('Log( Firingrate(Hz) )');
        title('Firingrate');
        legend('location','eastoutside');

        % plot stats
        stageindx = [zeros(1,length(stats{ipart}.isi_sleepstage{1}.isi{1})) ...
            ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{2})), ...
            ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{3}))*2, ...
            ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{4}))*3, ...
            ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{5}))*4, ...
            ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{6}))*5];
        stageISI = [stats{ipart}.isi_sleepstage{1}.isi{1}, ...
            stats{ipart}.isi_sleepstage{1}.isi{2}, ...
            stats{ipart}.isi_sleepstage{1}.isi{3}, ...
            stats{ipart}.isi_sleepstage{1}.isi{4}, ...
            stats{ipart}.isi_sleepstage{1}.isi{5}, ...
            stats{ipart}.isi_sleepstage{1}.isi{6}];


        mean_freq = [];
        stdev_freq = [];
        fano = [];
        CV = [];
        BI = [];
        X = [];
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            mean_freq   = [mean_freq; stats{ipart}.mean_freq{itemp}];
            stdev_freq  = [stdev_freq; stats{ipart}.stdev_freq{itemp}];
            fano        = [fano; stats{ipart}.fano_isi{itemp}];
            CV          = [CV; stats{ipart}.CV_isi{itemp}];
            BI          = [BI; stats{ipart}.burstindex{itemp}];
            X           = [X; 1,2,3,4,5,6];
        end

        colormap(C);

        subplot(3,4,9);
        b = bar(X',mean_freq','FaceColor','flat','Edgecolor','none');
        for k = 1:size(mean_freq,1)
            b(k).CData = k;
        end
        set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
        axis tight
        box off
        title('Freq');
        xlim([1.5,6.5]); % remove x, aka 'rest of data'

        subplot(3,4,10);
        b = bar(X',fano','FaceColor','flat','Edgecolor','none');
        for k = 1:size(fano',1)
            b(k).CData = k;
        end
        set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
        axis tight
        box off
        title('Fano');
        xlim([1.5,6.5]); % remove x, aka 'rest of data'

        subplot(3,4,11);
        b = bar(X',CV','FaceColor','flat','Edgecolor','none');
        for k = 1:size(CV',1)
            b(k).CData = k;
        end
        set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
        axis tight
        box off
        title('CV');
        xlim([1.5,6.5]); % remove x, aka 'rest of data'

        subplot(3,4,12);
        b = bar(X',BI','FaceColor','flat','Edgecolor','none');
        for k = 1:size(BI',1)
            b(k).CData = k;
        end
        set(gca,'Xtick', 1 : 6,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out');
        axis tight
        box off
        title('BI');
        xlim([1.5,6.5]); % remove x, aka 'rest of data'

        % print to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-all_clusters_hypnospikestats.pdf']),'-r600');
        close all
%
%         % Crosscorrelation between clusters per stage
%         for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'
%
%             cfgtemp             = [];
%             cfgtemp.trials      = find(SpikeTrials{ipart}.trialinfo.stage == istage);
%             cfgtemp.binsize     = 0.0005;
%             cfgtemp.maxlag      = 0.015;
%             cfgtemp.debias      = 'yes';
%             cfgtemp.method      = 'xcorr';
%             stats{ipart}.xcorr  = ft_spike_xcorr(cfgtemp,SpikeTrials{ipart});
%
%             fig = figure;
%             set(fig, 'units','normalized','position', [0 0 1 0.5]);
%             i = 1;
%             for ix = 1 : size(stats{ipart}.xcorr.xcorr,1)
%                 for iy = 1 : size(stats{ipart}.xcorr.xcorr,2)
%
%                     if ix > iy
%                         c = 'k';
%                     elseif ix < iy
%                         c = 'k';
%                     elseif ix == iy
%                         c = 'b';
%                     end
%
%                     x = stats{ipart}.xcorr.time;
%                     y = squeeze(stats{ipart}.xcorr.xcorr(ix,iy,:));
%                     if ~any(isnan(y))
%
%                         h = subplot(size(stats{ipart}.xcorr.xcorr,1),size(stats{ipart}.xcorr.xcorr,2),i); hold;
%
%                         bar(x,y,c,'EdgeColor','none','Barwidth',1);
%                         Lx = 1:length(x)/2;
%                         Rx = length(x)/2 : length(x);
%                         xintL = linspace(x(Lx(1)),x(Lx(end)),100)';
%                         yintL = spline(x(Lx),y(Lx),xintL);
%                         yintL = smooth1q(yintL,10);
%                         xintR = linspace(x(Rx(1)),x(Rx(end)),100)';
%                         yintR = spline(x(Rx),y(Rx),xintR);
%                         yintR = smooth1q(yintR,10);
%                         %                     plot(xintL,yintL,'r','linewidth',1);
%                         %                     plot(xintR,yintR,'r','linewidth',1);
%                         axis tight
%                         set(h,'yticklabel',{[]});
%                         set(h,'xticklabel',{[]});
%
%                     end
%                     pbaspect([1 1 1])
%                     i = i + 1;
%                 end
%             end
%
%             % print to file
%             fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
%             set(fig,'PaperOrientation','landscape');
%             set(fig,'PaperUnits','normalized');
%             set(fig,'PaperPosition', [0 0 1 1]);
%             print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-xcorr_SleepStage_',num2str(istage),'.pdf']),'-r600');
%             close all
%         end %istage

    end % ipart

    save(fname,'stats','-v7.3');

end % if file already exists
