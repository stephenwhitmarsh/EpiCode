function [stats] = spikeratestatsPSG(cfg, SpikeRaw, SpikeTrials, Hypnogram, force)

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

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikestatsSleepStages', cfg.circus.postfix, '.mat']);
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
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-ISI_all_data_sleepstage.pdf']),'-r600');
        close all

        % ISI per stage
        cfgtemp                                 = [];
        cfgtemp.trials                          = 'all';
        cfgtemp.outputunit                      = 'proportion';
        cfgtemp.bins                            = [0 : 0.0005 : 0.100]; %cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
        cfgtemp.param                           = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials                      = 'yes';
        stats{ipart}.isi_all_trials             = ft_spike_isi(cfgtemp,SpikeTrials{ipart});
        correct_stage_indx = 1 - min(unique(SpikeTrials{ipart}.trialinfo.stage)); %to make the first stage as the index 1

        % ISI descriptives per sleep stage
        for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'

            % ISI
            cfgtemp                                 = [];
            cfgtemp.trials                          = find(SpikeTrials{ipart}.trialinfo.stage == istage);
            cfgtemp.outputunit                      = 'proportion';
            cfgtemp.bins                            = [0 : 0.0005 : 0.100]; %cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
            cfgtemp.param                           = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
            stats{ipart}.isi_sleepstage{istage+correct_stage_indx}   = ft_spike_isi(cfgtemp,SpikeTrials{ipart});
            
            % create stats for each cluster x stage
            for itemp = 1 : length(SpikeTrials{ipart}.label)
                
                clear trialavg_isi
                isi_intraburst    = [];
                isi_interburst    = [];
                trials      = find(SpikeTrials{ipart}.trialinfo.stage == istage);
                i           = 1;
                isi_pooled  = [];

                for itrial = trials'

                    % get timings and ISIs per trial
                    indx            = SpikeTrials{ipart}.trial{itemp} == itrial;
                    t               = SpikeTrials{ipart}.time{itemp}(indx);
                    isi_all         = diff(t);

                    % counting bursts as in Colder et al. 1996, & Staba et al. 2002
                    indx            = isi_all < 0.01; % ask Stephane
                    burstindx       = zeros(size(indx));
                    toremove        = [];
                    
                    for burstlength = 1 : 10 %nr of consecutive ISI<0.01

                        pattern     = [false, true(1,burstlength), false];
                        bindx       = strfind(indx, pattern);

                        if ~isempty(bindx)
                            burstindx(bindx+1) = burstlength; % note +1 because pattern starts with zero
                            %fprintf('Found %d bursts of length %d in trial %d \n',length(bindx), burstlength, itrial);
                            
                            % add to list to correct for bursts
                            for ii = 1 : size(bindx,2)

                                % remove all but first spike (at +1)
                                toremove = [toremove, bindx(ii)+2:bindx(ii)+2+burstlength-1]; % burstlength = 1; 0 1 0 -> 0 1 x 0. %bindx begins by zero
                                
                                % add ISI within bursts
                                isi_intraburst = [isi_intraburst, isi_all(bindx(ii)+1:bindx(ii)+burstlength)]; %modif Paul : +1 because begins by zero, but then take all the isi of the burst (why was a -1 there ?)
                                
                            end
                        end
                        
                        stats{ipart}.burstsum{itemp}{istage+correct_stage_indx}(itrial, burstlength) = sum(length(bindx));  %+2 because some stages are -1

                    end

                    % concatinate ISIs, but only between bursts (not within)
                    t_interburst               = t(burstindx ~= 0);
                    isi_interburst             = [isi_interburst, diff(t_interburst)];

                    % remove subsequenct APs after first AP of a burst
                    t_corrected                 = t;
                    t_corrected(toremove)       = [];
                    isi_corrected               = diff(t_corrected);
                    
                    % basic descriptives
                    trialavg_isi(i)             = nanmean(isi_all);
                    trialfreq(i)                = 1/nanmean(isi_all);
                    spikecount(i)               = size(t,2);
                    spikecount_corrected(i)     = size(t_corrected,2);
                    
                    % according to Ponce-Alvarez, 2010
                    x                           = isi_corrected(1:end-1) ./ isi_corrected(2:end);
                    CV2_instant                 = 2 * abs(x - 1) ./ (x + 1);
                    CV2_trial(i)                = mean(CV2_instant);

                    x                           = isi_intraburst(1:end-1) ./ isi_intraburst(2:end);
                    CV2_intraburst_instant      = 2 * abs(x - 1) ./ (x + 1);
                    CV2_intraburst_trial(i)     = mean(CV2_intraburst_instant);

                    LV_instant                  = 3 * (x - 1).^2 ./ (x + 1).^2;
                    LV_trial(i)                 = mean(LV_instant);
                    IR_instant                  = abs(log(x));
                    IR_trial(i)                 = mean(IR_instant);
                    SI_instant                  = 0.5 * log((x+1).^2/(4*x));
                    SI_trial(i)                 = mean(SI_instant);

                    % concatinate ISIS over trials, corrected for bursts: for pooled CV
                    isi_pooled                  = [isi_pooled, isi_corrected];

                    % calculate CV per trial for averged CV
                    CV_trial(i)                 = nanstd(isi_corrected) / nanmean(isi_corrected);

                    % short vs long ISIs for BI
                    short(i)                    = sum(isi_all < 0.010);
                    long(i)                     = sum(isi_all < 0.100);

                    i = i + 1;
                end

                % get stats per sleep stage, over trials
                %modif temporaire Paul
                if strcmp(SpikeTrials{ipart}.label{itemp}, 'temp_8') && istage == 0
                    stats{ipart}.CV2_trialavg{itemp}(istage+correct_stage_indx)              = NaN;
                    stats{ipart}.CV2_intraburst_trialavg{itemp}(istage+correct_stage_indx)   = NaN;
                    stats{ipart}.isi_intraburst{itemp}{istage+correct_stage_indx}            = NaN;
                    stats{ipart}.isi_interburst{itemp}{istage+correct_stage_indx}            = NaN;
                    stats{ipart}.mean_freq{itemp}(istage+correct_stage_indx)                 = NaN;
                    stats{ipart}.stdev_freq{itemp}(istage+correct_stage_indx)                = NaN;
                    stats{ipart}.mean_isi{itemp}(istage+correct_stage_indx)                  = NaN;
                    stats{ipart}.burstindex{itemp}(istage+correct_stage_indx)                = NaN;
                    stats{ipart}.FF{itemp}(istage+correct_stage_indx)                        = NaN;
                    
                else
                    stats{ipart}.isi_intraburst{itemp}{istage+correct_stage_indx}            = isi_intraburst;
                    stats{ipart}.isi_interburst{itemp}{istage+correct_stage_indx}            = isi_interburst;
                    stats{ipart}.burst_trialsum{itemp}(istage+correct_stage_indx,:)          = sum(stats{ipart}.burstsum{itemp}{istage+correct_stage_indx});
                    stats{ipart}.mean_freq{itemp}(istage+correct_stage_indx)                 = nanmean(trialfreq);
                    stats{ipart}.stdev_freq{itemp}(istage+correct_stage_indx)                = nanstd(trialfreq);
                    [N,EDGES]                                               = histcounts(trialfreq,'BinWidth',0.5);
                    [M,I]                                                   = max(N);
                    stats{ipart}.mode_freq{itemp}(istage+correct_stage_indx)                 = mean(EDGES(I:I+1));
                    stats{ipart}.mean_isi{itemp}(istage+correct_stage_indx)                  = nanmean(trialavg_isi);
                    stats{ipart}.var_isi{itemp}(istage+correct_stage_indx)                   = nanstd(isi_pooled)^2;
                    stats{ipart}.burstindex{itemp}(istage+correct_stage_indx)                = sum(short) / sum(long);
                    stats{ipart}.FF{itemp}(istage+correct_stage_indx)                        = nanstd(spikecount_corrected)^2 / nanmean(spikecount_corrected);
                    stats{ipart}.spikecount{itemp}(istage+correct_stage_indx)                = sum(spikecount);
                    stats{ipart}.spikecount_corrected{itemp}(istage+correct_stage_indx)      = sum(spikecount_corrected);
                    stats{ipart}.CV_pooled{itemp}(istage+correct_stage_indx)                 = nanstd(isi_pooled)   / nanmean(isi_pooled);
                    stats{ipart}.CV_trialavg{itemp}(istage+correct_stage_indx)               = nanmean(CV_trial);
                    stats{ipart}.CV2_trialavg{itemp}(istage+correct_stage_indx)              = nanmean(CV2_trial);
                    stats{ipart}.CV2_intraburst_trialavg{itemp}(istage+correct_stage_indx)   = nanmean(CV2_intraburst_trial);
                    stats{ipart}.LV_trialavg{itemp}(istage+correct_stage_indx)               = nanmean(LV_trial);
                    stats{ipart}.IR_trialavg{itemp}(istage+correct_stage_indx)               = nanmean(IR_trial);
                    stats{ipart}.SI_trialavg{itemp}(istage+correct_stage_indx)               = nanmean(SI_trial);
                    stats{ipart}.nrwindows{itemp}(istage+correct_stage_indx)                 = sum(SpikeTrials{ipart}.trialinfo.stage == istage);
                end
            end
        end %istage
        
        close all
        
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            
            if ~isempty(Hypnogram)
                % plot hypnogram with spikerate for each cluster
                Hypnogram_part = Hypnogram(Hypnogram.part == ipart,:);
                
                
                fig = figure;
                C = linspecer(5); % colormap for nr. of sleep stages
                
                subplot(4,1,1); hold;
                
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
            end
            
            
            nr_stages = length(unique(SpikeTrials{ipart}.trialinfo.stage));
            fig = figure;
            C = linspecer(nr_stages); % colormap for nr. of sleep stages

            % calculate firingrate over time (here: trials)
            for itrial = SpikeTrials{ipart}.trialinfo.trialnr'
                indx    = SpikeTrials{ipart}.trial{itemp} == itrial;
                isi     = stats{ipart}.isi_all_trials.isi{itemp}(indx);
                stats{ipart}.trialavg_isi{itemp}(itrial) = nanmean(isi);
                stats{ipart}.trialavg_freq{itemp}(itrial) = 1 / nanmean(isi);
            end

            % plot spikerates
            subplot(4,1,2);
            x = (SpikeTrials{ipart}.trialinfo.starttime + (SpikeTrials{ipart}.trialinfo.endtime - SpikeTrials{ipart}.trialinfo.starttime)/2);
            plot(x,log(stats{ipart}.trialavg_freq{itemp}));
            %         plot(x,stats{ipart}.trialavg_freq{itemp}); %MODIF PAUL LOG
            %xlim(xl);
            %         ylabel('Firingrate(Hz)');%MODIF PAUL LOG
            ylabel('Log( Firingrate(Hz) )');
            title('Firingrate');
            
            %             % plot stats
            %             stageindx = [zeros(1,length(stats{ipart}.isi_sleepstage{1}.isi{1})) ...
            %                 ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{2})), ...
            %                 ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{3}))*2, ...
            %                 ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{4}))*3, ...
            %                 ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{5}))*4, ...
            %                 ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{6}))*5];
            %             stageISI = [stats{ipart}.isi_sleepstage{1}.isi{1}, ...
            %                 stats{ipart}.isi_sleepstage{1}.isi{2}, ...
            %                 stats{ipart}.isi_sleepstage{1}.isi{3}, ...
            %                 stats{ipart}.isi_sleepstage{1}.isi{4}, ...
            %                 stats{ipart}.isi_sleepstage{1}.isi{5}, ...
            %                 stats{ipart}.isi_sleepstage{1}.isi{6}];
            

            subplot(4,6,13);
            hold on
            errorbar(1:nr_stages,stats{ipart}.mean_freq{itemp},zeros(size(stats{ipart}.stdev_freq{itemp})),stats{ipart}.stdev_freq{itemp},'k','LineStyle','none');
            bar(1:nr_stages,stats{ipart}.mean_freq{itemp});
            %set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'x','W','1','2','3','R'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            %xlim([1.5,6.5]); % remove x, aka 'rest of data' %Paul set cfg
            title('Mean Freq+SD');
            
            %not for paul. add nr of seizures and length of pre inj in titles

            subplot(4,6,14);
            bar(1:5,stats{ipart}.nrwindows{itemp}(2:end));
            set(gca,'Xtick', 1 : 5,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
            axis tight
            box off
            title('Nr. windows');

            subplot(4,6,15);
            hold off
            plot(1:nr_stages,stats{ipart}.mode_freq{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('Mode Freq');

            subplot(4,6,16);
            plot(1:nr_stages,stats{ipart}.FF{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('Fano Factor');

            subplot(4,6,17);
            plot(1:nr_stages,stats{ipart}.CV_trialavg{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('CV trialavg');

            subplot(4,6,18);
            plot(1:nr_stages,stats{ipart}.CV2_trialavg{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('CV2 trialavg');

            subplot(4,6,19);
            plot(1:nr_stages,stats{ipart}.CV2_intraburst_trialavg{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('CV2 intraburst trialavg');

            subplot(4,6,20);
            plot(1:nr_stages,stats{ipart}.burstindex{itemp}(1:end),'.-');
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            axis tight
            box off
            title('Burstindex');

            subplot(4,6,21); hold;
            %Paul : revoir, pq que les bursts de 1 à 3 PA
            for iburstnr = 1:3
                plot(1:nr_stages,(stats{ipart}.burst_trialsum{itemp}(1:end,iburstnr) ./ stats{ipart}.spikecount{itemp}(1:end)') * 100','.-');
            end
            legend({'2','3','4'});
            set(gca,'Xtick', 1:nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out'); %Paul: ajouter cfg.stages.name
            ylabel('% bursts of total spikes');
            axis tight
            box off
            title('# bursts / # spikes');
            box off
            
            subplot(4,6,22); hold;
            %Paul à voir. PDF ?
            for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'+correct_stage_indx
                bar(stats{ipart}.isi_sleepstage{istage}.time*1000,stats{ipart}.isi_sleepstage{istage}.avg(itemp,:),1,'FaceColor',C(istage,:),'FaceAlpha',0.5);
            end
            [~, I] = max(stats{ipart}.isi_sleepstage{1}.avg(itemp,:));
            title(sprintf('ISI, max(W): %.1fms',stats{ipart}.isi_sleepstage{istage}.time(I)*1000));
            legend('-1','1','2','3','4','5','location','best','fontsize',4);

            ylabel('PDF');
            xlabel('ms');
            axis tight
            box off
            
            subplot(4,6,23); hold;
            for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'+correct_stage_indx
                histogram(stats{ipart}.isi_intraburst{itemp}{istage}*1000,'Normalization','pdf','FaceColor',C(istage,:),'EdgeColor','None','BinWidth',0.25);
            end
            [N,EDGES] = histcounts(stats{ipart}.isi_intraburst{itemp}{istage}*1000,1000);
            [~, I] = max(N);
            title(sprintf('ISI_intra, max(W): %.1fms', (EDGES(I)+EDGES(I+1))/2), 'Interpreter', 'none');
            axis tight
            xlabel('ms');
            ylabel('PDF');
            legend('-1','1','2','3','4','5','location','best','fontsize',4);

            box off
            
            subplot(4,6,24); hold;
            for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'+correct_stage_indx
                histogram(stats{ipart}.isi_interburst{itemp}{istage}*1000,1000,'Normalization','pdf','FaceColor',C(istage,:),'EdgeColor','None','BinWidth',10);
            end
            [N,EDGES] = histcounts(stats{ipart}.isi_interburst{itemp}{istage}*1000,1000);
            [~, I] = max(N);
            title(sprintf('ISI_inter, max(W): %.1fms', (EDGES(I)+EDGES(I+1))/2), 'Interpreter', 'none');
            legend('-1','1','2','3','4','5','location','best','fontsize',4);

            ylabel('PDF');
            xlabel('ms');
            xlim([0, 1000]);
            box off

            % print to file
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_',SpikeTrials{ipart}.label{itemp},'_stagespikestats.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_',SpikeTrials{ipart}.label{itemp},'_stagespikestats.png']),'-r600');
            close all
        end

        %% plot hypnogram with spikerate for all clusters at once, per night (part)
        %Hypnogram_part = Hypnogram(Hypnogram.part == ipart,:);
        
        fig = figure;
        %subplot(4,1,1); hold;
        %
        %     X = [];
        %     Y = [];
        %     for im = 1 : height(Hypnogram_part)
        %         if ~isempty(X)
        %             % if there's a gap, 'fill' with 0
        %             if Hypnogram_part.starttime(im) ~= X(end)
        %                 X = [X, X(end) Hypnogram_part.starttime(im)];
        %                 Y = [Y, 0,  0];
        %             end
        %         end
        %         X = [X, Hypnogram_part.starttime(im), Hypnogram_part.endtime(im)];
        %
        %         % height in hypnogram is based on order of config.hyp.contains
        %         switch cell2mat(Hypnogram_part.stagelabel(im))
        %             case 'AWAKE'
        %                 y = 5;
        %             case 'REM'
        %                 y = 4;
        %             case 'STAGE 1'
        %                 y = 3;
        %             case 'STAGE 2'
        %                 y = 2;
        %             case 'STAGE 3'
        %                 y = 1;
        %         end
        %         Y = [Y, y, y];
        %     end
        %
        %     for i = 1 : length(X)-1
        %         if Y(i) ~= 0 && Y(i+1) ~= 0
        %             if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
        %                 plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
        %             else
        %                 plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
        %             end
        %         end
        %     end
        %     set(gca,'Layer','top');
        %     set(gca,'Ytick', 1 : 5,'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'},'TickDir','out');
        %     title('Hypnogram');
        %     axis tight;
        %     xl = xlim;
        %     ylim([0 6]);
        %     legend('stages','location','eastoutside');
        %
      
        % plot spikerates
        subplot(4,1,2); hold;
        C = linspecer(length(SpikeTrials{ipart}.label));

        x = (SpikeTrials{ipart}.trialinfo.starttime + (SpikeTrials{ipart}.trialinfo.endtime - SpikeTrials{ipart}.trialinfo.starttime)/2);
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(x,log(stats{ipart}.trialavg_freq{itemp}),'color',C(itemp,:));
            %         plot(x,stats{ipart}.trialavg_freq{itemp},'color',C(itemp,:));%MODIF PAUL LOG
        end
        %xlim(xl);
        ylabel('Log( Firingrate(Hz) )');
        %     ylabel('Firingrate(Hz)');%MODIF PAUL LOG
        title('Firingrate');
        legend('location','eastoutside');
        
        %     % plot stats
        %     stageindx = [zeros(1,length(stats{ipart}.isi_sleepstage{1}.isi{1})) ...
        %         ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{2})), ...
        %         ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{3}))*2, ...
        %         ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{4}))*3, ...
        %         ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{5}))*4, ...
        %         ones(1, length(stats{ipart}.isi_sleepstage{1}.isi{6}))*5];
        %     stageISI = [stats{ipart}.isi_sleepstage{1}.isi{1}, ...
        %         stats{ipart}.isi_sleepstage{1}.isi{2}, ...
        %         stats{ipart}.isi_sleepstage{1}.isi{3}, ...
        %         stats{ipart}.isi_sleepstage{1}.isi{4}, ...
        %         stats{ipart}.isi_sleepstage{1}.isi{5}, ...
        %         stats{ipart}.isi_sleepstage{1}.isi{6}];
        %
     
        mean_freq = [];
        stdev_freq = [];
        fano = [];
        CV = [];
        BI = [];
        X = [];
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            mean_freq   = [mean_freq; stats{ipart}.mean_freq{itemp}];
            stdev_freq  = [stdev_freq; stats{ipart}.stdev_freq{itemp}];
            X           = [X; 1,2,3,4,5,6];
        end
        
        %Modif Paul Norm : scale firts instead of norm
        

        colormap(C);
        subplot(4,6,13); hold on
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.mean_freq{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('Freq');

        subplot(4,6,14); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.stdev_freq{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('SD Freq');

        subplot(4,6,15); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.mode_freq{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('Mode Freq');

        subplot(4,6,16); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.CV_pooled{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');

        ylabel('normalized');
        axis tight
        box off
        title('CV pooled');

        subplot(4,6,17); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.CV_trialavg{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');

        ylabel('normalized');
        axis tight
        box off
        title('CV trialavg');

        subplot(4,6,18); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.CV2_trialavg{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('CV2 trialavg');

        subplot(4,6,19); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(stats{ipart}.burstindex{itemp}(1:end),'scale','first'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('Burstindex');

        subplot(4,6,20); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,normalize(sum(stats{ipart}.burst_trialsum{itemp}(1:end,:),2) ./ stats{ipart}.spikecount_corrected{itemp}(1:end)'),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('# bursts / # spikes');

        % print to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-all_clusters_hypnospikestats_NORM.pdf']),'-r600');
        close all
        
        
        %Modif Paul : add plot without normalize
        fig = figure;
        
        colormap(C);
        subplot(4,6,13); hold on
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.mean_freq{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        ylim([0 30]); %PROVISOIRE PAUL
        box off
        title('Freq');
        
        subplot(4,6,14); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.stdev_freq{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('SD Freq');
        
        subplot(4,6,15); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.mode_freq{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('Mode Freq');
        
        %     subplot(4,6,16); hold on;
        %     for itemp = 1 : length(SpikeTrials{ipart}.label)
        %         plot(1:nr_stages,stats{ipart}.CV_pooled{itemp}(1:end),C(itemp,:));
        %     end
        %     set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        %     ylabel('normalized');
        %     axis tight
        %     box off
        %     title('CV pooled');
        
        subplot(4,6,17); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.CV_trialavg{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('CV trialavg');
        
        subplot(4,6,18); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.CV2_trialavg{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        %ylim([0.5 1.5]); %PAUL MODIF PROVISOIRE
        title('CV2 trialavg');
        
        subplot(4,6,19); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,stats{ipart}.burstindex{itemp}(1:end),'.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('Burstindex');
        
        subplot(4,6,20); hold on;
        for itemp = 1 : length(SpikeTrials{ipart}.label)
            plot(1:nr_stages,sum(stats{ipart}.burst_trialsum{itemp}(1:end,:),2) ./ stats{ipart}.spikecount_corrected{itemp}(1:end)','.-','color',C(itemp,:));
        end
        set(gca,'Xtick', 1 : nr_stages,'Xticklabels',{'-1','1','2','3','4','5'},'TickDir','out');
        ylabel('normalized');
        axis tight
        box off
        title('# bursts / # spikes');
        
        % print to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-all_clusters_hypnospikestats_notNORM.pdf']),'-r600');
        close all
        
        %
        %     %% Crosscorrelation between clusters per stage
        %     for istage = unique(SpikeTrials{ipart}.trialinfo.stage)'
        %
        %         cfgtemp             = [];
        %         cfgtemp.trials      = find(SpikeTrials{ipart}.trialinfo.stage == istage);
        %         cfgtemp.binsize     = 0.0005;
        %         cfgtemp.maxlag      = 0.015;
        %         cfgtemp.debias      = 'yes';
        %         cfgtemp.method      = 'xcorr';
        %         stats{ipart}.xcorr  = ft_spike_xcorr(cfgtemp,SpikeTrials{ipart});
        %
        %         fig = figure;
        %         set(fig, 'units','normalized','position', [0 0 1 0.5]);
        %         i = 1;
        %
        %         for ix = 1 : size(stats{ipart}.xcorr.xcorr,1)
        %             for iy = 1 : size(stats{ipart}.xcorr.xcorr,2)
        %
        %                 if ix > iy
        %                     c = 'k';
        %                 elseif ix < iy
        %                     c = 'k';
        %                 elseif ix == iy
        %                     c = 'b';
        %                 end
        %
        %                 x = stats{ipart}.xcorr.time;
        %                 y = squeeze(stats{ipart}.xcorr.xcorr(ix,iy,:));
        %                 if ~any(isnan(y))
        %
        %                     h = subplot(size(stats{ipart}.xcorr.xcorr,1),size(stats{ipart}.xcorr.xcorr,2),i); hold;
        %
        %                     bar(x,y,c,'EdgeColor','none','Barwidth',1);
        %                     Lx = 1:length(x)/2;
        %                     Rx = length(x)/2 : length(x);
        %                     xintL = linspace(x(Lx(1)),x(Lx(end)),100)';
        %                     yintL = spline(x(Lx),y(Lx),xintL);
        %                     yintL = smooth1q(yintL,10);
        %                     xintR = linspace(x(Rx(1)),x(Rx(end)),100)';
        %                     yintR = spline(x(Rx),y(Rx),xintR);
        %                     yintR = smooth1q(yintR,10);
        %                     %                     plot(xintL,yintL,'r','linewidth',1);
        %                     %                     plot(xintR,yintR,'r','linewidth',1);
        %                     axis tight
        %                     set(h,'yticklabel',{[]});
        %                     set(h,'xticklabel',{[]});
        %
        %                 else
        %                     test = test+1;
        %                 end
        %                 pbaspect([1 1 1])
        %                 i = i + 1;
        %             end
        %         end
        %
        %         % print to file
        %         fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        %         set(fig,'PaperOrientation','landscape');
        %         set(fig,'PaperUnits','normalized');
        %         set(fig,'PaperPosition', [0 0 1 1]);
        %         print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-xcorr_SleepStage_',num2str(istage),'.pdf']),'-r600');
        %         close all
        %     end %istage
  
    end % ipart

    save(fname,'stats','-v7.3');

end % if file already exists
