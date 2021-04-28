
function [stats, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(cfg,SpikeRaw,SpikeTrials,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [stats, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(cfg,SpikeRaw,SpikeTrials,force)
% 
% Creates spike stats and figures based on raw and epoched spike data 
%
% Necessary input:
% 
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files 
% cfg.circus.channel    = micro electrode names
%
% The following can be retreaved with readSpykingCircus.m:
%
% SpikeRaw = raw spike data in FieldTrip raw spike data structure
% SpikeTrials = spike data epoched in FieldTrip trial data structure
%
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fname = fullfile(cfg.datasavedir,[cfg.prefix,'all_data_spikedata_stats.mat']);
if exist(fname,'file') && force == false
    load(fname,'stats','stats_bar','sdf_orig_out','sdf_bar_out');
else
    
    temp                = dir(fullfile(cfg.datasavedir,[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.ncs']));
    hdr_fname           = fullfile(temp(1).folder,temp(1).name);
    hdr                 = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % redefine trials to 1-second windows for ISI
    cfgtemp             = [];
    cfgtemp.trl         = (1 : hdr.Fs : hdr.nSamples)';
    cfgtemp.trl(:,2)    = cfgtemp.trl(:,1) + hdr.Fs;
    cfgtemp.trl(:,3)    = zeros(size(cfgtemp.trl,1),1);
    cfgtemp.trl         = cfgtemp.trl(1:end-1,:);
    cfgtemp.trlunit     = 'samples';
    cfgtemp.hdr         = hdr;
    spiketrials_1s      = ft_spike_maketrials(cfgtemp,SpikeRaw);
    
    % ISI over 1-second windows
    cfgtemp             = [];
    cfgtemp.outputunit  = 'spikecount';
    cfgtemp.bins        = [0:0.0005:0.025];   % use bins of 0.5 milliseconds
    cfgtemp.param       = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
    stats.isi_1s        = ft_spike_isi(cfgtemp,spiketrials_1s);
       
    % plot ISI for each templates of 1-sec windows`
    nrtemplates =  length(SpikeRaw.label);
    
    fig = figure; hold;
    for itemp = 1 : nrtemplates
        % plot ISI fo 1-second windows
        subplot(round(nrtemplates/2+0.5),2,itemp);
        bar(stats.isi_1s.time*1000,stats.isi_1s.avg(itemp,:),1);
        [y,indx] = max(stats.isi_1s.avg(itemp,:));
        title(sprintf('Unit: %d, Max ISI: %.1fms',itemp,stats.isi_1s.time(indx)*1000));
        xlabel('ms');
        xticks(stats.isi_1s.time*1000);
        xtickangle(90);
        axis tight
        grid on
        set(gca,'fontsize',6);
    end
    
    % print ISI to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'ISI_1s_windows.pdf']),'-r600');
    
    % ISI as continuous data
    for itemp = 1 : nrtemplates
        stats.isi{itemp} = diff(SpikeRaw.sample{itemp}) / hdr.Fs * 1000;
    end
    
    % plot ISI histogram
    fig = figure; hold;
    for itemp = 1 : nrtemplates
        subplot(round(nrtemplates/2+0.25),2,itemp);
        histogram(stats.isi{itemp},'BinWidth',0.5,'BinLimits',[0,25]);
        xlabel('ms');
        xtickangle(90);
        axis tight
        grid on
        set(gca,'fontsize',6);
    end
    
    % print ISI to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'ISI_own.pdf']),'-r600');
    
    % cross-correlation between template (over 1 second trials)
    cfgtemp             = [];
    cfgtemp.binsize     = 0.0005;
    cfgtemp.maxlag      = 0.015;
    cfgtemp.debias      = 'yes';
    cfgtemp.method      = 'xcorr';
    stats.xcorr         = ft_spike_xcorr(cfgtemp,spiketrials_1s);
    
    % plot xcorr
    fig = figure;
    set(fig, 'units','normalized','position', [0 0 1 0.5]);
    i = 1;
    for ix = 1 : size(stats.xcorr.xcorr,1)
        for iy = 1 : size(stats.xcorr.xcorr,2)
            
            if ix > iy
                c = [0 0 0];
            end
            
            if ix < iy
                c = [0 0 0];
            end
            
            if ix == iy
                c = [0 0 0.8];
            end
            
            x = stats.xcorr.time;
            y = squeeze(stats.xcorr.xcorr(ix,iy,:));
            if ~any(isnan(y))
                
                h = subplot(size(stats.xcorr.xcorr,1),size(stats.xcorr.xcorr,2),i);
                hold;
                
                Lx = 1:length(x)/2;
                Rx = length(x)/2 : length(x);
                
                xintL = linspace(x(Lx(1)),x(Lx(end)),100)';
                yintL = spline(x(Lx),y(Lx),xintL);
                yintL = smooth(yintL,10);
                
                xintR = linspace(x(Rx(1)),x(Rx(end)),100)';
                yintR = spline(x(Rx),y(Rx),xintR);
                yintR = smooth(yintR,10);
                
                
                bar(x,y);
                plot(xintL,yintL,'r','linewidth',1);
                plot(xintR,yintR,'r','linewidth',1);
                axis tight
                ax = axis;
                ylim([0,ax(4)]);
                set(h,'yticklabel',{[]});
                t = sprintf('%dx%d',ix,iy);
                title(t);
                pbaspect([1 1 1])
                grid on
            end
            i = i + 1;
        end
    end
    
    % print to file
    fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'xcorr_1s_windows.pdf']),'-r600');

    % stats per pattern
    for ilabel = 1 : size(SpikeTrials,2)
        
        % ISI over 1-second windows
        cfgtemp                         = [];
        cfgtemp.outputunit              = 'proportion';
        cfgtemp.bins                    = [0:0.0005:0.025];   % use bins of 0.5 milliseconds
        cfgtemp.param                   = 'coeffvar';       % compute the coefficient of variation (sd/mn of isis)
        cfgtemp.keeptrials              = 'yes';
        stats.isi_pattern_all{ilabel}   = ft_spike_isi(cfgtemp,SpikeTrials{ilabel});
        
        cfgtemp.latency                 = [-2 -0.15];
        stats.isi_pattern_bl{ilabel}    = ft_spike_isi(cfgtemp,SpikeTrials{ilabel});
        
        cfgtemp.latency                 = [-0.15 0.15];
        stats.isi_pattern_ac{ilabel}    = ft_spike_isi(cfgtemp,SpikeTrials{ilabel});
        
        % plot ISI per channel, per latency, per template
        fig = figure; hold;
        for itemp = 1 : nrtemplates
            
            subplot(nrtemplates,3,(itemp-1)*3+1);
            bar(stats.isi_pattern_all{ilabel}.time*1000,stats.isi_pattern_all{ilabel}.avg(itemp,:),1);
            [y,indx] = max(stats.isi_pattern_all{ilabel}.avg(itemp,:));
            title(sprintf('Unit: %d, Max ISI: all%.1fms',itemp,stats.isi_pattern_all{ilabel}.time(indx)*1000));
            xlabel('ms');
            xticks(stats.isi_pattern_all{ilabel}.time*1000);
            xtickangle(90);
            axis tight
            grid on
            set(gca,'fontsize',6);
            
            subplot(nrtemplates,3,(itemp-1)*3+2);
            bar(stats.isi_pattern_bl{ilabel}.time*1000,stats.isi_pattern_bl{ilabel}.avg(itemp,:),1);
            [y,indx] = max(stats.isi_pattern_bl{ilabel}.avg(itemp,:));
            title(sprintf('Unit: %d, Max ISI baseline: %.1fms',itemp,stats.isi_pattern_bl{ilabel}.time(indx)*1000));
            xlabel('ms');
            xticks(stats.isi_pattern_bl{ilabel}.time*1000);
            xtickangle(90);
            axis tight
            grid on
            set(gca,'fontsize',6);
            
            subplot(nrtemplates,3,(itemp-1)*3+3);
            bar(stats.isi_pattern_ac{ilabel}.time*1000,stats.isi_pattern_ac{ilabel}.avg(itemp,:),1);
            [y,indx] = max(stats.isi_pattern_ac{ilabel}.avg(itemp,:));
            title(sprintf('Unit: %d, Max ISI active: %.1fms',itemp,stats.isi_pattern_ac{ilabel}.time(indx)*1000));
            xlabel('ms');
            xticks(stats.isi_pattern_ac{ilabel}.time*1000);
            xtickangle(90);
            axis tight
            grid on
            set(gca,'fontsize',6);
        end
        
        % print ISIs to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'ISI_',cfg.name{ilabel},'.pdf']),'-r600');
        
        % spike density function, with smoothed version
        cfgtemp                         = [];
        cfgtemp.fsample                 = cfg.spike.resamplefs;   % sample at 1000 hz
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.latency                 = [cfg.epoch.toi{ilabel}(1), cfg.epoch.toi{ilabel}(2)];
        cfgtemp.timwin                  = cfg.spike.toispikerate{ilabel} / 4;
        sdf_orig                        = ft_spikedensity(cfgtemp,SpikeTrials{ilabel});
        cfgtemp.timwin                  = cfg.spike.toispikerate{ilabel};
        sdf_smooth                      = ft_spikedensity(cfgtemp,SpikeTrials{ilabel});
        
        % prepare dummy data with baseline value per trial for stats
        slim(1)                         = find(sdf_orig.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
        slim(2)                         = find(sdf_orig.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
        sdf_bl                          = sdf_orig;
        sdf_bl.trial                    = ones(size(sdf_orig.trial)) .* nanmean(sdf_orig.trial(:,:,slim(1):slim(2)),3); % replace with mean
        
        % cluster stats per template
        for itemp = 1 : nrtemplates
            
            % statistics
            cfgtemp = [];
            cfgtemp.statistic                   = 'ft_statfun_depsamplesT';
            cfgtemp.alpha                       = cfg.stats.alpha;
            cfgtemp.clusteralpha                = 0.01;
            cfgtemp.method                      = 'montecarlo';
            cfgtemp.computestat                 = 'yes';
            cfgtemp.correctm                    = 'cluster';
            cfgtemp.latency                     = [cfg.stats.bltoi{ilabel}(2) sdf_orig.time(end)]; % active period starts after baseline
            cfgtemp.ivar                        = 1;
            cfgtemp.uvar                        = 2;
            cfgtemp.design(1,:)                 = [ones(1,size(sdf_orig.trial,1)) ones(1,size(sdf_orig.trial,1)) *2];
            cfgtemp.design(2,:)                 = [1 : size(sdf_orig.trial,1) 1 : size(sdf_orig.trial,1)];
            cfgtemp.numrandomization            = 1000;
            cfgtemp.channel                     = itemp;
            stats.clusterstat{ilabel}{itemp}    = ft_timelockstatistics(cfgtemp,sdf_smooth,sdf_bl);
            
            % calculate baseline
            slim(1) = find(sdf_orig.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
            slim(2) = find(sdf_orig.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
            stats.clusterstat{ilabel}{itemp}.bl.avg        = nanmean(sdf_orig.avg(itemp,slim(1):slim(2)),2);
            stats.clusterstat{ilabel}{itemp}.bl.var        = nanmean(sdf_orig.var(itemp,slim(1):slim(2)),2);
            stats.clusterstat{ilabel}{itemp}.bl.dof        = nanmean(sdf_orig.dof(itemp,slim(1):slim(2)),2);
            stats.clusterstat{ilabel}{itemp}.bl.trialavg   = nanmean(sdf_orig.trial(:,itemp,slim(1):slim(2)),3);
        end
        
        % plot cluster stats per patterns
        for itemp = 1 : nrtemplates
            
            clear x y si sel lag
            fig             = figure;
            fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
            % plot positive clusters
            
            subplot(3,3,[1 2]); hold;
            
            % plot firing rate as line
            plot(sdf_smooth.time,sdf_smooth.avg(itemp,:),'b');
            plot(sdf_orig.time,sdf_orig.avg(itemp,:),'k');
            
            lag = -(sdf_orig.time(1) - cfg.stats.bltoi{ilabel}(2)) * cfg.spike.resamplefs;
            
            if isfield(stats.clusterstat{ilabel}{itemp},'posclusters')
                for ipos = 1 : size(stats.clusterstat{ilabel}{itemp}.posclusters,2)
                    if stats.clusterstat{ilabel}{itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        sel = find(stats.clusterstat{ilabel}{itemp}.posclusterslabelmat == ipos);
                        plot(stats.clusterstat{ilabel}{itemp}.time(sel),sdf_smooth.avg(itemp,sel+lag),'g','linewidth',2); % smoothed
                        si = sel+lag;
                        [Y,I] = max(sdf_smooth.avg(itemp,sel+lag));
                        
                        x = sdf_smooth.time(si(I));
                        y = sdf_smooth.avg(itemp,si(I));
                        plot(x,y+0.25,'v','markersize',10,'color',[0 1 0],'markerfacecolor',[0 1 0]);
                        d = (sdf_smooth.avg(itemp,si(I)) / stats.clusterstat{ilabel}{itemp}.bl.avg) * 100;
                        text(x+0.05,y+0.25,sprintf('%.1f%%\n',d),'HorizontalAlignment','left','VerticalAlignment','middle');
                    end
                end
            end
            
            % plot negative clusters
            if isfield(stats.clusterstat{ilabel}{itemp},'negclusters')
                for ineg = 1 : size(stats.clusterstat{ilabel}{itemp}.negclusters,2)
                    if stats.clusterstat{ilabel}{itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        sel = find(stats.clusterstat{ilabel}{itemp}.negclusterslabelmat == ineg);
                        plot(stats.clusterstat{ilabel}{itemp}.time(sel),sdf_smooth.avg(itemp,sel+lag),'r','linewidth',2); % smoothed
                        si = sel+lag;
                        [Y,I] = min(sdf_smooth.avg(itemp,sel+lag));
                        x = sdf_smooth.time(si(I));
                        y = sdf_smooth.avg(itemp,si(I));
                        plot(x,y-0.25,'^','markersize',10,'color',[1 0 0],'markerfacecolor',[1 0 0]);
                        d = (sdf_smooth.avg(itemp,si(I)) / stats.clusterstat{ilabel}{itemp}.bl.avg) * 100;
                        text(x+0.05,y-0.25,sprintf('%.1f%%\n',d),'HorizontalAlignment','left','VerticalAlignment','middle');
                    end
                end
            end
            
            % plot baseline
            plot(cfg.epoch.toi{ilabel},[stats.clusterstat{ilabel}{itemp}.bl.avg, stats.clusterstat{ilabel}{itemp}.bl.avg],':k');
            
            % determine zero crossing
            zci             = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
            indx            = zci(sdf_smooth.avg(itemp,:) - stats.clusterstat{ilabel}{itemp}.bl.avg);
            times           = sdf_smooth.time(indx);
            remove          = find(diff(times) < 0.15); % minimal duration
            times(remove)   = [];
            indx(remove)    = [];
            
            % plot durations between zero crossings
            dtimes = diff(times);
            plot(times,sdf_smooth.avg(itemp,indx),'.','markersize',10,'color',[0 0 1],'markerfacecolor',[1 0 0]);
            for itxt = 1 : size(dtimes,2)
                x = (times(itxt) + times(itxt+1))/2;
                y = stats.clusterstat{ilabel}{itemp}.bl.avg;
                text(x,y,sprintf('%.0fms',dtimes(itxt)*1000),'HorizontalAlignment','center','VerticalAlignment','bottom','color',[0 0 1],'fontsize',6);
            end
            axis tight;
            
            % plot baseline patch
            x = [cfg.stats.bltoi{ilabel}(1) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(1)];
            ax = axis;
            y = [ax(3) ax(3) ax(4) ax(4)];
            patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
            
            % plot baseline text
            x = (cfg.stats.bltoi{ilabel}(1) + cfg.stats.bltoi{ilabel}(2))/2;
            y = ax(4)*0.9;
            d = stats.clusterstat{ilabel}{itemp}.bl.avg;
            text(x,y,sprintf('%.1f p/s',d),'HorizontalAlignment','center');
            xlabel('Time (sec)');
            ylabel('Spikerate (Hz)');
            
            % plot firing rate as line
            plot(sdf_smooth.time,sdf_smooth.avg(itemp,:),'b');
            plot(sdf_orig.time,sdf_orig.avg(itemp,:),'k');
            
            % prepare data for stats on bar graph
            [n,e,b] = histcounts(sdf_orig.time,200);
            binsize = diff(e);
            x = e(1:end-1) + binsize/2;
            for i = 1 : size(n,2)
                for itrial = 1 : size(sdf_orig.trial,1)
                    y(itrial,i) = mean(sdf_orig.trial(itrial,itemp,b == i));
                end
            end
            
            sel                     = find(~any(isnan(y)));
            sdf_bar                 = [];
            sdf_bar.trial(:,1,:)    = y(:,sel);
            sdf_bar.time            = x(sel);
            sdf_bar.label{1}        = sdf_smooth.label{itemp};
            sdf_bar.dimord          = 'rpt_chan_time';
            sdf_bar.avg             = squeeze(mean(sdf_bar.trial,1));
            
            sdf_bar_out{ilabel}{itemp} = sdf_bar;
            sdf_orig_out{ilabel}{itemp} = sdf_orig;
            
            % calculate baseline for dummy stats
            slim(1)             = find(sdf_bar.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
            slim(2)             = find(sdf_bar.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
            sdf_bar_bl          = sdf_bar;
            sdf_bar_bl.trial    = ones(size(sdf_bar_bl.trial)) .* nanmean(sdf_bar.trial(:,slim(1):slim(2)),2); % replace with mean
            
            % clusterstats on bargraph
            cfgtemp = [];
            cfgtemp.statistic                       = 'ft_statfun_depsamplesT';
            cfgtemp.alpha                           = cfg.stats.alpha;
            cfgtemp.clusteralpha                    = 0.1;
            cfgtemp.method                          = 'montecarlo';
            cfgtemp.computestat                     = 'yes';
            cfgtemp.correctm                        = 'cluster';
            cfgtemp.latency                         = [cfg.stats.bltoi{ilabel}(2) sdf_orig.time(end)]; % active perio starts after baseline
            cfgtemp.ivar                            = 1;
            cfgtemp.uvar                            = 2;
            cfgtemp.design(1,:)                     = [ones(1,size(sdf_bar.trial,1)) ones(1,size(sdf_bar.trial,1)) *2];
            cfgtemp.design(2,:)                     = [1 : size(sdf_bar.trial,1) 1 : size(sdf_bar.trial,1)];
            cfgtemp.numrandomization                = 1000;
            cfgtemp.channel                         = 1;
            stats_bar.clusterstat{ilabel}{itemp}    = ft_timelockstatistics(cfgtemp,sdf_bar,sdf_bar_bl);
            
            % calculate baseline
            slim(1) = find(sdf_bar.time > cfg.stats.bltoi{ilabel}(1), 1, 'first');
            slim(2) = find(sdf_bar.time < cfg.stats.bltoi{ilabel}(2), 1, 'last');
            stats_bar.clusterstat{ilabel}{itemp}.bl.avg        = nanmean(sdf_bar.avg(slim(1):slim(2)));
            stats_bar.clusterstat{ilabel}{itemp}.bl.trialavg   = nanmean(sdf_bar.trial(:,1,slim(1):slim(2)),3);
            
            % plot firing rate as bar graph
            subplot(3,3,[4 5]); hold;
            clear y
            [n,e,b] = histcounts(sdf_orig.time,200);
            for i = 1 : size(n,2)
                y(i) = mean(sdf_orig.avg(itemp,b == i));
            end
            binsize = diff(e);
            x = e(1:end-1) + binsize/2;
            bar(x,y,1,'edgecolor','none','facecolor',[0 0 0]);
            axis tight;
            hold on
            
            % plot positive clusters
            lag = size(sdf_bar.avg,1) - size(stats_bar.clusterstat{ilabel}{itemp}.mask,2);
            
            if isfield(stats_bar.clusterstat{ilabel}{itemp},'posclusters')
                for ipos = 1 : size(stats_bar.clusterstat{ilabel}{itemp}.posclusters,2)
                    if stats_bar.clusterstat{ilabel}{itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        sel = find(stats_bar.clusterstat{ilabel}{itemp}.posclusterslabelmat == ipos);
                        bar(stats_bar.clusterstat{ilabel}{itemp}.time(sel),sdf_bar.avg(sel+lag),1,'facecolor','g','edgecolor','g');
                        
                        % plot percentage
                        si = sel+lag;
                        [Y,I] = max(sdf_bar.avg(sel+lag));
                        x = sdf_bar.time(si(I));
                        y = sdf_bar.avg(si(I));
                        d = (sdf_bar.avg(si(I)) / stats_bar.clusterstat{ilabel}{itemp}.bl.avg) * 100 - 100;
                        stats_bar.clusterstat{ilabel}{itemp}.maxcluster.perc{ipos} = d;
                        stats_bar.clusterstat{ilabel}{itemp}.maxcluster.x{ipos} = x;
                        stats_bar.clusterstat{ilabel}{itemp}.maxcluster.y{ipos} = y;
                        text(x+0.05,y+0.25,sprintf('+%.1f%%\n',d),'HorizontalAlignment','center','VerticalAlignment','middle');  
                    end
                end
            end
            
            % plot negative clusters
            if isfield(stats_bar.clusterstat{ilabel}{itemp},'negclusters')
                for ineg = 1 : size(stats_bar.clusterstat{ilabel}{itemp}.negclusters,2)
                    if stats_bar.clusterstat{ilabel}{itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        sel = find(stats_bar.clusterstat{ilabel}{itemp}.negclusterslabelmat == ineg);
                        bar(stats_bar.clusterstat{ilabel}{itemp}.time(sel),sdf_bar.avg(sel+lag),1,'facecolor','r','edgecolor','r');
                        
                        % plot percentage
                        si = sel+lag;
                        [Y,I] = min(sdf_bar.avg(sel+lag));
                        x = sdf_bar.time(si(I));
                        y = sdf_bar.avg(si(I));
                        d = 100 - (sdf_bar.avg(si(I)) / stats_bar.clusterstat{ilabel}{itemp}.bl.avg) * 100;
                        stats_bar.clusterstat{ilabel}{itemp}.mincluster.perc{ineg} = d;
                        stats_bar.clusterstat{ilabel}{itemp}.mincluster.x{ineg} = x;
                        stats_bar.clusterstat{ilabel}{itemp}.mincluster.y{ineg} = y;                        
                        text(x+0.05,y-0.25,sprintf('-%.1f%%\n',d),'HorizontalAlignment','center','VerticalAlignment','middle');
                    end
                end
            end
            
            % plot baseline patch
            x = [cfg.stats.bltoi{ilabel}(1) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(2) cfg.stats.bltoi{ilabel}(1)];
            ax = axis;
            y = [ax(3) ax(3) ax(4) ax(4)];
            patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
            
            % plot baseline
            plot(cfg.epoch.toi{ilabel},[stats.clusterstat{ilabel}{itemp}.bl.avg, stats.clusterstat{ilabel}{itemp}.bl.avg],':k');
            
            % plot raster
            subplot(3,3,[7,8]); hold;
            cfgtemp                 = [];
            cfgtemp.spikechannel    = itemp;
            cfgtemp.latency         = [cfg.epoch.toi{ilabel}(1), cfg.epoch.toi{ilabel}(2)];
            cfgtemp.trialborders    = 'no';
            ft_spike_plot_raster(cfgtemp,SpikeTrials{ilabel});
            
            % plot ISI patch in raster
            ax = axis;
            patch([-2 -0.15 -0.15 -2],[ax(3) ax(3) ax(4) ax(4)],'b','facealpha',0.2,'edgecolor','none');
            patch([-0.15 0.15 0.15 -0.15],[ax(3) ax(3) ax(4) ax(4)],'r','facealpha',0.2,'edgecolor','none');
            
            % plot ISI
            subplot(3,3,3); hold;
            barh = bar(stats.isi_1s.time*1000,[stats.isi_1s.avg(itemp,:); stats.isi_pattern_all{ilabel}.avg(itemp,:); stats.isi_pattern_bl{ilabel}.avg(itemp,:); stats.isi_pattern_ac{ilabel}.avg(itemp,:);]',1,'grouped','edgecolor','none');
            barh(1).FaceColor = [0 0 0];
            barh(2).FaceColor = [0 1 0];
            barh(3).FaceColor = [0 0 1];
            barh(4).FaceColor = [1 0 0];
            legend({'1 sec all data','wh0le trial','baseline trial','active trial'});
            xticks(stats.isi_pattern_ac{ilabel}.time*1000);
            xtickangle(90);
            axis tight
            grid on
            set(gca,'fontsize',4);
            ylabel('Count');
            xlabel('ISI (ms)');
            
            % peak width accoridng to Gast et. al
            subplot(3,3,6); hold;
            temp        = dir(fullfile(cfg.datasavedir,[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.ncs']));
            hdr_fname   = fullfile(temp(1).folder,temp(1).name);
            hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
            tempsel = squeeze(SpikeRaw.template(itemp,SpikeRaw.template_maxchan(itemp),:));
            temptime = ((1:size(SpikeRaw.template,3))/hdr.Fs*1000)';
            
            % interpolate template
            temptime_int = linspace(temptime(1),temptime(end),10000);
            tempsel_int = pchip(temptime,tempsel,temptime_int);
            plot(temptime_int,tempsel_int,'k');
            
            axis tight
            [Ypos,Xpos] = findpeaks(tempsel_int,temptime_int,'NPeaks',1,'SortStr','descend');
            [Yneg,Xneg] = findpeaks(-tempsel_int,temptime_int,'NPeaks',2,'SortStr','descend','MinPeakDistance',0.5);
            plot([Xpos,Xneg(1)],[Yneg(1)*0.3, Yneg(1)*0.3],'-o','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
            plot([Xpos,Xneg(2)],-[Yneg(2)*0.3, Yneg(2)*0.3],'-o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
            
            x = (Xpos + Xneg(1))/2;
            y = Yneg(1)*0.1;
            text(x,y,sprintf('%.0fms',abs(Xpos-Xneg(1))*1000),'HorizontalAlignment','center');
            
            x = (Xpos + Xneg(2))/2;
            y = -Yneg(2)*0.1;
            text(x,y,sprintf('%.0fms',abs(Xpos-Xneg(2))*1000),'HorizontalAlignment','center');
            xlabel('time')
            
            midline = min(tempsel_int) + (max(tempsel_int) - min(tempsel_int)) / 2 ;
            indx = zci(tempsel_int - midline);
            plot(temptime_int(indx),[midline, midline],'-o','Color',[0 1 0],'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
            
            x = sum(temptime_int(indx))/length(indx);
            y = midline*1.1;
            text(x,y,sprintf('%.0fms',(temptime_int(indx(2))-temptime_int(indx(1)))*1000),'HorizontalAlignment','center');
            
            % print to file
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'spikerates_template_',num2str(itemp),'_pattern_',cfg.name{ilabel},'.pdf']),'-r600');
            set(fig,'PaperOrientation','portrait');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'spikerates_template_',num2str(itemp),'_pattern_',cfg.name{ilabel},'.png']),'-r600');
            close all
        end % itemplate
        
    end % ilabel
    
    save(fname,'stats','stats_bar','sdf_orig_out','sdf_bar_out','-v7.3');
    
end % if file already exists
