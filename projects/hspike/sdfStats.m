function [sdf] = sdfStats(cfg, SpikeTrials, force)

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'sdfStats.mat']);
if exist(fname, 'file') && force == false
    load(fname, 'sdf');
    return
end
    
for ipart = 1 : size(SpikeTrials, 2)
    
    if ipart > size(SpikeTrials,2)
        continue
    end
    if isempty(SpikeTrials{ipart})
        continue
    end
    for markername = string(fields(SpikeTrials{ipart}))'
        
        % spike density function, with smoothed version
        cfgtemp                         = [];
        cfgtemp.fsample                 = cfg.spike.resamplefs;   % sample at 1000 hz
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.latency                 = [cfg.spike.toi.(markername)(1), cfg.spike.toi.(markername)(2)];
        cfgtemp.timwin                  = [-0.05 0.05];
        sdf_orig{ipart}                 = ft_spikedensity(cfgtemp, SpikeTrials{ipart}.(markername));
        
        
%         figure;        
%         cfgtemp                 = [];
%         ft_spike_plot_raster(cfgtemp,  SpikeTrials{ipart}.(markername), sdf_orig{ipart});
%         
%         
        % prepare dummy data with baseline value per trial for stats
        slim(1)                         = find(sdf_orig{ipart}.time > cfg.stats.bl.(markername)(1), 1, 'first');
        slim(2)                         = find(sdf_orig{ipart}.time < cfg.stats.bl.(markername)(2), 1, 'last');
        sdf_bl{ipart}                   = sdf_orig{ipart};
        sdf_bl{ipart}.trial             = ones(size(sdf_orig{ipart}.trial)) .* nanmean(sdf_orig{ipart}.trial(:, :, slim(1):slim(2)), 3); % replace with mean
        
        % cluster stats per template
        for itemp = 1 : nrtemplates
            
            % statistics
            cfgtemp = [];
            cfgtemp.statistic                                   = 'ft_statfun_depsamplesT';
            cfgtemp.alpha                                       = cfg.stats.alpha;
            cfgtemp.clusteralpha                                = 0.01;
            cfgtemp.method                                      = 'montecarlo';
            cfgtemp.computestat                                 = 'yes';
            cfgtemp.correctm                                    = 'cluster';
            cfgtemp.latency                                     = [cfg.stats.bl.(markername)(2) sdf_orig{ipart}.time(end)]; % active period starts after baseline
            cfgtemp.ivar                                        = 1;
            cfgtemp.uvar                                        = 2;
            cfgtemp.design(1, :)                                 = [ones(1, size(sdf_smooth{ipart}.trial, 1)) ones(1, size(sdf_bl{ipart}.trial, 1)) *2];
            cfgtemp.design(2, :)                                 = [1 : size(sdf_smooth{ipart}.trial, 1) 1 : size(sdf_bl{ipart}.trial, 1)];
            cfgtemp.numrandomization                            = 100;
            cfgtemp.channel                                     = itemp;
            stats_firingrate{ipart}.(markername).clusterstat{itemp} = ft_timelockstatistics(cfgtemp, sdf_smooth{ipart}, sdf_bl{ipart});
            
            % calculate baseline
            slim(1) = find(sdf_orig{ipart}.time > cfg.stats.(markername).bl(1), 1, 'first');
            slim(2) = find(sdf_orig{ipart}.time < cfg.stats.(markername).bl(2), 1, 'last');
            stats{ipart}.(markername).clusterstat{itemp}.bl.avg        = nanmean(sdf_orig{ipart}.avg(itemp, slim(1):slim(2)), 2);
            stats_firingrate{ipart}.(markername).clusterstat{itemp}.bl.var        = nanmean(sdf_orig{ipart}.var(itemp, slim(1):slim(2)), 2);
            stats_firingrate{ipart}.(markername).clusterstat{itemp}.bl.dof        = nanmean(sdf_orig{ipart}.dof(itemp, slim(1):slim(2)), 2);
            stats_firingrate{ipart}.(markername).clusterstat{itemp}.bl.trialavg   = nanmean(sdf_orig{ipart}.trial(:, itemp, slim(1):slim(2)), 3);
        end
        
    end
    
end