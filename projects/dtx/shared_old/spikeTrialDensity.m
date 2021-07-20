function [stats] = spikeTrialDensity(cfg, SpikeTrials, force)

% SPIKERATESTATSEVENTS calculates spike statistics
%
% use as
%   [stats] = spikeTrialDensity(cfg, SpikeTrials, force)

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

cfg.spike.part_list     = ft_getopt(cfg.spike, 'part_list', 'all');

if strcmp(cfg.spike.part_list, 'all')
    cfg.spike.part_list = 1:size(cfg.directorylist, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'spikeTrialDensity.mat']);
if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname, 'stats');
    return
end

stats = {};

for ipart = cfg.spike.part_list
    
    if ipart > size(SpikeTrials, 2)
        continue
    end
    if isempty(SpikeTrials{ipart})
        continue
    end
    
    for markername = string(fields(SpikeTrials{ipart}))'
        
        if isempty(SpikeTrials{ipart}.(markername))
            continue
        end
        
        %spike psth (no smoothing)
        cfgtemp                         = [];
        if isfield(cfg.spike, 'psthbin')
            cfgtemp.binsize             = cfg.spike.psthbin.(markername);
        end
        cfgtemp.keeptrials              = 'yes';
        psth_event                      = ft_spike_psth(cfgtemp,SpikeTrials{ipart}.(markername));
        stats{ipart}.psth.(markername)  = psth_event;
        
        % spike density function, with smoothed version
        cfgtemp                         = [];
        cfgtemp.fsample                 = cfg.spike.resamplefs.(markername);   % sample at 1000 hz
        cfgtemp.keeptrials              = 'yes';
        if isfield(cfg.spike, 'sdftimwin')
            cfgtemp.timwin              = cfg.spike.sdftimwin.(markername);
        end
        sdf_line{ipart}                 = ft_spikedensity(cfgtemp, SpikeTrials{ipart}.(markername));
        
        % prepare data for stats on bar graph
        [n, e, b] = histcounts(sdf_line{ipart}.time, cfg.spike.nrsdfbins.(markername));
        binsize = diff(e);
        y = nan(size(sdf_line{ipart}.trial, 1), size(sdf_line{ipart}.label, 2), size(n, 2));
        x = e(1:end-1) + binsize/2;
        for itemp = 1 : size(sdf_line{ipart}.label, 2)
            for i = 1 : size(n, 2)
                for itrial = 1 : size(sdf_line{ipart}.trial, 1)
                    y(itrial, itemp, i) = mean(sdf_line{ipart}.trial(itrial, itemp, b == i), 3, 'omitnan');
                end
            end
        end
        
        sdf_bar{ipart}                      = [];
        sdf_bar{ipart}.trial                = y;
        sdf_bar{ipart}.time                 = x;
        sdf_bar{ipart}.label                = sdf_line{ipart}.label;
        sdf_bar{ipart}.dimord               = 'rpt_chan_time';
        sdf_bar{ipart}.avg                  = squeeze(mean(sdf_bar{ipart}.trial, 1, 'omitnan'));
        
        % contain spike density in stat data
        stats{ipart}.sdf_bar.(markername)       = sdf_bar{ipart};
        stats{ipart}.sdf_lin.(markername)       = sdf_line{ipart};
        
        % calculate baseline for dummy stats
        slim(1)                             = find(sdf_bar{ipart}.time > cfg.stats.bl.(markername)(1), 1, 'first');
        slim(2)                             = find(sdf_bar{ipart}.time < cfg.stats.bl.(markername)(2), 1, 'last');
        sdf_bar_bl{ipart}                   = sdf_bar{ipart};
        
        % baseline is mean during baseline period
        bl                                  = nanmean(sdf_bar{ipart}.trial(:, :, slim(1):slim(2)), 3);
        sdf_bar_bl{ipart}.trial             = ones(size(sdf_bar_bl{ipart}.trial)) .* bl;
        
        % clusterstats on bargraph, separate for each unit
        for itemp = 1 : size(sdf_bar{ipart}.label, 2)
            
            cfgtemp                                         = [];
            cfgtemp.channel                                 = itemp;
            cfgtemp.statistic                               = 'ft_statfun_depsamplesT';
            cfgtemp.alpha                                   = cfg.stats.alpha;
            cfgtemp.clusteralpha                            = 0.05;
            cfgtemp.method                                  = 'montecarlo';
            cfgtemp.computestat                             = 'yes';
            cfgtemp.correctm                                = 'cluster';
            cfgtemp.latency                                 = [cfg.stats.bl.(markername)(2) sdf_bar{ipart}.time(end)]; % active perio starts after baseline
            cfgtemp.ivar                                    = 1;
            cfgtemp.uvar                                    = 2;
            cfgtemp.design(1, :)                            = [ones(1, size(sdf_bar{ipart}.trial, 1)) ones(1, size(sdf_bar{ipart}.trial, 1)) *2];
            cfgtemp.design(2, :)                            = [1 : size(sdf_bar{ipart}.trial, 1) 1 : size(sdf_bar{ipart}.trial, 1)];
            cfgtemp.numrandomization                        = 1000;
            try
                stats{ipart}.stat.(markername){itemp}           = ft_timelockstatistics(cfgtemp, sdf_bar{ipart}, sdf_bar_bl{ipart});
                stats{ipart}.stat.(markername){itemp}.baseline  = bl;
            end
            
        end
        
    end
end

save(fname, 'stats', '-v7.3');
