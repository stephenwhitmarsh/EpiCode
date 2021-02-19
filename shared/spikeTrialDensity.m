function [stats] = spikeTrialDensity(cfg, SpikeTrials, force)

% SPIKERATESTATSEVENTS calculates spike statistics
%
% use as
%   [stats] = spikeTrialDensity(cfg, SpikeTrials, force)
% or, if need to only load precomputed data :
%   [stats] = spikeTrialDensity(cfg)
% 
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

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'spikeTrialDensity.mat']);

if nargin == 1
    if exist(fname, 'file')
        fprintf('Reading %s\n', fname);
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname, 'stats');
            catch ME
                err_count = err_count + 1;
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        stats = {};
        return
    end
end

if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    
    count = 0;
    err_count = 0;
    while count == err_count
        try
            load(fname, 'stats');
        catch ME
            err_count = err_count + 1;
        end
        count = count + 1;
    end
    return
end

cfg.spike.part_list     = ft_getopt(cfg.spike, 'part_list', 'all');

if strcmp(cfg.spike.part_list, 'all')
    cfg.spike.part_list = 1:size(cfg.directorylist, 2);
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
        
        %spike pssth (no smoothing)
        cfgtemp                         = [];
        if isfield(cfg.spike, 'psthbin')
            cfgtemp.binsize             = cfg.spike.psthbin.(markername);
        end
        cfgtemp.keeptrials              = 'yes';
%         cfgtemp.trials                  = ~SpikeTrials{ipart}.(markername).trialinfo.artefact;        
        psth_event                      = ft_spike_psth(cfgtemp,SpikeTrials{ipart}.(markername));
        stats{ipart}.psth.(markername)  = psth_event;
        
        % spike density function, with smoothed version
        cfgtemp                         = [];
        cfgtemp.fsample                 = cfg.spike.resamplefs.(markername);   % sample at 1000 hz
%         cfgtemp.trials                  = ~SpikeTrials{ipart}.(markername).trialinfo.artefact;
        cfgtemp.keeptrials              = 'yes';
        if isfield(cfg.spike, 'sdftimwin')
            cfgtemp.timwin              = cfg.spike.sdftimwin.(markername);
        end
        sdf_line{ipart}                 = ft_spikedensity(cfgtemp, SpikeTrials{ipart}.(markername));
        
        % prepare data for stats on bar graph
        [n, e, b] = histcounts(sdf_line{ipart}.time, cfg.spike.nrsdfbins);
        binsize = diff(e);
        y = nan(size(sdf_line{ipart}.trial, 1), size(sdf_line{ipart}.label, 2), size(n, 2));
        x = e(1:end-1) + binsize/2;
        for itemp = 1 : size(sdf_line{ipart}.label, 2)
            for i = 1 : size(n, 2)
                for itrial = 1 : size(sdf_line{ipart}.trial, 1)
                    y(itrial, itemp, i) = mean(sdf_line{ipart}.trial(itrial, itemp, b == i));
                end
            end
        end
        
        sdf_bar{ipart}                      = [];
        sdf_bar{ipart}.trial                = y;
        sdf_bar{ipart}.time                 = x;
        sdf_bar{ipart}.label                = sdf_line{ipart}.label;
        sdf_bar{ipart}.dimord               = 'rpt_chan_time';
        sdf_bar{ipart}.avg                  = squeeze(mean(sdf_bar{ipart}.trial, 1));
        
        % contain spike density in stat data
        stats{ipart}.sdf_bar.(markername)   = sdf_bar{ipart};
        stats{ipart}.sdf_lin.(markername)   = sdf_line{ipart};
        
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
            
            % do stats on data without artefacts
            cleanindx                                       = ~SpikeTrials{ipart}.(markername).trialinfo.artefact;
            trialcount                                      = sum(cleanindx);
            cfgtemp.design(:, 1)                            = [ones(1, trialcount) ones(1, trialcount) * 2];
            cfgtemp.design(:, 2)                            = [1 : trialcount 1 : trialcount];
            cfgtemp.numrandomization                        = 1000;
            dat_sel                                         = sdf_bar{ipart};
            dat_sel.trial                                   = dat_sel.trial(cleanindx, :, :);
            dat_bl_sel                                      = sdf_bar_bl{ipart};
            dat_bl_sel.trial                                = sdf_bar_bl{ipart}.trial(cleanindx, :, :);

            stats{ipart}.stat.(markername){itemp}           = ft_timelockstatistics(cfgtemp, dat_sel, dat_bl_sel);
            stats{ipart}.stat.(markername){itemp}.baseline  = bl;
            clear dat_sel dat_bl_sel
        end 
    end
end

save(fname, 'stats', '-v7.3');
