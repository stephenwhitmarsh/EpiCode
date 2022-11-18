function [stats] = spikePSTH(cfg, SpikeTrials, force, LFP)

% SPIKERATESTATSEVENTS calculates spike statistics
%
% use as
%   [stats] = spikeTrialDensity(cfg, SpikeTrials, force)
% or, if need to only load precomputed data :
%   [stats] = spikeTrialDensity(cfg)
% 
% cfg.spike.psthbin.(markername) : bin size of the psth in seconds
% cfg.spike.psthoutputunit : output of the psth, either 'rate' (default), 
% 'spikecount' or 'proportion' (see ft_spike_psth for detailed documentation)
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
                disp('Something went wrong loading the file. Trying again...')    
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
            disp('Something went wrong loading the file. Trying again...')            
        end
        count = count + 1;
    end
    return
end

cfg.spike.part_list      = ft_getopt(cfg.spike, 'part_list', 'all');
cfg.spike.psthoutputunit = ft_getopt(cfg.spike, 'psthoutputunit', 'rate');

if strcmp(cfg.spike.part_list, 'all')
    cfg.spike.part_list = 1:size(cfg.directorylist, 2);
end
stats = {};

if nargin < 4
    % read LFP for correlation later on
    cfgtemp             = cfg;
    cfgtemp.LFP.name    = cfg.spike.name;
    LFP                 = readLFP(cfgtemp);
end

for ipart = cfg.spike.part_list
    
    if ipart > size(SpikeTrials, 2)
        continue
    end
    if isempty(SpikeTrials{ipart})
        continue
    end
    
    for markername = string(cfg.spike.name)
               
        if isempty(SpikeTrials{ipart}.(markername))
            continue
        end
        
        % spike peristimulus time histogram
        cfgtemp                         = [];
        if isfield(cfg.spike, 'psthbin')
            cfgtemp.binsize             = cfg.spike.psthbin.(markername);
        end
        cfgtemp.outputunit              = cfg.spike.psthoutputunit;
        cfgtemp.keeptrials              = 'yes';
        psth_event                      = ft_spike_psth(cfgtemp, SpikeTrials{ipart}.(markername));
        stats{ipart}.psth.(markername)  = psth_event;
        
%         % spike density function NOT USED ANYMORE
%         cfgtemp                         = [];
%         cfgtemp.fsample                 = cfg.spike.resamplefs.(markername);
%         cfgtemp.keeptrials              = 'yes';
%         if isfield(cfg.spike, 'sdftimwin')
%             cfgtemp.timwin              = cfg.spike.sdftimwin.(markername);
%         end
%         stats{ipart}.sdf.(markername)   = ft_spikedensity(cfgtemp, SpikeTrials{ipart}.(markername));
        
        % prepare data for stats on bar graph
%         [n, e, b] = histcounts(stats{ipart}.sdf_lin.(markername).time, cfg.spike.nrsdfbins);
%         binsize = diff(e);
%         y = nan(size(stats{ipart}.sdf_lin.(markername).trial, 1), size(stats{ipart}.sdf_lin.(markername).label, 2), size(n, 2));
%         x = e(1:end-1) + binsize/2;
%         for itemp = 1 : size(stats{ipart}.sdf_lin.(markername).label, 2)
%             for i = 1 : size(n, 2)
%                 for itrial = 1 : size(stats{ipart}.sdf_lin.(markername).trial, 1)
%                     y(itrial, itemp, i) = mean(stats{ipart}.sdf_lin.(markername).trial(itrial, itemp, b == i));
%                 end
%             end
%         end
%         stats{ipart}.sdf_bar.(markername).trial     = y;
%         stats{ipart}.sdf_bar.(markername).time      = x;
%         stats{ipart}.sdf_bar.(markername).label     = stats{ipart}.sdf_lin.(markername).label;
%         stats{ipart}.sdf_bar.(markername).dimord    = 'rpt_chan_time';
%         stats{ipart}.sdf_bar.(markername).avg       = squeeze(mean(stats{ipart}.sdf_bar.(markername).trial, 1));
%           

        % calculate baseline for dummy stats
        slim(1)                             = find(stats{ipart}.psth.(markername).time > cfg.stats.bl.(markername)(1), 1, 'first');
        slim(2)                             = find(stats{ipart}.psth.(markername).time < cfg.stats.bl.(markername)(2), 1, 'last');
        psth_bl                             = stats{ipart}.psth.(markername);
        
        % baseline is mean during baseline period
        bl                                  = nanmean(stats{ipart}.psth.(markername).trial(:, :, slim(1):slim(2)), 3);
        psth_bl.trial                       = ones(size(psth_bl.trial)) .* bl;
        
        % clusterstats
        for itemp = 1 : size(stats{ipart}.psth.(markername).label, 2)
            
            cfgtemp                                         = [];
            cfgtemp.channel                                 = itemp;
            cfgtemp.statistic                               = 'ft_statfun_depsamplesT';
            cfgtemp.alpha                                   = cfg.stats.alpha;
            cfgtemp.clusteralpha                            = 0.01; %0.05
            cfgtemp.method                                  = 'montecarlo';
            cfgtemp.computestat                             = 'yes';
            cfgtemp.correctm                                = 'cluster';
            cfgtemp.latency                                 = [cfg.stats.bl.(markername)(2) stats{ipart}.psth.(markername).time(end)]; % active period starts after baseline
            cfgtemp.ivar                                    = 1;
            cfgtemp.uvar                                    = 2;
            
            % do stats on data without artefacts
            cleanindx                                       = ~SpikeTrials{ipart}.(markername).trialinfo.BAD_cnt;
            trialcount                                      = sum(cleanindx);
            
            if trialcount <= 1
                stats{ipart}.stat.(markername){itemp}       = [];
                fprintf('Not enough trials in %s in part %s for statistics', markername, ipart);
                continue
            end
            
            cfgtemp.design(:, 1)                            = [ones(1, trialcount) ones(1, trialcount) * 2];
            cfgtemp.design(:, 2)                            = [1 : trialcount 1 : trialcount];
            cfgtemp.numrandomization                        = 10000;
            dat_sel                                         = stats{ipart}.psth.(markername);
            dat_sel.trial                                   = dat_sel.trial(cleanindx, :, :);
            dat_bl_sel                                      = psth_bl;
            dat_bl_sel.trial                                = psth_bl.trial(cleanindx, :, :);

            stats{ipart}.stat.(markername){itemp}           = ft_timelockstatistics(cfgtemp, dat_sel, dat_bl_sel);
            stats{ipart}.stat.(markername){itemp}.baseline  = bl;

            % note if unit responds statistically
            stats{ipart}.stat.(markername){itemp}.responsive_pos   = false;
            stats{ipart}.stat.(markername){itemp}.responsive_neg   = false;
            stats{ipart}.stat.(markername){itemp}.responsive       = false;
            
            % TODO: simpify (lis this really needed?)
            if isfield(stats{ipart}.stat.(markername){itemp}, 'posclusters')
                for ipos = 1 : size(stats{ipart}.stat.(markername){itemp}.posclusters, 2)
                    if stats{ipart}.stat.(markername){itemp}.posclusters(ipos).prob < cfg.stats.alpha
                        stats{ipart}.stat.(markername){itemp}.responsive_pos = true;
                        stats{ipart}.stat.(markername){itemp}.responsive = true;
                    end
                end
            end
            if isfield(stats{ipart}.stat.(markername){itemp}, 'negclusters')
                for ineg = 1 : size(stats{ipart}.stat.(markername){itemp}.negclusters, 2)
                    if stats{ipart}.stat.(markername){itemp}.negclusters(ineg).prob < cfg.stats.alpha
                        stats{ipart}.stat.(markername){itemp}.responsive_neg = true;
                        stats{ipart}.stat.(markername){itemp}.responsive = true;
                    end
                end
            end
            
            % TODO: define latency independently
            % correlate firing rate with LFP and pick largest abs(rho)
            if ~isempty(LFP{ipart}.(markername))
                LFP_avg = ft_timelockanalysis([], LFP{ipart}.(markername));
                
                % resample to same time-axis
                cfgtemp         = [];
                cfgtemp.time{1} = stats{ipart}.psth.(markername).time;
                LFP_ds          = ft_resampledata(cfgtemp, LFP_avg);
                
                [corr_rho, corr_pval] = ...
                    corr(LFP_ds.avg(:, :)', stats{ipart}.psth.(markername).avg', 'type', 'pearson');
                               
                for iunit = 1 : size(stats{ipart}.psth.(markername).label, 2)
                    stats{ipart}.psth.(markername).corr_chan{iunit}     = LFP_ds.label;
                    stats{ipart}.psth.(markername).corr_rho(iunit, :)   = corr_rho(:, iunit);
                    stats{ipart}.psth.(markername).corr_pval(iunit, :)  = corr_pval(:, iunit);
                end
            end
            clear dat_sel dat_bl_sel
        end % itemp
    end % markername
end % ipart

save(fname, 'stats', '-v7.3');
