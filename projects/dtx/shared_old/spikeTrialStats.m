function [stats] = spikeTrialStats(cfg, SpikeTrials, force, postfix)

% SPIKERATESTATSEVENTS calculates and plots spike statistics
%
% use as
%   [stats] = spikeTrialStats(cfg, SpikeTrials, force, postfix)
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

if nargin < 4
    postfix = [];
end
cfg.spike.RPV = ft_getopt(cfg.spike, 'RPV', 0.001);

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'spikestats-', postfix, '.mat']);
if exist(fname, 'file') && force == false
    fprintf('Reading %s\n', fname);
    load(fname, 'stats');
    return
end

hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

for ipart = 1 : size(SpikeTrials, 2)
    
    for markername = string(fields(SpikeTrials{ipart}))'
        
        if isempty(SpikeTrials{ipart}.(markername)); continue; end
        
        if isfield(SpikeTrials{ipart}.(markername).trialinfo, 'hyplabel')
            SpikeTrials{ipart}.(markername).trialinfo.hyplabel(SpikeTrials{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
        end
        
        % ISI over conditions
        cfgtemp                                     = [];
        cfgtemp.outputunit                          = 'spikecount';
        cfgtemp.bins                                =  cfg.spike.ISIbins;%0 : 0.0005 : 0.200;   % use bins of 0.5 milliseconds
        cfgtemp.param                               = 'coeffvar';         % compute the coefficient of variation (sd/mn of isi)
        isi_temp                                    = ft_spike_isi(cfgtemp, SpikeTrials{ipart}.(markername));
        
%         % Xcorr over conditions
%         cfgtemp                                     = [];
%         cfgtemp.binsize                             = 0.001; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
%         cfgtemp.maxlag                              = 0.200;
%         cfgtemp.outputunit                          = 'proportion';        
%         xcorr_temp                                  = ft_spike_xcorr(cfgtemp, SpikeTrials{ipart}.(markername));
        
        % because cfg.keeptrials doesn't work:
        if  isfield(SpikeTrials{ipart}.(markername).trialinfo, 'hyplabel')
            for hyplabel = hyplabels
                
                trials = SpikeTrials{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                
                % ISI per condition
                cfgtemp                                     = [];
                cfgtemp.outputunit                          = 'proportion';
                cfgtemp.bins                                = 0 : 0.0005 : 0.200; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
                cfgtemp.param                               = 'coeffvar';         % compute the coefficient of variation (sd/mn of isi)
                cfgtemp.trials                              = trials;
                isi_hyp_temp.(hyplabel)                     = ft_spike_isi(cfgtemp, SpikeTrials{ipart}.(markername));
                
                %             % Xcorr per condition
                %             cfgtemp                                     = [];
                %             cfgtemp.binsize                             = 0.001; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
                %             cfgtemp.maxlag                              = 0.200;
                %             cfgtemp.trials                              = trials;
                %             xcorr_hyp_temp.(hyplabel)                   = ft_spike_xcorr(cfgtemp, SpikeTrials{ipart}.(markername));
            end
        end
        
        for itemp = 1 : size(SpikeTrials{ipart}.(markername).label, 2)
            
            stats{ipart}.(markername){itemp}.isi            = isi_temp.isi{itemp};
            stats{ipart}.(markername){itemp}.isi_avg        = isi_temp.avg(itemp, :);
            stats{ipart}.(markername){itemp}.isi_avg_time   = isi_temp.time;
            stats{ipart}.(markername){itemp}.label          = isi_temp.label{itemp};
            refr = sum(stats{ipart}.(markername){itemp}.isi < cfg.spike.RPV);
            tot  = length((stats{ipart}.(markername){itemp}.isi));
            stats{ipart}.(markername){itemp}.RPV = refr/tot;
            
            % spike autocorr
            cfgtemp              = [];
            cfgtemp.spikechannel = SpikeTrials{ipart}.(markername).label{itemp};
            spike_temp           = ft_spike_select(cfgtemp,SpikeTrials{ipart}.(markername));
            stats{ipart}.(markername){itemp}.autocorr   = ft_spike_xcorr(cfgtemp, spike_temp);
            stats{ipart}.(markername){itemp}.autocorr   = rmfield(stats{ipart}.(markername){itemp}.autocorr, 'cfg');
            %             stats{ipart}.(markername){itemp}.xcorr          = xcorr_temp.xcorr;
            %             stats{ipart}.(markername){itemp}.xcorr_time     = xcorr_temp.time;
            
            if isfield(SpikeTrials{ipart}.(markername).trialinfo, 'hyplabel')
                for hyplabel = hyplabels
                    stats{ipart}.(markername){itemp}.(hyplabel).isi            = isi_hyp_temp.(hyplabel).isi{itemp};
                    stats{ipart}.(markername){itemp}.(hyplabel).isi_avg        = isi_hyp_temp.(hyplabel).avg(itemp, :);
                    stats{ipart}.(markername){itemp}.(hyplabel).isi_avg_time   = isi_hyp_temp.(hyplabel).time;
                    stats{ipart}.(markername){itemp}.(hyplabel).label          = isi_hyp_temp.(hyplabel).label{itemp};
                    %                 stats{ipart}.(markername){itemp}.(hyplabel).xcorr          = xcorr_hyp_temp.(hyplabel).xcorr;
                    %                 stats{ipart}.(markername){itemp}.(hyplabel).xcorr_time     = xcorr_hyp_temp.(hyplabel).time;
                end
            end
            
            % find bursts
            isi_intraburst = [];
            isi_interburst = [];
            
            % add trialinfo so that trials can later be selected
            stats{ipart}.(markername){itemp}.trialinfo = SpikeTrials{ipart}.(markername).trialinfo;
            
            for itrial = 1 : size(SpikeTrials{ipart}.(markername).trialinfo, 1)
                
                % get timings and ISIs per trial
                indx            = SpikeTrials{ipart}.(markername).trial{itemp} == itrial;
                t               = SpikeTrials{ipart}.(markername).time{itemp}(indx);
                isi_all         = diff(t);
                amps            = SpikeTrials{ipart}.(markername).amplitude{itemp}(indx);
                
                % counting bursts as in Colder et al. 1996, & Staba et al. 2002
                indx            = isi_all < 0.03; %modif paul set 30 ms as in intracellular data
                burstindx       = zeros(size(indx));
                toremove        = [];
                for burstlength = 1 : 10
                    pattern     = [false, true(1, burstlength), false];
                    bindx       = strfind(indx, pattern);
                    if ~isempty(bindx)
                        burstindx(bindx+1) = burstlength; % note +1 index because pattern starts with zero
                        
                        % add to list to correct for bursts
                        for ii = 1 : size(bindx, 2)
                            
                            % remove all but first spike (at +1)
                            toremove = [toremove, bindx(ii) + 2 : bindx(ii) + 2 + burstlength - 1]; % burstlength = 1; 0 1 0 -> 0 1 x 0
                            
                            % add ISI within bursts
                            isi_intraburst = [isi_intraburst, isi_all(bindx(ii)+1:bindx(ii)+burstlength-1)];
                            
                        end
                    end
                    stats{ipart}.(markername){itemp}.burstsum(itrial, burstlength) = length(bindx);
                end
                stats{ipart}.(markername){itemp}.burst_trialsum(itrial)            = sum(stats{ipart}.(markername){itemp}.burstsum(itrial, :));
                
                % concatinate ISIs, but only between bursts (not within)
                stats{ipart}.(markername){itemp}.t_interburst{itrial}              = t(burstindx ~= 0);
                stats{ipart}.(markername){itemp}.isi_interburst{itrial}            = diff(stats{ipart}.(markername){itemp}.t_interburst{itrial});
                
                % remove subsequenct APs after first AP of a burst
                stats{ipart}.(markername){itemp}.t_corrected{itrial}               = t;
                stats{ipart}.(markername){itemp}.t_corrected{itrial}(toremove)     = [];
                stats{ipart}.(markername){itemp}.isi_corrected{itrial}             = diff(stats{ipart}.(markername){itemp}.t_corrected{itrial});
                
                % basic descriptives
                stats{ipart}.(markername){itemp}.trialavg_isi(itrial)              = nanmean(isi_all);
                stats{ipart}.(markername){itemp}.trialfreq(itrial)                 = 1/nanmean(isi_all);
                stats{ipart}.(markername){itemp}.trialfreq_corrected(itrial)       = 1/nanmean(stats{ipart}.(markername){itemp}.isi_corrected{itrial});                
                stats{ipart}.(markername){itemp}.spikecount(itrial)                = size(t, 2);
                stats{ipart}.(markername){itemp}.spikecount_corrected(itrial)      = size(stats{ipart}.(markername){itemp}.t_corrected{itrial}, 2);
                
                % according to Ponce-Alvarez, 2010
                x                                                                  = stats{ipart}.(markername){itemp}.isi_corrected{itrial}(1:end-1) ./ stats{ipart}.(markername){itemp}.isi_corrected{itrial}(2:end);
                CV2_instant                                                        = 2 * abs(x - 1) ./ (x + 1);
                stats{ipart}.(markername){itemp}.CV2_trial(itrial)                 = nanmean(CV2_instant);
                x                                                                  = isi_intraburst(1:end-1) ./ isi_intraburst(2:end);
                CV2_intraburst_instant                                             = 2 * abs(x - 1) ./ (x + 1);
                stats{ipart}.(markername){itemp}.CV2_intraburst_trial(itrial)      = nanmean(CV2_intraburst_instant);
                LV_instant                                                         = 3 * (x - 1).^2 ./ (x + 1).^2;
                stats{ipart}.(markername){itemp}.LV_trial(itrial)                  = nanmean(LV_instant);
                IR_instant                                                         = abs(log(x));
                stats{ipart}.(markername){itemp}.IR_trial(itrial)                  = nanmean(IR_instant);
                SI_instant                                                         = 0.5 * log((x+1).^2/(4*x));
                stats{ipart}.(markername){itemp}.SI_trial(itrial)                  = nanmean(SI_instant);
                
                % calculate CV per trial for averged CV
                stats{ipart}.(markername){itemp}.CV_trial(itrial)                  = nanstd(stats{ipart}.(markername){itemp}.isi_corrected{itrial}) / nanmean(stats{ipart}.(markername){itemp}.isi_corrected{itrial});
                
                % short vs long ISIs for BI
                stats{ipart}.(markername){itemp}.short(itrial)                     = sum(isi_all < 0.010);
                stats{ipart}.(markername){itemp}.long(itrial)                      = sum(isi_all < 0.100);
                
                stats{ipart}.(markername){itemp}.amplitude(itrial)                 = nanmean(amps);
                
                if isnan(nanmean(amps))       
                    stats{ipart}.(markername){itemp}.burst_trialsum(itrial) = nan;
                end
            end
        end
    end
end

save(fname, 'stats', '-v7.3');

