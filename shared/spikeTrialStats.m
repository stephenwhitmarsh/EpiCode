function [stats] = spikeTrialStats(cfg, SpikeTrials, force)

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

cfg.spike.RPV       = ft_getopt(cfg.spike, 'RPV', 0.002);
cfg.spike.postfix   = ft_getopt(cfg.spike, 'postfix', '');
fname               = fullfile(cfg.datasavedir, [cfg.prefix, 'spikestats', cfg.spike.postfix, '.mat']);

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
    load(fname, 'stats');
    return
end

hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

for ipart = 1 : size(SpikeTrials, 2)
    
    for markername = string(cfg.spike.name)
        
        if isempty(SpikeTrials{ipart}.(markername)); continue; end
        
        if isfield(SpikeTrials{ipart}.(markername).trialinfo, 'hyplabel')
            SpikeTrials{ipart}.(markername).trialinfo.hyplabel(SpikeTrials{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
        end
        
        % ISIs independent of hyplabels
        cfgtemp                                     = [];
        cfgtemp.outputunit                          = 'spikecount';
        cfgtemp.bins                                = cfg.spike.ISIbins;%0 : 0.0005 : 0.200;   % use bins of 0.5 milliseconds
        cfgtemp.param                               = 'coeffvar';         % compute the coefficient of variation (sd/mn of isi)
        isi_temp_all                                = ft_spike_isi(cfgtemp, SpikeTrials{ipart}.(markername));
        
        %         % Xcorr over conditions
        %         cfgtemp                                     = [];
        %         cfgtemp.binsize                             = 0.001; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
        %         cfgtemp.maxlag                              = 0.200;
        %         cfgtemp.outputunit                          = 'proportion';
        %         cfgtemp.channelcmb                          = [SpikeTrials{ipart}.(markername).label', SpikeTrials{ipart}.(markername).label'];
        %         xcorr_temp                                  = ft_spike_xcorr(cfgtemp, SpikeTrials{ipart}.(markername));
        
        % ISIs dependent of hyplabels
        if any(contains(SpikeTrials{ipart}.(markername).trialinfo.hyplabel, hyplabels))
            for hyplabel = hyplabels
                
                trials = SpikeTrials{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                
                if ~any(trials)
                    continue
                end
                
                % ISI per condition
                cfgtemp                                     = [];
                cfgtemp.outputunit                          = 'proportion';
                cfgtemp.bins                                = 0 : 0.0005 : 0.200; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
                cfgtemp.param                               = 'coeffvar';         % compute the coefficient of variation (sd/mn of isi)
                cfgtemp.trials                              = trials;
                isi_temp_label.(hyplabel)                   = ft_spike_isi(cfgtemp, SpikeTrials{ipart}.(markername));
                
                % Xcorr per condition
                %                 cfgtemp                                     = [];
                %                 cfgtemp.binsize                             = 0.001; % cfg.spike.ISIbins;   % use bins of 0.5 milliseconds
                %                 cfgtemp.maxlag                              = 0.200;
                %                 cfgtemp.trials                              = trials;
                %                 cfgtemp.channelcmb                          = [SpikeTrials{ipart}.(markername).label', SpikeTrials{ipart}.(markername).label'];
                %                 xcorr_hyp_temp.(hyplabel)                   = ft_spike_xcorr(cfgtemp, SpikeTrials{ipart}.(markername));
            end
        end
        
        for itemp = 1 : size(SpikeTrials{ipart}.(markername).label, 2)
            
            stats{ipart}.(markername){itemp}.isi            = isi_temp_all.isi{itemp};
            stats{ipart}.(markername){itemp}.isi_avg        = isi_temp_all.avg(itemp, :);
            stats{ipart}.(markername){itemp}.isi_avg_time   = isi_temp_all.time;
            stats{ipart}.(markername){itemp}.label          = isi_temp_all.label{itemp};
            
            % percentage refractory period violations
            refr = sum(    stats{ipart}.(markername){itemp}.isi < cfg.spike.RPV);
            tot  = length((stats{ipart}.(markername){itemp}.isi));
            stats{ipart}.(markername){itemp}.RPV = refr/tot;
            
            % spike autocorr
            %             cfgtemp              = [];
            %             cfgtemp.spikechannel = SpikeTrials{ipart}.(markername).label{itemp};
            %             spike_temp           = ft_spike_select(cfgtemp,SpikeTrials{ipart}.(markername));
            %             stats{ipart}.(markername){itemp}.autocorr   = ft_spike_xcorr(cfgtemp, spike_temp);
            %             stats{ipart}.(markername){itemp}.autocorr   = rmfield(stats{ipart}.(markername){itemp}.autocorr, 'cfg');
            %             stats{ipart}.(markername){itemp}.xcorr      = xcorr_temp.xcorr;
            %             stats{ipart}.(markername){itemp}.xcorr_time = xcorr_temp.time;
            
            if any(contains(SpikeTrials{ipart}.(markername).trialinfo.hyplabel, hyplabels))
                for hyplabel = hyplabels
                    if ~isfield(isi_temp_label, hyplabel)
                        continue
                    end
                    stats{ipart}.(markername){itemp}.(hyplabel).isi            = isi_temp_label.(hyplabel).isi{itemp};
                    stats{ipart}.(markername){itemp}.(hyplabel).isi_avg        = isi_temp_label.(hyplabel).avg(itemp, :);
                    stats{ipart}.(markername){itemp}.(hyplabel).isi_avg_time   = isi_temp_label.(hyplabel).time;
                    stats{ipart}.(markername){itemp}.(hyplabel).label          = isi_temp_label.(hyplabel).label{itemp};
                    %                     stats{ipart}.(markername){itemp}.(hyplabel).xcorr          = xcorr_hyp_temp.(hyplabel).xcorr;
                    %                     stats{ipart}.(markername){itemp}.(hyplabel).xcorr_time     = xcorr_hyp_temp.(hyplabel).time;
                end
            end
            
            ft_progress('init','text', fprintf('Starting on "%s", part %d of %d, unit %d of %d', markername, ipart, size(SpikeTrials, 2), itemp, size(SpikeTrials{ipart}.(markername).label, 2)));
              
            % add trialinfo so that trials can later be selected
            stats{ipart}.(markername){itemp}.trialinfo = SpikeTrials{ipart}.(markername).trialinfo;
            
            for itrial = 1 : size(SpikeTrials{ipart}.(markername).trialinfo, 1)
                
                ft_progress(itrial / size(SpikeTrials{ipart}.(markername).trialinfo, 1), 'Computing stats for trial %d of %d', itrial, size(SpikeTrials{ipart}.(markername).trialinfo, 1));

                % get timings and ISIs per trial
                indx            = SpikeTrials{ipart}.(markername).trial{itemp} == itrial;
                trial_length    = SpikeTrials{ipart}.(markername).trialtime(itrial, 2) - SpikeTrials{ipart}.(markername).trialtime(itrial, 1);
                t               = SpikeTrials{ipart}.(markername).time{itemp}(indx);
                isi_all         = diff(t);
                amps            = SpikeTrials{ipart}.(markername).amplitude{itemp}(indx);
                
                % counting bursts as in Colder et al. 1996, & Staba et al. 2002
                indx            = isi_all < 0.01;
                burstindx       = zeros(size(indx));
                toremove        = [];
                isi_intraburst  = []; 
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
                stats{ipart}.(markername){itemp}.trialfreq(itrial)                 = length(t) / trial_length;%1/nanmean(isi_all);
                stats{ipart}.(markername){itemp}.trialfreq_corrected(itrial)       = 1/nanmean(stats{ipart}.(markername){itemp}.isi_corrected{itrial});
                stats{ipart}.(markername){itemp}.spikecount(itrial)                = size(t, 2);
                stats{ipart}.(markername){itemp}.spikecount_corrected(itrial)      = size(stats{ipart}.(markername){itemp}.t_corrected{itrial}, 2);
                
                % regularity metrics according to Ponce-Alvarez, 2010
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
                
                % average spike amplitude
                stats{ipart}.(markername){itemp}.amplitude(itrial)                 = nanmean(amps);
                if isnan(nanmean(amps))
                    stats{ipart}.(markername){itemp}.burst_trialsum(itrial) = nan;
                end
                
            end % itrial
            
            ft_progress('close');
            
        end % itemp
        
        %%%%%%%%%%%
        %% SPIKY %%
        %%%%%%%%%%%
        
        if ~exist('SPIKY_check_spikes', 'file')
            disp('SPIKY is not in your path. Skipping SPIKY calculations');
            continue
        end
        
        % preallocate spike structure, with shuffled control
        
        clear spikedata        
        spikedata{size(SpikeTrials{ipart}.(markername).trialinfo, 1)}{size(SpikeTrials{ipart}.(markername).label, 2)} = [];
        spikedata_ctrl = spikedata;
        
        for iunit = 1:size(SpikeTrials{ipart}.(markername).label, 2)
            shuffled = SpikeTrials{ipart}.(markername).trial{iunit}(randperm(size(SpikeTrials{ipart}.(markername).trial{iunit}, 2)));
            for itrial = 1:size(SpikeTrials{ipart}.(markername).trialinfo, 1)
                spikedata{itrial}{iunit}        = SpikeTrials{ipart}.(markername).time{iunit}(SpikeTrials{ipart}.(markername).trial{iunit} == itrial);
                spikedata_ctrl{itrial}{iunit}   = SpikeTrials{ipart}.(markername).time{iunit}(shuffled == itrial);
            end
        end
        
        ft_progress('init','text', sprintf('Compute spike synchrony with Spiky \n%s p%d : %s', cfg.prefix(1:end-1), ipart, markername));
        for itrial = 1 : length(spikedata)
            
            ft_progress(0, 'processing trial %d from %d', itrial, length(spikedata));
            
            try
                
                % For trial data
                if length(spikedata{itrial}) < 2 
                    %need at least 2 spike trains to compute synchrony
                    stats{ipart}.(markername){itemp}.dist(:, itrial)            = nan;
                    stats{ipart}.(markername){itemp}.dist_label                 = nan;
                    stats{ipart}.(markername){itemp}.dist_perm(:, itrial)       = nan;
                    stats{ipart}.(markername){itemp}.dist_spikenr(:, itrial)    = nan;
                    stats{ipart}.(markername){itemp}.dist_isi(:, itrial)        = nan;
                    stats{ipart}.(markername){itemp}.dist_pooled(:, itrial)     = nan;
                    stats{ipart}.(markername){itemp}.dist_psth(:, itrial)       = nan;
                    continue
                end
                spikes                  = spikedata{itrial};
                ori_spikes              = spikes; %used to create control spikes data
                para.tmin               = SpikeTrials{ipart}.(markername).trialtime(itrial,1);
                para.tmax               = SpikeTrials{ipart}.(markername).trialtime(itrial,2);
                para.dts                = 1/SpikeTrials{ipart}.(markername).hdr.Fs;
                para.select_measures    = [0 1 0 0 0 0 0 0];  % {'ISI';'SPIKE';'RI_SPIKE';'SPIKE_realtime';'SPIKE_forward';'SPIKE_synchro';'SPIKE_order';'PSTH'};
                para.num_trains         = length(spikes);
                d_para                  = para;
                SPIKY_check_spikes
                para                    = d_para;
                spiky_ori               = SPIKY_loop_f_distances(spikes, para);
                
                % for control data with permuted trials
                % use same parameters
                spikes       = spikedata_ctrl{itrial};
                d_para       = para;
                SPIKY_check_spikes
                para         = d_para;
                ctrl_trials  = SPIKY_loop_f_distances(spikes, para);
                
                % for control data with random distribution, same spike numbers
                para.choice  = 1;
                spikes       = SPIKY_f_spike_train_surrogates(ori_spikes, para);
                ctrl_spikenr = SPIKY_loop_f_distances(spikes, para);
                
                % for control data with same ISI distribution
                para.choice = 2;
                spikes      = SPIKY_f_spike_train_surrogates(ori_spikes, para);
                ctrl_isi = SPIKY_loop_f_distances(spikes, para);
                
                % for control data with same pool of spikes
                para.choice = 3;
                spikes      = SPIKY_f_spike_train_surrogates(ori_spikes, para);
                ctrl_pooled = SPIKY_loop_f_distances(spikes, para);
                
                % for control data with same psth
                para.choice = 4;
                spikes      = SPIKY_f_spike_train_surrogates(ori_spikes, para);
                ctrl_psth   = SPIKY_loop_f_distances(spikes, para);
            catch
                warning('Something went wrong calculating Spiky values\n');
            end
            
            try
                % reshape to unit-by-unit, and retain only spike distance
                for itemp = 1 : size(SpikeTrials{ipart}.(markername).label, 2)
                    sel         = 1 : size(SpikeTrials{ipart}.(markername).label, 2);
                    sel(itemp)  = [];
                    stats{ipart}.(markername){itemp}.dist(:, itrial)            = spiky_ori.SPIKE.matrix(itemp, sel);
                    stats{ipart}.(markername){itemp}.dist_label                 = SpikeTrials{ipart}.(markername).label(sel);
                    stats{ipart}.(markername){itemp}.dist_perm(:, itrial)       = ctrl_trials.SPIKE.matrix(itemp, sel);
                    stats{ipart}.(markername){itemp}.dist_spikenr(:, itrial)    = ctrl_spikenr.SPIKE.matrix(itemp, sel);
                    stats{ipart}.(markername){itemp}.dist_isi(:, itrial)        = ctrl_isi.SPIKE.matrix(itemp, sel);
                    stats{ipart}.(markername){itemp}.dist_pooled(:, itrial)     = ctrl_pooled.SPIKE.matrix(itemp, sel);
                    stats{ipart}.(markername){itemp}.dist_psth(:, itrial)       = ctrl_psth.SPIKE.matrix(itemp, sel);
                end
            catch
                warning('Something went wrong calculating Spiky values\n');
            end
            
            clear spikes ori_spikes
        end
        ft_progress('close');
        
    end % markername
end % ipart

save(fname, 'stats', '-v7.3');
