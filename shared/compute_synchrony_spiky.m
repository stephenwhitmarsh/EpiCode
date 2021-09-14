function [spiky_results] = compute_synchrony_spiky(cfg, SpikeTrials, force, suffix)

% COMPUTE_SYNCHRONY_SPIKY computes spike synchrony with the SPIKY toolbox,
% on spike data created with readSpikeTrials_MuseMarkers.m or
% readSpikeTrials_windowed, and outputs the results.
% 
% Use as
%   [spiky_results] = compute_synchrony_spiky(cfg, SpikeTrials, force, suffix)
%
% SPIKY toolbox is available here : 
% http://www.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/SPIKY.html
%
% the 'SPIKY' folder must be added to path, including 'SPIKY/SPIKE_MEX' 
%
% The 'suffix' argument is optional
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
    suffix = string.empty;
end

fname = fullfile(cfg.datasavedir, sprintf('%sspiky_results%s.mat',cfg.prefix, suffix));
if exist(fname, 'file') && force == false
    fprintf('Reading %s \n', fname);
    load(fname, 'spiky_results');
    return
end

for ipart = 1:size(SpikeTrials, 2)
    for markername = string(fieldnames(SpikeTrials{ipart}))'

        %create spike structure to use with spiky
        clear spikedata
        for iunit = 1:size(SpikeTrials{ipart}.(markername).label, 2)
            for itrial = 1:size(SpikeTrials{ipart}.(markername).trialinfo, 1)
                idx                      = SpikeTrials{ipart}.(markername).trial{iunit} == itrial;
                spikedata{itrial}{iunit} = SpikeTrials{ipart}.(markername).time{iunit}(SpikeTrials{ipart}.(markername).trial{iunit} == itrial);
            end
        end
        
        %control : permute units from different trials
        %each unit come from a different trial. There are not several units of
        %the same trial. Except if there are less trials than units
        % TODO: use PERMUTE function
        clear spikedata_ctrl
        for itrial = 1:length(spikedata)
            trial_list = 1:length(spikedata);
            trial_list(trial_list==itrial) = []; % remove the current trial of the random selection of trials
            for i_unit = 1:length(spikedata{itrial})
                if isempty(trial_list) %if more trials than units, start again the random selection 
                    trial_list = 1:length(spikedata);
                    trial_list(trial_list==itrial) = [];
                end
                rand_trial = trial_list(randi(length(trial_list)));
                trial_list(trial_list==rand_trial) = []; %remove the selected trial for the next unit
                spikedata_ctrl{itrial}{i_unit} = spikedata{rand_trial}{i_unit};
            end
        end
               
        trialtime   = SpikeTrials{ipart}.(markername).trialtime;
        
        %% spiky loop
        ft_progress('init','text', sprintf('Compute spike synchrony with Spiky \n%s p%d : %s', cfg.prefix(1:end-1), ipart, markername));
        for itrial = 1:length(spikedata)
            
            ft_progress(0, 'processing trial %d from %d', itrial, length(spikedata));
            
            % For trial data
            spikes      = spikedata{itrial};
            ori_spikes  = spikes; %used to create control spikes data
            para.tmin   = trialtime(itrial,1);
            para.tmax   = trialtime(itrial,2);
            para.dts    = 1/SpikeTrials{ipart}.(markername).hdr.Fs;
            para.select_measures        = [1 1 0 0 0 0 0 1];  % Select measures (0-do not calculate,1-calculate)
            
            %TODO: add in config of function
            m_para.all_measures_string  = {'ISI';'SPIKE';'RI_SPIKE';'SPIKE_realtime';'SPIKE_forward';'SPIKE_synchro';'SPIKE_order';'PSTH'};  % order of select_measures
            para.num_trains             = length(spikes);
            
            %SPIKY setting of params
            d_para = para;
            SPIKY_check_spikes
            para = d_para;
            % TODO: add in config of function  
            spiky_results{ipart}.(markername).trials{itrial} = SPIKY_loop_f_distances(spikes,para); 
            
            % for control data with trials merged
            %same parameters
            spikes = spikedata_ctrl{itrial};
            d_para = para;
            SPIKY_check_spikes
            para=d_para;
            spiky_results{ipart}.(markername).ctrl_trials{itrial} = SPIKY_loop_f_distances(spikes,para);
            
            % for control data with random distribution, same spike numbers
            para.choice = 1;%1;%spike number;%4; %select psth conservation
            spikes      = SPIKY_f_spike_train_surrogates(ori_spikes,para);
            spiky_results{ipart}.(markername).ctrl_spikenr{itrial} = SPIKY_loop_f_distances(spikes,para);

            % for control data with same ISI distribution
            para.choice = 2;%1;%spike number;%4; %select psth conservation
            spikes      = SPIKY_f_spike_train_surrogates(ori_spikes,para);
            spiky_results{ipart}.(markername).ctrl_isi{itrial} = SPIKY_loop_f_distances(spikes,para);
            
            % for control data with same pool of spikes
            para.choice = 3;%1;%spike number;%4; %select psth conservation
            spikes      = SPIKY_f_spike_train_surrogates(ori_spikes,para);
            spiky_results{ipart}.(markername).ctrl_pooled{itrial} = SPIKY_loop_f_distances(spikes,para);
         
            % for control data with same psth
            para.choice = 4;%1;%spike number;%4; %select psth conservation
            spikes      = SPIKY_f_spike_train_surrogates(ori_spikes,para);
            spiky_results{ipart}.(markername).ctrl_psth{itrial}    = SPIKY_loop_f_distances(spikes,para);
            
            clear spikes ori_spikes
        end
        ft_progress('close');
        
%         %gather results to average
%         clear temp
%         for ifield = ["trials","ctrl_trials","ctrl_psth", "ctrl_pooled", "ctrl_isi", "ctrl_spikenr"]
%             for idata = ["ISI", "SPIKE", "PSTH"]
%                 for itrial = 1:length(spiky_results{ipart}.(markername).ctrl_psth)
%                     temp.(ifield).(idata).time{itrial}  = spiky_results{ipart}.(markername).(ifield){itrial}.(idata).time;% - results{itrial}.SPIKE.time(end); %timelock to the end
%                     temp.(ifield).(idata).trial{itrial} = spiky_results{ipart}.(markername).(ifield){itrial}.(idata).profile;
%                 end
%             end
%         end
%         
%         %average profiles and matrix
%         for ifield = ["trials", "ctrl_trials", "ctrl_psth", "ctrl_pooled", "ctrl_isi", "ctrl_spikenr"]
%             for idata = ["ISI", "SPIKE", "PSTH"]
%                 %profiles
%                 fprintf('average %s dissimilarity profiles an matrix for %s\n', idata, ifield);
%                 
%                 sel = ~cellfun(@(c) all(c==0), temp.(ifield).(idata).trial);
%                 
%                 [spiky_results{ipart}.(markername).profileavg.(idata).(ifield).time, spiky_results{ipart}.(markername).profileavg.(idata).(ifield).avg] = ...
%                     SPIKY_f_average_pi(temp.(ifield).(idata).time(sel), temp.(ifield).(idata).trial(sel), para.dts);
%                 %matrix
%                 if ~strcmp(idata, 'PSTH')
%                     for itrial = 1:length(spiky_results{ipart}.(markername).ctrl_psth)
%                         mat{itrial} = spiky_results{ipart}.(markername).(ifield){itrial}.(idata).matrix;
%                     end
%                     mat_concat = cat(3, mat{:});
%                     spiky_results{ipart}.(markername).matrixavg.(idata).(ifield) = mean(mat_concat, 3, 'omitnan');
%                     %remove 1 in the diagonal
%                     spiky_results{ipart}.(markername).matrixavg.(idata).(ifield)(spiky_results{ipart}.(markername).matrixavg.(idata).(ifield) == 1) = NaN;
%                 end
%             end
%         end
%         
        %keep infos for later analysis
        for ifield = ["trialinfo", "trialtime", "label", "hdr", "channelname"]
            spiky_results{ipart}.(markername).(ifield) = SpikeTrials{ipart}.(markername).(ifield);
        end
        spiky_results{ipart}.(markername) = orderfields(spiky_results{ipart}.(markername));
    end
end

fprintf('Save Spiky data to %s \n', fname);
save(fname,'spiky_results','-v7.3');