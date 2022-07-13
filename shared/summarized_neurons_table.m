function neurons_table = summarized_neurons_table(cfg, markername, force, spikestats, waveformstats, spiketrials)

% Use as : 
%       neurons_table = summarized_neurons_table(cfg, ipart, markername, force, spikestats, waveformstats, spiketrials)

% summarized_neurons_table creates a table with the summary of the main
% descriptive measurements computed for each unit.
% 
% Input : 
% - markername    : name of the analysis to use to create the table
% - force         : whether to redo analyses or read previous save (true/false)
% - spikestats    : output from spikeTrialStats.m
% - waveformstats : output from spikeWaveformStats.m
% - spiketrials   : output from readSpikeTrials.m
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
% 

for ipart = 1 : size(cfg.directorylist, 2)
    
    fname = fullfile(cfg.tablesavedir,[cfg.prefix, 'p', num2str(ipart), '-neurons_table_', char(markername), '.xlsx']);
    
    if exist(fname, 'file') && force == false
        fprintf('Reading %s\n', fname);
        neurons_table{ipart} = readtable(fname);
        continue
    end
    
    neurons_table{ipart} = table.empty;
    
    for itemp = 1:size(spikestats{ipart}.(markername), 2)
        
        %check that the different structure have the same order :
        assert(strcmp(waveformstats{ipart}.label{itemp}, spikestats{ipart}.(markername){itemp}.label));
        assert(strcmp(waveformstats{ipart}.label{itemp}, spiketrials{ipart}.(markername).label{itemp}));
        
        neurons_table{ipart}.patient_ID{itemp}                 = cfg.prefix(1:end-1);
        neurons_table{ipart}.clusterID{itemp}                  = spikestats{ipart}.(markername){itemp}.label;
        neurons_table{ipart}.cluster_group{itemp}              = strrep(waveformstats{ipart}.cluster_group{itemp}, ' ', '');
        maxchan       = spiketrials{ipart}.(markername).template_maxchan(itemp) + 1; %+1 because starts at zero
        neurons_table{ipart}.channel{itemp}                    = cfg.circus.channel{maxchan};
        if ~isfield(cfg.circus, 'channelname')
            a = split(neurons_table{ipart}.clusterID{itemp}, '_');
            neurons_table{ipart}.electrode_bundle{itemp}       = a{3};
        end
        neurons_table{ipart}.RPV{itemp}                        = spikestats{ipart}.(markername){itemp}.RPV * 100;
        neurons_table{ipart}.n_trials{itemp}                   = size(spikestats{ipart}.(markername){itemp}.trialinfo, 1);
        neurons_table{ipart}.n_PAs{itemp}                      = length(spikestats{ipart}.(markername){itemp}.isi);
        
        neurons_table{ipart}.freq_Hz_mean{itemp}               = nanmean(spikestats{ipart}.(markername){itemp}.trialfreq);
        neurons_table{ipart}.freq_Hz_std{itemp}                = nanstd(spikestats{ipart}.(markername){itemp}.trialfreq);
        neurons_table{ipart}.amplitude_uV_mean{itemp}          = nanmean(spikestats{ipart}.(markername){itemp}.amplitude);
        neurons_table{ipart}.amplitude_uV_std{itemp}           = nanstd(spikestats{ipart}.(markername){itemp}.amplitude);
        neurons_table{ipart}.CV2_mean{itemp}                   = nanmean(spikestats{ipart}.(markername){itemp}.CV2_trial);
        neurons_table{ipart}.CV2_std{itemp}                    = nanstd(spikestats{ipart}.(markername){itemp}.CV2_trial);
        neurons_table{ipart}.burst_per_min_mean{itemp}         = nanmean(spikestats{ipart}.(markername){itemp}.burst_trialsum);
        neurons_table{ipart}.burst_per_min_std{itemp}          = nanstd(spikestats{ipart}.(markername){itemp}.burst_trialsum);
        
        %spike distance (synchrony)
        dist = [];
        if ~isfield(cfg.circus, 'channelname')
            dist = mean(spikestats{ipart}.(markername){itemp}.dist, 2);
        else
            % use only  units on the same bundle of electrodes :
            for iunit = 1:size(spikestats{ipart}.(markername){itemp}.dist_label, 2)
                unitname = spikestats{ipart}.(markername){itemp}.dist_label{iunit};
                a = split(unitname, '_');
                b = split(neurons_table{ipart}.clusterID{itemp}, '_');
                if strcmp(a{3}, b{3})
                    dist = [dist, mean(spikestats{ipart}.(markername){itemp}.dist(iunit, :), 2)];
                end
            end
        end
        neurons_table{ipart}.spike_distance_mean{itemp} = mean(dist);
        neurons_table{ipart}.spike_distance_std{itemp}  = std(dist);
        
        neurons_table{ipart}.halfwidth_ms{itemp}             = waveformstats{ipart}.halfwidth.val(itemp) * 1000;
        neurons_table{ipart}.peaktrough_ms{itemp}            = waveformstats{ipart}.peaktrough.val(itemp) * 1000;
        neurons_table{ipart}.troughpeak_ms{itemp}            = waveformstats{ipart}.troughpeak.val(itemp) * 1000;
    end
    
    delete(fname);
    writetable(neurons_table{ipart}, fname);
end
