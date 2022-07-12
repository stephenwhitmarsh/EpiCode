function data = remove_too_short_trials(cfg, data)
 
% Use as : 
%       data = remove_too_short_trials(cfg, data)

% Remove all trials which are too short.
% 
% Input : 
% - cfg.min_trial_length.(markername) : duration, in second, below which trials 
%                                       are considered too short and are removed (default = 0)
% - data                              : raw data or spike data, in the "EpiCode" format
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

cfg.min_trial_length = ft_getopt(cfg, 'min_trial_length', []);

for ipart = 1 :size(data, 2)    
    for markername = string(fieldnames(data{ipart}))'
        
        if isfield(data{ipart}.(markername), 'fsample')
            fs = data{ipart}.(markername).fsample;
        elseif isfield(data{ipart}.(markername), 'hdr')
            fs = data{ipart}.(markername).hdr.Fs;
        else
            error('cannot find sampling frequency in the data structure');
        end
        
        cfg.min_trial_length.(markername) = ft_getopt(cfg.min_trial_length, char(markername), 0);
        trials_length = (data{ipart}.(markername).trialinfo.endsample - data{ipart}.(markername).trialinfo.begsample) / fs;
        
        datatype = ft_datatype(data{ipart}.(markername));
        
        cfgtemp = [];
        cfgtemp.trials = trials_length >= cfg.min_trial_length.(markername);
        
        if strcmp(datatype, 'raw')
            data{ipart}.(markername) = ft_selectdata(cfgtemp, data{ipart}.(markername));
        elseif strcmp(datatype, 'spike')
            data{ipart}.(markername) = ft_spike_select_rmfulltrials(cfgtemp, data{ipart}.(markername));
        else
            error('cannot find data type');
        end
        
    end
end
