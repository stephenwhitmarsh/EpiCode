function data = remove_artefacted_trials(cfg, data)

% Use as
%       data = remove_artefacted_trials(cfg, data)
% 
% Remove artefacted trials in raw and spike data, based on trialinfo.BAD_sec.
% 
% Input : 
% - cfg.minbadtime.(markername) : minimum duration, in second, of artefacts in
%                                 trialinfo.BAD_sec to consider as an artefact 
%                                 and remove it from the data structure (default = 0)
% - data                        : raw data or spike data, in the "EpiCode" format
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


if isempty(data)
    return
end

cfg.minbadtime = ft_getopt(cfg, 'minbadtime', []);

for ipart = 1:size(data,2)
    
    fprintf('For part %d\n', ipart);
    
    for markername = string(fieldnames(data{ipart}))'
        
        fprintf('*** For %s ***\n', markername);
        
        cfg.minbadtime.(markername) = ft_getopt(cfg.minbadtime, char(markername), 0);
        
        if isempty(data{ipart}.(markername))
            continue
        end
        
        type = ft_datatype(data{ipart}.(markername));
       
        sel = data{ipart}.(markername).trialinfo.BAD_sec <= cfg.minbadtime.(markername);
        cfgtemp        = [];
        cfgtemp.trials = sel;
        
        if strcmp(type, 'raw')
            data{ipart}.(markername) = ft_selectdata(cfgtemp, data{ipart}.(markername));
            
        elseif strcmp(type, 'spike')
            data{ipart}.(markername) = ft_spike_select_rmfulltrials(cfgtemp, data{ipart}.(markername));
            
        else
            error('datatype %s is not supported', type);
        end
       
    end %markername
end %ipart

end