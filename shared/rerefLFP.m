function [LFP] = rerefLFP(cfg, MuseStruct, force)

% READLFP rereferences epoched data according to details in the config.
% Use as
%   [LFP] = rerefLFP(cfg, LFP, force)
%
% cfg.LFP.demean              = 'yes';
% cfg.LFP.baselinewindow      = [-0.5 0];
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

if nargin < 3
    force = false;
end

write               = ft_getopt(cfg.LFP, 'write', true);
cfg.LFP.postfix     = ft_getopt(cfg.LFP, 'postfix', []);

% if not, loop over markers
for markername = string(cfg.LFP.name)
    
    fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'LFPreref_', markername, cfg.LFP.postfix, '.mat'));
    
    if exist(fname, 'file') && ~force
        fprintf('Reading %s\n', fname);
        count = 0;
        err_count = 0;
        while count == err_count
            try
                temp = load(fname);
                for ipart = 1 : size(cfg.directorylist, 2)
                    LFP{ipart}.(markername) = temp.LFP{ipart}.(markername);
                end
            catch ME
                err_count = err_count + 1;
            end
            count = count + 1;
        end
    else
        fprintf('Will be (re-) computing LFP data for %s\n', markername);
        
        fname_orig = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'LFP_', markername, cfg.LFP.postfix, '.mat'));
        if exist(fname_orig, 'file')
            fprintf('Reading %s\n', fname_orig);
            
            temp = load(fname_orig);
            for ipart = 1 : size(MuseStruct, 2)
                LFP{ipart}.(markername) = temp.LFP{ipart}.(markername);
            end
        else
            fprintf('You still need to calculate the LFP!');
            for ipart = 1 : size(MuseStruct, 2)              
                LFP{ipart}.(markername) = [];
            end
        end
        
        if strcmp(cfg.template.reref, 'yes')
            
            for ipart = 1 : size(MuseStruct, 2)
                
                if isempty(LFP{ipart}.(markername))
                    continue
                end
                
                if ~contains(LFP{ipart}.(markername).label{1}, '-')
                    
                    labels_nonum    = regexprep(LFP{ipart}.(markername).label, '[0-9_]', '');
                    [~, ~, indx]    = unique(labels_nonum);
                    
                    % average per part (day) then reref
                    clear group
                    for i = 1 : max(indx)
                        cfgtemp             = [];
                        cfgtemp.reref       = 'yes';
                        cfgtemp.refmethod   = 'bipolar';
                        cfgtemp.demean      = 'yes';
                        cfgtemp.baselinewindow = [-0.3, -0.1];
                        cfgtemp.channel     = LFP{ipart}.(markername).label(indx==i);
                        group{i}            = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
                    end
                    LFP{ipart}.(markername) = ft_appenddata([], group{:});
                else
                    sprintf('Already rereferenced');
                end
            end
        end
        
        if write
            fprintf('Saving rereferenced LFP data for %s\n', markername);
            saveMarker_LFP(LFP, markername, fname)
        end
    end

end
