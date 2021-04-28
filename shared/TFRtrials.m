function [TFR] = TFRtrials(cfg, LFP, force)

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

fname_out = fullfile(cfg.datasavedir,[cfg.prefix, 'TFRtrials.mat']);

if nargin == 1
    if exist(fname_out, 'file')
        fprintf('Reading %s\n', fname_out);
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname_out, 'TFR');

            catch ME
                err_count = err_count + 1;
            end
            count = count + 1;
        end
        return   
    else
        fprintf('Computing TFR\n');
    end
end

if exist(fname_out, 'file') && force == false
    fprintf('Loading TFR\n');
    load(fname_out, 'TFR');
    return
end

fprintf('Computing TFR\n');

for ipart = 1 : size(LFP, 2)
    
    for markername = string(fields(LFP{ipart}))'
        
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all'; 
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'hanning';
        cfgtemp.pad                     = 'nextpow2'; 
        cfgtemp.keeptrials              = cfg.TFR.keeptrials;
        cfgtemp.foi                     = cfg.TFR.foi.(markername);
        cfgtemp.t_ftimwin               = cfg.TFR.t_ftimwin.(markername);
        cfgtemp.toi                     = cfg.TFR.toi.(markername);
        cfgtemp.feedback                = 'off';
        TFR{ipart}.(markername)         = ft_freqanalysis(cfgtemp, LFP{ipart}.(markername));
        
        % to save memory, remove if not already done
        try
            TFR{ipart}                  = rmfield(TFR{ipart}, 'cfg');
        catch
        end
        
    end % markername
end % ipart
    
save(fname_out, 'TFR', '-v7.3');

