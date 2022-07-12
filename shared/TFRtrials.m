function [TFR] = TFRtrials(cfg, force, Trialdata)

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


cfg.TFR.postfix     = ft_getopt(cfg.TFR, 'postfix', []);
write               = ft_getopt(cfg.TFR, 'write', true);

if nargin == 1
    for markername = string(cfg.TFR.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'TFR_', markername, cfg.TFR.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
            for ipart = 1 : size(cfg.directorylist, 2)
                try
                    TFR{ipart}.(markername) = temp.TFR{ipart}.(markername);
                catch
                end
            end
        else
            fprintf('(re-) computing TFR data for %s\n', markername);
        end
    end
    return
    
elseif ~force
    missing = [];
    for markername = string(cfg.TFR.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'TFR_', markername, cfg.TFR.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
            for ipart = 1 : size(cfg.directorylist, 2)
                try
                    TFR{ipart}.(markername) = temp.TFR{ipart}.(markername);
                catch
                end
            end
        else
            fprintf('(re-) computing TFR data for %s\n', markername);
            missing = [missing, markername];
        end
    end
    cfg.TFR.name = missing;
end

if nargin < 3 && ~isempty(cfg.TFR.name) 
    Trialdata = readLFP(cfg);
end

% loop over markers
for markername = string(cfg.TFR.name)
    
    fname_out = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'TFR_', markername, '.mat'));
    
    if exist(fname_out, 'file') && force == false
        fprintf('Loading precomputed TFR data for %s\n', markername);
        
        temp = load(fname_out, 'TFR');
        for ipart = 1 : size(MuseStruct, 2)
            TFR{ipart}.(markername) = temp.TFR{ipart}.(markername);
        end
        continue
        
    else
        fprintf('(re-) computing TFR data for %s\n', markername);
    end
    
    % loop over parts within subject
    for ipart = 1 : size(cfg.directorylist, 2)
        
        if ~isempty(Trialdata{ipart}.(markername))
            cfgtemp                         = [];
            cfgtemp.channel                 = 'all';
            cfgtemp.method                  = 'mtmconvol';
            cfgtemp.output                  = 'pow';
            cfgtemp.taper                   = 'hanning';
%             cfgtemp.taper                   = 'dpss';
%             cfgtemp.tapsmofrq               = 5;
            cfgtemp.pad                     = 'nextpow2';
            cfgtemp.keeptrials              = cfg.TFR.keeptrials;
            cfgtemp.foi                     = cfg.TFR.foi.(markername);
            cfgtemp.t_ftimwin               = cfg.TFR.t_ftimwin.(markername);
            cfgtemp.toi                     = cfg.TFR.toi.(markername);
            cfgtemp.feedback                = 'off';
            TFR{ipart}.(markername)         = ft_freqanalysis(cfgtemp, Trialdata{ipart}.(markername));
            
            % to save memory, remove if not already done
            if isfield(TFR{ipart},'cfg')
                TFR{ipart} = rmfield(TFR{ipart}, 'cfg');
            end
        else
            sprintf('No LFP for %s part %d', markername, ipart);
        end
    end % markername
    if write
        fprintf('Saving TFR data for %s\n', markername);
        saveMarker_TFR(TFR, markername, fname_out)
    end
end % ipart

