function [FFT] = FFTtrials(cfg, force, LFP)
% 
% Use as : 
%       [FFT] = FFTtrials(cfg, force, LFP)
% 
% The LFP argument is optional. If LFP is not in the input, then it is 
% loaded with the function readLFP.m.
% 
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

cfg.FFT             = ft_getopt(cfg, 'FFT', []);
cfg.FFT.postfix     = ft_getopt(cfg.FFT, 'postfix', []);
cfg.FFT.keeptrials  = ft_getopt(cfg.FFT, 'keeptrials', 'yes');
write               = ft_getopt(cfg.FFT, 'write', true);

if nargin == 1
    for markername = string(cfg.FFT.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'FFT_', markername, cfg.FFT.postfix, '.mat'));
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
                    FFT{ipart}.(markername) = temp.FFT{ipart}.(markername);
                catch
                end
            end
        else
            fprintf('(re-) computing FFT data for %s\n', markername);
        end
    end
    return
    
elseif ~force
    missing = [];
    for markername = string(cfg.FFT.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'FFT_', markername, cfg.FFT.postfix, '.mat'));
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
                    FFT{ipart}.(markername) = temp.FFT{ipart}.(markername);
                catch
                end
            end
        else
            fprintf('(re-) computing FFT data for %s\n', markername);
            missing = [missing, markername];
        end
    end
    cfg.FFT.name = missing;
end

if nargin < 3
    load_LFP = true;
end

if ~isempty(cfg.FFT.name)
    
    if load_LFP
        cfg.LFP.name = cfg.FFT.name;
        LFP = readLFP(cfg);
    end
    
    % loop over markers
    for markername = string(cfg.FFT.name)
        
        fname_out = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'FFT_', markername, cfg.FFT.postfix, '.mat'));

        % loop over parts within subject
        for ipart = 1 : size(cfg.directorylist, 2)
            
            fprintf('Processing "%s", part %d of %d\n', markername, ipart, size(cfg.directorylist, 2))

            if isempty(LFP{ipart}.(markername))
                fprintf('No LFP for %s part %d', markername, ipart);
                continue
            end
                                    
            cfgtemp                         = [];
            cfgtemp.channel                 = 'all';
            cfgtemp.method                  = 'mtmfft';
            cfgtemp.output                  = 'pow';
            cfgtemp.taper                   = 'hanning';
            cfgtemp.keeptrials              = cfg.FFT.keeptrials;
            cfgtemp.foi                     = cfg.FFT.foi.(markername);
            cfgtemp.pad                     = 'nextpow2';
            FFT{ipart}.(markername)         = ft_freqanalysis(cfgtemp, LFP{ipart}.(markername));
            
            % to save memory, remove cfg if not already done
            if isfield(FFT{ipart},'cfg')
                FFT{ipart} = rmfield(FFT{ipart}, 'cfg');
            end
                      
        end % ipart
    end % markername
    
    if write
        fprintf('Saving FFT data for %s\n', markername);
        saveMarker_FFT(FFT, markername, fname_out)
    end
end