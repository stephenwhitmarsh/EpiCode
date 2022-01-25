function [TFR] = dtx_TFRtrials(cfg, Trialdata, force)

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

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'TFRtrials.mat']);

if exist(fname_out,'file') && force == false
    fprintf('Loading precomputed TFR data\n');
    load(fname_out,'TFR');
    return
end
fprintf('(re-) computing TFR data\n');

TFR = [];
for ipart = 1 : size(Trialdata,2)
    for markername = string(fields(Trialdata{ipart}))'
        fprintf('For markername %s\n', markername);
        switch markername
            case "Interictal"
                cfgtemp            = [];
                cfgtemp.channel    = 'all';
                cfgtemp.method     = 'mtmconvol';
                cfgtemp.output     = 'pow';
                cfgtemp.taper      = 'dpss';
                cfgtemp.tapsmofrq  = 2; %number, the amount of spectral smoothing through multi-tapering. Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
                cfgtemp.pad        = 'nextpow2';
                cfgtemp.keeptrials = 'yes';
                cfgtemp.foi        = 5:1:100;
                cfgtemp.t_ftimwin  = ones(size(cfgtemp.foi)); %setparam
                cfgtemp.toi        = cfg.epoch.toi.(markername)(1) : 2 : cfg.epoch.toi.(markername)(2);
                TFR{ipart}.(markername)   = ft_freqanalysis(cfgtemp,Trialdata{ipart}.(markername));
            otherwise
                cfgtemp                         = [];
                cfgtemp.channel                 = 'all';
                cfgtemp.method                  = 'mtmconvol';
                cfgtemp.output                  = 'pow';
                cfgtemp.taper                   = 'hanning';
                cfgtemp.pad                     = 'nextpow2';
                cfgtemp.keeptrials              = 'yes';
                cfgtemp.foi                     = 5:1:100;
                cfgtemp.t_ftimwin               = 7./cfgtemp.foi;
                cfgtemp.toi                     = cfg.epoch.toi.(markername)(1) : 0.012 : cfg.epoch.toi.(markername)(2);
                TFR{ipart}.(markername)         = ft_freqanalysis(cfgtemp,Trialdata{ipart}.(markername));
        end
    end
end % ipart
save(fname_out,'TFR','-v7.3');