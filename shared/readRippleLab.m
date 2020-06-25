function [dat, tbl] = readRippleLab(fname)

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

if ~exist(fname,'file')
    fprintf('%s does not exist!\n',fname);
    return
end

% Load data
dat_rl      = load(fname,'-mat'); %

% Extract header
f           = fields(dat_rl);
mont        = f{find(~strcmp(f,'st_FileData'))};
temp        = dat_rl.(f{find(~strcmp(f,'st_FileData'))});
trials      = temp.v_Intervals;
Fs          = dat_rl.st_FileData.v_SampleRate;
label       = dat_rl.st_FileData.v_Labels;

evenements      = temp.st_HFOInfo.m_EvtLims;                            % absolute limit of event (number of sample)
evenements_rel  = temp.st_HFOInfo.m_Rel2IntLims;                        % relative limit of event (number of sample)
debut           = temp.s_StartSample;                                   %
Fs              = temp.st_HFOInfo.s_Sampling;                           % sampling rate

% Create fieldtrip structure
dat                     = [];
dat.label               = label;
dat.trialinfo           = temp.st_HFOInfo.v_EvType; % name of the event (Ripples,FastRipple,..)
dat.trialinfo(:,2:3)    = evenements_rel;
dat.sampleFs            = Fs;

for itrial = 1 : size(trials,1)
    dat.trial{itrial}   = trials{itrial}';
    dat.time{itrial}    = linspace(-length(trials{itrial}) / Fs / 2, length(trials{itrial}) / Fs / 2, length(trials{itrial}) );
end

% fout = fullfile(analysisdatadir,fout_data{ifile});
% fprintf('Saving FieldTrip data to: %s\n',fout)
% save(fout,'dat')

% create table with results
tbl             = table;
tbl.Eventnr     = (1:length(evenements(:,1)))';
tbl.Eventtype   = temp.st_HFOInfo.v_EvType;
tbl.EventStart  = ((debut/Fs)+(evenements(:,1)/Fs));
tbl.Duration    = (evenements(:,2)-evenements(:,1))/(Fs/1000);    % duratin of the event
