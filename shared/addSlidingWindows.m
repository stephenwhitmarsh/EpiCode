function [config, MuseStruct] = addSlidingWindows(config, MuseStruct)

% [config, MuseStruct] = addSlidingWindows(config, MuseStruct)
% Add a "window" field to the config and MuseStruct, allowing sliding
% time-window analyses for both LFP, TFR, FFT and spike data.
%
% ### Necessary input:
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% MuseStruct            = info (e.g. events, files) of original data,
%                         used to segment the spikes into trials
%
% config.window.length  = 60 (seconds)
% config.window.overlap = 0.5 (fraction)
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

% get the default cfg options
config.window                  = ft_getopt(config,        'window', []);
config.window.length           = ft_getopt(config.window, 'length', 60);
config.window.overlap          = ft_getopt(config.window, 'overlap', 0);

% window is symmetrical around 0
length = config.window.length;

for ipart = 1 : size(config.directorylist, 2)
    for idir = 1 : size(config.directorylist{ipart}, 2)
        sprintf('Adding part %d, dir %d\n', ipart, idir)
        temp = dir(fullfile(config.rawdir, config.directorylist{ipart}{idir}, ['*', config.circus.channel{1}, '.ncs']));
        hdr  = ft_read_header(fullfile(config.rawdir, config.directorylist{ipart}{idir}, temp.name));
        MuseStruct{ipart}{idir}.markers.window__START__.synctime = 0 : (length - length * config.window.overlap) : (hdr.nSamples / hdr.Fs - length / 2);
        MuseStruct{ipart}{idir}.markers.window__START__.clock    = seconds(MuseStruct{ipart}{idir}.markers.window__START__.synctime) + MuseStruct{ipart}{idir}.starttime;
        MuseStruct{ipart}{idir}.markers.window__START__.events   = size(MuseStruct{ipart}{idir}.markers.window__START__.synctime, 2);
        MuseStruct{ipart}{idir}.markers.window__END__.synctime   = MuseStruct{ipart}{idir}.markers.window__START__.synctime + length;
        MuseStruct{ipart}{idir}.markers.window__END__.clock      = MuseStruct{ipart}{idir}.markers.window__START__.clock + seconds(length);
        MuseStruct{ipart}{idir}.markers.window__END__.events     = size(MuseStruct{ipart}{idir}.markers.window__END__.synctime, 2);
    end
end

config.muse.startmarker.window  = 'window__START__';
config.muse.endmarker.window    = 'window__END__';

config.epoch.toi.window         = [0 config.window.length];
config.epoch.pad.window         = 0;

% add window configuration to settings 

if isfield(config, 'LFP')
    if isfield(config.LFP, 'name')
        config.LFP.toi.window           = [0 0];
        config.LFP.pad.window           = 0;
        config.epoch.toi.window         = [0 0];
        config.epoch.pad.window         = 0;        
    end
end

if isfield(config, 'FFT')
    if isfield(config.FFT, 'name')
        config.FFT.toi.window           = [0 0];
        config.FFT.pad.window           = 0;
    end
end

if isfield(config, 'TFR')
    if isfield(config.TFR, 'name')
        config.TFR.toi.window           = [0 0];
        config.TFR.pad.window           = 0;
    end
end

if isfield(config, 'spike')
    if isfield(config.spike, 'name')
        config.spike.toi.window         = [0 0];
        config.spike.pad.window         = 0;
    end
end


