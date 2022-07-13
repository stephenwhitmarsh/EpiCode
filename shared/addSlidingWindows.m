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
config.window.name             = ft_getopt(config.window, 'name', {'window'});
config.window.length           = ft_getopt(config.window, 'length', []);
config.window.overlap          = ft_getopt(config.window, 'overlap', []);

for markername = string(config.window.name)
    
    % window is symmetrical around 0
    length = config.window.length.(markername);
    start_name = sprintf('%s__START__', markername);
    end_name   = sprintf('%s__END__', markername);
    for ipart = 1 : size(config.directorylist, 2)
        for idir = 1 : size(config.directorylist{ipart}, 2)
            
            ft_info('Adding %s in part %d, dir %d\n', markername, ipart, idir)
            
            temp = dir(fullfile(config.rawdir, config.directorylist{ipart}{idir}, ['*', config.circus.channel{1}, '.ncs']));
            hdr  = ft_read_header(fullfile(config.rawdir, config.directorylist{ipart}{idir}, temp.name));
            
            MuseStruct{ipart}{idir}.markers.(start_name).synctime = 0 : (length - length * config.window.overlap.(markername)) : (hdr.nSamples / hdr.Fs - length / 2);
            MuseStruct{ipart}{idir}.markers.(start_name).clock    = seconds(MuseStruct{ipart}{idir}.markers.(start_name).synctime) + MuseStruct{ipart}{idir}.starttime;
            MuseStruct{ipart}{idir}.markers.(start_name).events   = size(MuseStruct{ipart}{idir}.markers.(start_name).synctime, 2);
            MuseStruct{ipart}{idir}.markers.(end_name).synctime   = MuseStruct{ipart}{idir}.markers.(start_name).synctime + length;
            MuseStruct{ipart}{idir}.markers.(end_name).clock      = MuseStruct{ipart}{idir}.markers.(start_name).clock + seconds(length);
            MuseStruct{ipart}{idir}.markers.(end_name).events     = size(MuseStruct{ipart}{idir}.markers.(end_name).synctime, 2);
            
        end
    end
    
    config.muse.startmarker.(markername)  = start_name;
    config.muse.endmarker.(markername)    = end_name;
    
    config.epoch.toi.(markername)         = [0 config.window.length.(markername)];
    config.epoch.pad.(markername)         = 0;
    
    % add window configuration to settings
    
    if isfield(config, 'LFP')
        if isfield(config.LFP, 'name')
            config.LFP.toi.(markername)           = [0 0];
            config.LFP.pad.(markername)           = 0;
            config.epoch.toi.(markername)         = [0 0];
            config.epoch.pad.(markername)         = 0;
        end
    end
    
    if isfield(config, 'FFT')
        if isfield(config.FFT, 'name')
            config.FFT.toi.(markername)           = [0 0];
            config.FFT.pad.(markername)           = 0;
        end
    end
    
    if isfield(config, 'TFR')
        if isfield(config.TFR, 'name')
            config.TFR.toi.(markername)           = [0 0];
            config.TFR.pad.(markername)           = 0;
        end
    end
    
    if isfield(config, 'spike')
        if isfield(config.spike, 'name')
            config.spike.toi.(markername)         = [0 0];
            config.spike.pad.(markername)         = 0;
        end
    end
end

if isempty(markername)
    warning('No sliding window was added');
end
