function MuseStruct = merge_MuseMarkers(cfg, MuseStruct)

% Use as : 
%           MuseStruct = merge_MuseMarkers(cfg, MuseStruct)
% 
% merge_MuseMarkers merges several different muse markers into a new one
% 
% Input :
% cfg.markers_to_merge = {'marker1', 'marker2', 'marker3'};
% cfg.newmarker        = 'new_name';
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


for ipart = 1:size(MuseStruct, 2)
    for idir = 1:size(MuseStruct{ipart}, 2)
        synctime = [];
        clock = [];
        for markername = string(cfg.markers_to_merge)
            if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'synctime')
                continue
            end
            synctime = [synctime, MuseStruct{ipart}{idir}.markers.(markername).synctime];
            clock    = [clock, MuseStruct{ipart}{idir}.markers.(markername).clock];
        end
        synctime = sort(synctime);
        clock    = sort(clock);
        MuseStruct{ipart}{idir}.markers.(cfg.newmarker).synctime = synctime;
        MuseStruct{ipart}{idir}.markers.(cfg.newmarker).clock = clock;
    end
end