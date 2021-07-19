function [MuseStruct_new] = editMuseMarkers(cfg, MuseStruct_orig)

% EDITMUSEMARKERS adds, removes or renames markers in each MuseMarker structure
% 
% Usage: Use readMuseMarkers to obtain MuseStruct, and writeMuseMarkers to 
% write new MuseStruct to Muse event files.
% 
% Optional fields:
% 
% cfg.editmarkerfile.torename       = {'old1', 'new1'; 'old2', 'new2'; 'old3', 'new3'};
% cfg.editmarkerfile.toremove       = {'remove1', 'remove2', 'remove3'};
% cfg.editmarkerfile.toadd          = {'add1', 'add2', 'add3'};
% cfg.editmarkerfile.toadd_color    = {'#00ff00','#00ff00','#00ff00'};
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

MuseStruct_new = MuseStruct_orig;

% get default parameters :
cfg.editmarkerfile              = ft_getopt(cfg, 'editmarkerfile', []);
cfg.editmarkerfile.torename     = ft_getopt(cfg.editmarkerfile, 'torename', {});
cfg.editmarkerfile.toremove     = ft_getopt(cfg.editmarkerfile, 'toremove', {});
cfg.editmarkerfile.toadd        = ft_getopt(cfg.editmarkerfile, 'toadd', {});
cfg.editmarkerfile.toadd_color  = ft_getopt(cfg.editmarkerfile, 'toadd_color', {});

nr_renamed = 0;
nr_removed = 0;
nr_added   = 0;

for ipart = 1 : size(MuseStruct_orig, 2)
    
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        
        % rename marker
        for imarker = 1 : size(cfg.editmarkerfile.torename,1)
            if isfield(MuseStruct_orig{ipart}{idir}.markers, cfg.editmarkerfile.torename{imarker,1})
                nr_renamed = nr_renamed + size(MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,1}).synctime, 2);
                
                % add to existing markers if already exists
                if isfield(MuseStruct_new{ipart}{idir}.markers, cfg.editmarkerfile.torename{imarker,2})   
                    MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}).synctime = [MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}).synctime, MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,1}).synctime];
                    MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}).clock    = [MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}).clock,    MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,1}).clock];                    
                else
                    MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}) = MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,1});
                end
                MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers, cfg.editmarkerfile.torename{imarker,1});
            end
        end
        
        % remove marker
        for imarker = 1 : size(cfg.editmarkerfile.toremove,2)
            if isfield(MuseStruct_orig{ipart}{idir}.markers,cfg.editmarkerfile.toremove{imarker})
                nr_removed = nr_removed + 1;%size(MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.toremove{imarker,1}).synctime, 2);    
                
                MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,cfg.editmarkerfile.toremove{imarker});
            end
        end
        
        % add marker
        for imarker = 1 : size(cfg.editmarkerfile.toadd,2)
            if ~isfield(MuseStruct_orig{ipart}{idir}.markers,cfg.editmarkerfile.toadd{imarker})
                nr_added = nr_added + 1;
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).events          = 0;
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).comment         = 'created with editMuseMarkers.m';
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).editable        = 'Yes';
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).classid         = sprintf('+%d',numel(fieldnames(MuseStruct_new{ipart}{idir}.markers)));
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).classgroupid    = '+3';
                if isempty(cfg.editmarkerfile.toadd_color)
                    MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).color       = '00ff00'; %default Muse marker color
                else
                    MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).color       = cfg.editmarkerfile.toadd_color{imarker};
                end
            end
        end
    end
end

ft_info('Renamed %d markers\nRemoved %d markers\nAdded %d (empty) markers\n', nr_renamed, nr_removed, nr_added);