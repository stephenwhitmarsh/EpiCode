function [MuseStruct] = addMuseBAD(cfg,MuseStruct)

% ADDMUSEBAD converts MuseStruct markers into BAD MuseStruct markers
%
% Use as
%    [MuseStruct] = addMuseBAD(cfg, MuseStruct)
%
% For example, convert seizures identified by 'Crise_Start' and 'Crise_End'
% markers into BAD__START__ and BAD__END__ markers, so those periods are
% then considered as artefacts, and can be ignored with the appropriate
% scripts (writeSpykingCircus.m,  removetrials_MuseMarkers.m)
%
% Mandatory input :
% MuseStruct        : structure with all the marker timings created by Muse
%                     (see readMuseMarkers.m)
%
% ## necessary cfg fields :
% cfg.bad.markerStart : name of the Muse marker used to identify timings to
%                     convert into BAD__START__. Can be 'begin', in this
%                     case the BAD__START__ marker is put at the begining
%                     of the dir.
% cfg.bad.markerEnd : name of the Muse marker used to identify timings to
%                     convert into BAD__END__. The marker used is the next
%                     occurence of cfg.bad.markerEnd after the selected
%                     cfg.bad.markerStart.Can be 'end', in this case the
%                     BAD__END__ marker is put at the end of the dir.
%
% Optional cfg fields :
% cfg.bad.part_list : list of parts to analyse. Can be 'all', 'last', or any
%                     array of indexes. Default = 'all'.
% cfg.bad.dir_list  : list of directories to analyse. Can be 'all', 'last',
%                     or any array of indexes (note that in this case, the
%                     input indexes must exist in all the 'parts' of
%                     MuseStruct). Default = 'all'.
% cfg.bad.sample_list : indication of the samples indexes of cfg.bad.markerStart
%                     to use to add BAD markers. Can be 'all', 'last', or any
%                     array of indexes (note that in this  case, the samples
%                     indexes must exist in all the parts and dirs
%                     selected). Default = 'all'.
% cfg.bad.time_from_begin : delay from cfg.bad.markerStart where the new
%                     BAD__START__ markers will be defined. In seconds, can
%                     be positive or negative. Default = 0.
% cfg.bad.time_from_end : delay from cfg.bad.markerEnd where the new BAD__END__
%                     markers will be defined. In seconds, can be positive
%                     or negative. Default = 0.
%
% ## OUTPUT :
% MuseStruct        : The input structure, with the requested BAD markers
%                     added.
%
% Notes :
% - if no cfg.bad.markerEnd is found after cfg.bad.markerStart, then the end
%   of the file is used to put the BAD__END__ marker
% - if markersStart = 'begin' and no cfg.bad.markerEnd is found on all the
%   file, nothing is done
% - it is not possible for now to put at the same time cfg.bad.markerStart =
%   'begin' and markeEnd = 'end'
% - all BAD markers are sorted in ascending order after all were added, to
%   avoid bug with Spyking-Circus

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

% get the default options
cfg.bad.part_list               = ft_getopt(cfg.bad, 'part_list', 'all');
cfg.bad.dir_list                = ft_getopt(cfg.bad, 'dir_list', 'all');
cfg.bad.sample_list             = ft_getopt(cfg.bad, 'sample_list', 'all');
cfg.bad.time_from_begin         = ft_getopt(cfg.bad, 'time_from_begin', 0);
cfg.bad.time_from_end           = ft_getopt(cfg.bad, 'time_from_end', 0);

if strcmp(cfg.bad.markerStart, 'begin') && strcmp(cfg.bad.markerEnd, 'end')
    error('It is not possible for now to put at the same time cfg.bad.markerStart = ''begin'' and markeEnd = ''end''');
end

if strcmp(cfg.bad.part_list, 'all')
    cfg.bad.part_list = 1:length(MuseStruct);
elseif strcmp(cfg.bad.part_list, 'last')
    cfg.bad.part_list = length(MuseStruct);
end

for ipart = cfg.bad.part_list
    
    if strcmp(cfg.bad.dir_list, 'all')
        cfg.bad.dir_list = 1:length(MuseStruct{ipart});
    elseif strcmp(cfg.bad.dir_list, 'last')
        cfg.bad.dir_list = length(MuseStruct{ipart});
    end
    
    for idir = cfg.bad.dir_list
        
        %find start samples
        samples = [];
        if strcmp(cfg.bad.sample_list, 'all') || strcmp(cfg.bad.sample_list, 'last') %find all samples if required
            if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.bad.markerStart)
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart), 'clock')
                continue
            end
            switch cfg.bad.sample_list
                case 'all'
                    samples = 1:size(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart).clock,2);
                case 'last'
                    samples = size(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart).clock,2);
            end
        else
            samples = cfg.bad.sample_list; %take input sample list if required
        end
        if strcmp(cfg.bad.markerStart, 'begin'), samples = 1; end %if cfg.bad.markerStart = begin, use sample = 1
        
        for isample = samples
            
            %find idx of cfg.bad.markerEnd
            if strcmp(cfg.bad.markerStart, 'begin')
                start = 0;
            else
                start   = round(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart).synctime(isample));
            end
            
            end_idx = [];
            if isfield(MuseStruct{ipart}{idir}.markers, cfg.bad.markerEnd)
                if isfield(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerEnd), 'synctime')
                    end_idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.bad.markerEnd).synctime) >= start,1,'first');
                end
            end
            
            if strcmp(cfg.bad.markerStart, 'begin') && isempty(end_idx)
                warning('With cfg.bad.markerStart = ''begin'', part %d dir %d : no cfg.bad.markerEnd found, nothing is done', ipart, idir);
                continue
            end
            
            %create BAD start and end field if they do not exist
            MuseStruct{ipart}{idir}.markers.BAD__START__            = ft_getopt(MuseStruct{ipart}{idir}.markers, 'BAD__START__', []);
            MuseStruct{ipart}{idir}.markers.BAD__END__              = ft_getopt(MuseStruct{ipart}{idir}.markers, 'BAD__END__', []);
            MuseStruct{ipart}{idir}.markers.BAD__START__.clock      = ft_getopt(MuseStruct{ipart}{idir}.markers.BAD__START__, 'clock', datetime.empty);
            MuseStruct{ipart}{idir}.markers.BAD__END__.clock        = ft_getopt(MuseStruct{ipart}{idir}.markers.BAD__END__, 'clock', datetime.empty);
            MuseStruct{ipart}{idir}.markers.BAD__START__.synctime   = ft_getopt(MuseStruct{ipart}{idir}.markers.BAD__START__, 'synctime', []);
            MuseStruct{ipart}{idir}.markers.BAD__END__.synctime     = ft_getopt(MuseStruct{ipart}{idir}.markers.BAD__END__, 'synctime', []);
            
            %Bad start
            if strcmp(cfg.bad.markerStart, 'begin')
                MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1)       = MuseStruct{ipart}{idir}.starttime;
                MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1)    = 0;
            else
                MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1) = ...
                    MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart).clock(isample)  +  seconds(cfg.bad.time_from_begin);
                
                MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1) = ...
                    MuseStruct{ipart}{idir}.markers.(cfg.bad.markerStart).synctime(isample)  +  cfg.bad.time_from_begin;
            end
            
            %Bad end
            if ~isempty(end_idx) && ~strcmp(cfg.bad.markerEnd, 'end')
                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1) = ...
                    MuseStruct{ipart}{idir}.markers.(cfg.bad.markerEnd).clock(end_idx)  +  seconds(cfg.bad.time_from_end);
                
                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1) = ...
                    MuseStruct{ipart}{idir}.markers.(cfg.bad.markerEnd).synctime(end_idx)  +  cfg.bad.time_from_end;
            else % if no cfg.bad.markerEnd event, or if cfg.bad.markerEnd = 'end', put BAD__END__ at the end of the file
                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1)     = ...
                    MuseStruct{ipart}{idir}.endtime;
                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1)  = ...
                    seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
            end
            
        end%isample
    end %idir
end %ipart

%% sort BAD markers in ascending order to avoid bug in Spyking Circus (06.2020 : will be corrected in the next Spyking Circus release)
for ipart = 1:size(MuseStruct,2)
    for idir = 1:size(MuseStruct{ipart},2)
        
        if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
            continue
        end
        
        if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'clock')
            continue
        end
        
        MuseStruct{ipart}{idir}.markers.BAD__START__.clock     = sort(MuseStruct{ipart}{idir}.markers.BAD__START__.clock);
        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime  = sort(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime);
        MuseStruct{ipart}{idir}.markers.BAD__END__.clock       = sort(MuseStruct{ipart}{idir}.markers.BAD__END__.clock);
        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime    = sort(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime);
        
    end
end
