function [MuseStruct] = addMuseBAD(MuseStruct, part_list, dir_list, markerStart, markerEnd, sample_list, time_from_begin, time_from_end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [MuseStruct] = addMuseBAD(MuseStruct, part_list, dir_list, markerStart, markerEnd, sample_list, time_from_begin, time_from_end)
% 
% Converts Muse marker timings into BAD timings. 
% For example, convert seizures identified by 'Crise_Start' and 'Crise_End' 
% markers into BAD__START__ and BAD__END__ markers, so those periods are
% then considered as artefacts, and can be ignored with the appropriate
% scripts (writeSpykingCircus.m,  removetrials_MuseMarkers.m)
%
% ### INPUT
% MuseStruct        : structure with all the marker timings created by Muse 
%                     (see readMuseMarkers.m)
% part_list         : list of parts to analyse. Can be 'all', 'last', or any
%                     array of indexes.
% dir_list          : list of directories to analyse. Can be 'all', 'last', 
%                     or any array of indexes (note that in this case, the 
%                     input indexes must exist in all the 'parts' of 
%                     MuseStruct).
% markerStart       : name of the Muse marker used to identify timings to 
%                     convert into BAD__START__. Can be 'begin', in this 
%                     case the BAD__START__ marker is put at the begining 
%                     of the dir.
% markerEnd         : name of the Muse marker used to identify timings to 
%                     convert into BAD__END__. The marker used is the next
%                     occurence of markerEnd after the selected markerStart.
%                     Can be 'end', in this case the BAD__END__ marker is 
%                     put at the end of the dir.
% sample_list       : indication of the samples indexes of markerStart to 
%                     use to add BAD markers. Can be 'all', 'last', or any 
%                     array of indexes (note that in this  case, the samples 
%                     indexes must exist in all the parts and dirs selected)
% time_from_begin   : delay from markerStart where the new BAD__START__
%                     markers will be defined. In seconds, can be positive 
%                     or negative.
% time_from_end     : delay from markerEnd where the new BAD__END__
%                     markers will be defined. In seconds, can be positive 
%                     or negative.
% 
% ## OUTPUT :
% MuseStruct        : The input structure, with the requested BAD markers
%                     added.
% 
% Notes :
% - if no markerEnd is found after markerStart, then the end of the file is
%   used to put the BAD__END__ marker
% - if markersStart = 'begin' and no markerEnd is found on all the file,
%   nothing is done
% - it is not possible for now to put at the same time markerStart =
%   'begin' and markeEnd = 'end'
% - all BAD markers are sorted in ascending order after all were added, to
%   avoid bug with Spyking-Circus
% 
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(markerStart, 'begin') && strcmp(markerEnd, 'end')
    error('It is not possible for now to put at the same time markerStart = ''begin'' and markeEnd = ''end''');
end

if strcmp(part_list, 'all')
    part_list = 1:length(MuseStruct);
elseif strcmp(part_list, 'last')
    part_list = length(MuseStruct);
end

for ipart = part_list
    
    if strcmp(dir_list, 'all')
        dir_list = 1:length(MuseStruct{ipart});
    elseif strcmp(dir_list, 'last')
        dir_list = length(MuseStruct{ipart});
    end
    
    for idir = dir_list
        
        %find start samples
        samples = [];
        if strcmp(sample_list, 'all') || strcmp(sample_list, 'last') %find all samples if required
            if isfield(MuseStruct{ipart}{idir}.markers, markerStart)
                if isfield(MuseStruct{ipart}{idir}.markers.(markerStart), 'clock')
                    switch sample_list
                        case 'all'
                            samples = 1:size(MuseStruct{ipart}{idir}.markers.(markerStart).clock,2);
                        case 'last'
                            samples = size(MuseStruct{ipart}{idir}.markers.(markerStart).clock,2);
                    end
                end
            end
        else
            samples = sample_list; %take input sample list if required
        end
        if strcmp(markerStart, 'begin'), samples = 1; end %if markerstart = begin, use sample = 1
        
        for isample = samples
            
            %find idx of markerEnd
            if strcmp(markerStart, 'begin')
                start = 0;
            else
                start   = round(MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample));
            end
            
            end_idx = [];
            if isfield(MuseStruct{ipart}{idir}.markers, markerEnd)
                if isfield(MuseStruct{ipart}{idir}.markers.(markerEnd), 'synctime')
                    end_idx = find(round(MuseStruct{ipart}{idir}.markers.(markerEnd).synctime) >= start,1,'first');
                end
            end
            
            if strcmp(markerStart, 'begin') && isempty(end_idx)
                warning('With markerStart = ''begin'', part %d dir %d : no markerEnd found, nothing is done', ipart, idir);
            else
                
                %Add BAD__START__ and BAD__End__ new timings
                
                %if there is BAD__START__ marker, add sample to the end of the fields
                if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                    
                    if strcmp(markerStart, 'begin')
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1)       = MuseStruct{ipart}{idir}.starttime;
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1)    = 0;
                    else
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                    end
                    
                else
                    % if there are no BAD__START markers, create field
                    MuseStruct{ipart}{idir}.markers.BAD__START__.clock = datetime.empty;
                    MuseStruct{ipart}{idir}.markers.BAD__START__.synctime = [];
                    
                    if strcmp(markerStart, 'begin')
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(1)       = MuseStruct{ipart}{idir}.starttime;
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(1)    = 0;
                    else
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                    end
                end
                
                %Bad end
                if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__END__')
                    
                    if ~isempty(end_idx) && ~strcmp(markerEnd, 'end')
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).clock(end_idx)  +  seconds(time_from_end);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(end_idx)  +  time_from_end;
                    else % if no markerend event, or if markerEnd = 'end', put BAD__END__ at the end of the file
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1)     = ...
                            MuseStruct{ipart}{idir}.endtime;
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1)  = ...
                            seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
                    end
                    
                else
                    
                    MuseStruct{ipart}{idir}.markers.BAD__END__.clock = datetime.empty;
                    MuseStruct{ipart}{idir}.markers.BAD__END__.synctime = [];
                    
                    if ~isempty(end_idx)&& ~strcmp(markerEnd, 'end')
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).clock(end_idx)  +  seconds(time_from_end);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(1) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(end_idx)  +  time_from_end;
                    else % if no markerend event, or if markerEnd = 'end', put BAD__END__ at the end of the file
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(1)     = ...
                            MuseStruct{ipart}{idir}.endtime;
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(1)  = ...
                            seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
                    end
                    
                end
            end
        end%isample
    end %idir
end %ipart

%% sort BAD markers in ascending order to avoid bug in Spyking Circus (06.2020 : will be corrected in the next Spyking Circus release)
for ipart = 1:size(MuseStruct,2)
    for idir = 1:size(MuseStruct{ipart},2)
        if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
            if isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'clock')
                MuseStruct{ipart}{idir}.markers.BAD__START__.clock     = sort(MuseStruct{ipart}{idir}.markers.BAD__START__.clock);
                MuseStruct{ipart}{idir}.markers.BAD__START__.synctime  = sort(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime);
                MuseStruct{ipart}{idir}.markers.BAD__END__.clock       = sort(MuseStruct{ipart}{idir}.markers.BAD__END__.clock);
                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime    = sort(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime);
            end
        end
    end
end

