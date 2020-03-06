function [MuseStruct] = addMuseBAD(MuseStruct, ipart, idir, markerStart, markerEnd, isample, time_from_begin, time_from_end)
% add BAD__START__ at markerStart(isample) with a delay of time_from_begin 
% add BAD__END__ at markerEnd(isample) with a delay of time_from_end
% isample can be a single index or an array of indexes
% isample can be 'all' for all samples of the markerStart/markerEnd
% time_before and time_after are expressed in seconds. Can be positive or
% negative
% ipart and idir : possible to input 'all' to consider all the parts and
% dirs


if strcmp(ipart, 'all')
    part_list = 1:length(MuseStruct);
else
    part_list = ipart;
end

if strcmp(isample, 'all')
    isample = ':';
end

for ipart = part_list
    
    if strcmp(idir, 'all')
        dir_list = 1:length(MuseStruct{ipart});
    else
        dir_list = idir;
    end
    
    for idir = dir_list
        
        if isfield(MuseStruct{ipart}{idir}.markers, markerStart)
            if isfield(MuseStruct{ipart}{idir}.markers, markerEnd)
                if length(MuseStruct{ipart}{idir}.markers.(markerStart)) - length(MuseStruct{ipart}{idir}.markers.(markerEnd)) == 0
                    
                    n_BAD_to_add{ipart}{idir} = length(MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample));
                    
                    %Bad start
                    % if there are already BAD__START__ markers, add it to
                    % the end
                    if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1:end+n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1:end+n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                        
                    else
                        % if there are no BAD__START markers, create field
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock = datetime.empty;
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime = [];
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.clock(1:n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(1:n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                        
                    end
                    
                    
                    %Bad end
                    if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__END__')
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1:end+n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).clock(isample)  +  seconds(time_from_end);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1:end+n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(isample)  +  time_from_end;
                        
                    else
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock = datetime.empty;
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime = [];
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.clock(1:n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).clock(isample)  +  seconds(time_from_end);
                        
                        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(1:n_BAD_to_add{ipart}{idir}) = ...
                            MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(isample)  +  time_from_end;
                    end
                else
                    error('Not same amout of %s and %s markers in part %d dir %d', markerStart, markerEnd, ipart, idir);
                    %this could not be an error. I let this as an error for
                    %now, to be sure that if it happens, it is voluntarily
                end
            end
        end
        
    end %idir
    %fprintf('%d BAD Start and End added for part %d\n', sum(cell2mat(n_BAD_to_add{ipart})), ipart);
end

end

