function [MuseStruct] = addMuseBAD(MuseStruct, part_list, dir_list, markerStart, markerEnd, sample_list, time_from_begin, time_from_end)
% add BAD__START__ at markerStart(isample) with a delay of time_from_begin
% add BAD__END__ at the next marker markerEnd (>= to markerstart) with a delay of time_from_end
% if there is no markerend after markerstart(isample) in the dir,
% BAD__END__ is put at the end of the dir
% isample can be a single index or an array of indexes
% isample can be 'all' for all samples of the markerStart/markerEnd
% time_before and time_after are expressed in seconds. Can be positive or
% negative
% ipart and idir : can be 'all'

if strcmp(part_list, 'all')
    part_list = 1:length(MuseStruct);
end

for ipart = part_list
    
    if strcmp(dir_list, 'all')
        dir_list = 1:length(MuseStruct{ipart});
    end
    
    for idir = dir_list
        
        if isfield(MuseStruct{ipart}{idir}.markers, markerStart)
            if isfield(MuseStruct{ipart}{idir}.markers, markerEnd)
                if isfield(MuseStruct{ipart}{idir}.markers.(markerStart), 'clock') && isfield(MuseStruct{ipart}{idir}.markers.(markerEnd), 'clock')
                    
                    
                    if strcmp(sample_list, 'all')
                        samples = 1:size(MuseStruct{ipart}{idir}.markers.(markerStart).clock,2);
                    else
                        samples = sample_list;
                    end
                    
                    for isample = samples
                        
                        %find idx of markerEnd
                        start   = round(MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample));
                        end_idx = find(round(MuseStruct{ipart}{idir}.markers.(markerEnd).synctime) >= start,1,'first');
                        
                        %if there is BAD__START__ marker, add sample to the end
                        if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                            
                            MuseStruct{ipart}{idir}.markers.BAD__START__.clock(end+1) = ...
                                MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                            
                            MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(end+1) = ...
                                MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                            
                        else
                            % if there are no BAD__START markers, create field
                            MuseStruct{ipart}{idir}.markers.BAD__START__.clock = datetime.empty;
                            MuseStruct{ipart}{idir}.markers.BAD__START__.synctime = [];
                            
                            MuseStruct{ipart}{idir}.markers.BAD__START__.clock(1) = ...
                                MuseStruct{ipart}{idir}.markers.(markerStart).clock(isample)  +  seconds(time_from_begin);
                            
                            MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(1) = ...
                                MuseStruct{ipart}{idir}.markers.(markerStart).synctime(isample)  +  time_from_begin;
                            
                        end
                        
                        %Bad end
                        if isfield(MuseStruct{ipart}{idir}.markers, 'BAD__END__')
                            
                            if ~isempty(end_idx)
                                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1) = ...
                                    MuseStruct{ipart}{idir}.markers.(markerEnd).clock(end_idx)  +  seconds(time_from_end);
                                
                                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1) = ...
                                    MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(end_idx)  +  time_from_end;
                            else
                                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(end+1)     = ...
                                    MuseStruct{ipart}{idir}.endtime;
                                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(end+1)  = ...
                                    seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
                            end
                                
                        else
                            
                            MuseStruct{ipart}{idir}.markers.BAD__END__.clock = datetime.empty;
                            MuseStruct{ipart}{idir}.markers.BAD__END__.synctime = [];
                            
                            if ~isempty(end_idx)
                                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(1) = ...
                                    MuseStruct{ipart}{idir}.markers.(markerEnd).clock(end_idx)  +  seconds(time_from_end);
                                
                                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(1) = ...
                                    MuseStruct{ipart}{idir}.markers.(markerEnd).synctime(end_idx)  +  time_from_end;
                            else
                                MuseStruct{ipart}{idir}.markers.BAD__END__.clock(1)     = ...
                                    MuseStruct{ipart}{idir}.endtime;
                                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(1)  = ...
                                    seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);                            
                            end
                                
                        end
                    end
                end
            end
        end
        
    end %idir
    %fprintf('%d BAD Start and End added for part %d\n', sum(cell2mat(n_BAD_to_add{ipart})), ipart);
end

end

