function [marker_clock, marker_synctime] = concatenateMuseMarkers(cfg, MuseStruct, ipart, markerName)
% Concatenate markers over all the dirs of the MuseStruct part.
% First line of the output array : clock or synctime
% Second line of the output array : nr of the dir

marker_clock = [];
marker_synctime = [];
length_previous = 0;

%concatenate clock times
for idir = 1:size(MuseStruct{ipart},2) 
    if isfield(MuseStruct{ipart}{idir}.markers, markerName)
        if isfield(MuseStruct{ipart}{idir}.markers.(markerName), 'clock')
            
            marker_temp = MuseStruct{ipart}{idir}.markers.(markerName).clock;
            
            for isample = 1:size(marker_temp,2)
                marker_clock{1,isample+length_previous} = marker_temp(isample); 
                marker_clock{2,isample+length_previous} = idir;
            end
            
            length_previous = length_previous + size(marker_temp,2); 
            
            
        else
            %fprintf('Marker %s exists but is empty in %s\n', markerName, MuseStruct{ipart}{idir}.directory);
        end
        
        
    else
        %fprintf('Marker %s does not exist in %s\n',markerName, MuseStruct{ipart}{idir}.directory);
    end
    
end

fprintf('%d times found for %s in %s \n',size(marker_clock,2),markerName, cfg.prefix(1:end-1));

% convert clock time to synctime
for isample = 1:size(marker_clock,2) 
    
    marker_synctime{1, isample} = seconds(marker_clock{1, isample} - MuseStruct{ipart}{1}.starttime);
    marker_synctime{2, isample} = marker_clock{2, isample};
    
end

end


