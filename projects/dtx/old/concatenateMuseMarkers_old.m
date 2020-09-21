function marker = concatenateMuseMarkers_old(MuseStruct, ipart, markerName)
% Concatenate markers over all the dirs of one MuseStruct part.
% First line of the output array : clock or synctime
% Second line of the output array : nr of the dir
% clock time : real time, in datetime format
% synctime : time in seconds since the begining of the first file
% dir : integer number of the data directiry (see cfg.directorylist)

marker.clock    = datetime.empty;
marker.synctime = [];
marker.dir      = [];
length_previous = 0;

%concatenate clock times
for idir = 1:size(MuseStruct{ipart},2) 
    if isfield(MuseStruct{ipart}{idir}.markers, markerName)
        if isfield(MuseStruct{ipart}{idir}.markers.(markerName), 'clock')
            
            marker_temp = MuseStruct{ipart}{idir}.markers.(markerName).clock;
            
            for isample = 1:size(marker_temp,2)
                marker.clock(isample+length_previous) = marker_temp(isample); 
                marker.dir(isample+length_previous) = idir;
            end
            
            length_previous = length_previous + size(marker_temp,2); 
            
            
        else
            %fprintf('Marker %s exists but is empty in %s\n', markerName, MuseStruct{ipart}{idir}.directory);
        end
        
        
    else
        %fprintf('Marker %s does not exist in %s\n',markerName, MuseStruct{ipart}{idir}.directory);
    end
    
end

fprintf('%d times found for %s\n',size(marker.clock,2),markerName);

% convert clock time to synctime
for isample = 1:size(marker.clock,2) 
    marker.synctime(isample) = seconds(marker.clock(isample) - MuseStruct{ipart}{1}.starttime);    
end

end