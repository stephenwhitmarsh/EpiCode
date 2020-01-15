function muse2SC

%% Place artefact file for Spyking-Circus in each data directory
for idir = 1 : size(MuseStruct,2)
    Starttime = [];
    Endtime = [];
    
    try
        if size(MuseStruct{idir}.markers.BAD__START__.events,2) > 0
            for ievent = 1 : size(MuseStruct{idir}.markers.BAD__START__.events,2)
                Starttime   = [Starttime;   MuseStruct{idir}.markers.BAD__START__.events(ievent).time];
                Endtime     = [Endtime;     MuseStruct{idir}.markers.BAD__END__.events(ievent).time];
            end
            
            filename = fullfile(MuseStruct{idir}.directory,'artifact.txt');
            fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
            
            % has to be in ms instead of seconds
            dlmwrite(filename,[Starttime*1000,Endtime*1000],'delimiter',' ','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense
        else
            fprintf('No artifact markers found in: %s\n',filename);
        end
    catch
        fprintf('No markerfile found in: %s\n',filename);
        
    end
end