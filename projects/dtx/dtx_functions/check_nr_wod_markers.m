function check_nr_wod_markers(cfg, MuseStruct, ipart)
%check that their are the same number of each marker which should be
%check that ponctual marker are set as ponctual, and not 'START' 'END'

for idir = 1:length(MuseStruct{ipart})

    if isfield(MuseStruct{ipart}{idir}.markers.Vent_Off,'synctime')
        try
            if isfield(MuseStruct{ipart}{idir}.markers.Vent_On,'synctime')
                if length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime) ~= length(MuseStruct{ipart}{idir}.markers.Vent_On.synctime)
                    fprintf('dir %d %s, %d Vent_Off and %d Vent_On\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime), length(MuseStruct{ipart}{idir}.markers.Vent_On.synctime));
                end
            end
        catch
            fprintf('dir %d %s, %d Vent_Off and 0 Vent_On\n', idir, cfg.directorylist{ipart}{idir}, size(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime,2));
        end
        
        try
            if isfield(MuseStruct{ipart}{idir}.markers.WOD,'synctime')
                if length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime) ~= length(MuseStruct{ipart}{idir}.markers.WOD.synctime)
                    fprintf('dir %d %s, %d Vent_Off and %d WOD\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime), length(MuseStruct{ipart}{idir}.markers.WOD.synctime));
                end
            end
        catch
            fprintf('dir %d %s, %d Vent_Off and 0 WOD\n', idir, cfg.directorylist{ipart}{idir}, size(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime,2));
        end
        
        try
            if isfield(MuseStruct{ipart}{idir}.markers.WOR,'synctime')
                if length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime) ~= length(MuseStruct{ipart}{idir}.markers.WOR.synctime)
                    fprintf('dir %d %s, %d Vent_Off and %d WOR\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime), length(MuseStruct{ipart}{idir}.markers.WOR.synctime));
                end
            end
        catch
            fprintf('dir %d %s, %d Vent_Off and 0 WOR\n', idir, cfg.directorylist{ipart}{idir}, size(MuseStruct{ipart}{idir}.markers.Vent_Off.synctime,2));
        end
        
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'BAD__START__')
        if isfield(MuseStruct{ipart}{idir}.markers.BAD__START__,'synctime')
            try
                if length(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime) ~= length(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime)
                    fprintf('dir %d %s, %d BAD__START__ and %d BAD__END__\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime), length(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime));
                end
            catch
                fprintf('dir %d %s, %d BAD__START__ and 0 BAD__END__\n', idir, cfg.directorylist{ipart}{idir}, size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2));
            end
        end
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'WOD__START__')
        fprintf('dir %d %s, marker error WOD__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'WOR__START__')
        fprintf('dir %d %s, marker error WOR__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'Vent_On__START__')
        fprintf('dir %d %s, marker error Vent_On__START____START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'Vent_Off__START__')
        fprintf('dir %d %s, marker error Vent_Off__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
end

end