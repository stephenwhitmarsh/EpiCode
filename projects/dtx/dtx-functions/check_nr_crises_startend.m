function check_nr_crises_startend(cfg, MuseStruct, ipart)

for idir = 1:length(MuseStruct{ipart})
    
    if isfield(MuseStruct{ipart}{idir}.markers.Crise_Start,'synctime') && isfield(MuseStruct{ipart}{idir}.markers.Crise_End,'synctime')
        if ~(length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime) == length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime))
            fprintf('dir %d %s, %d Crise_Start and %d Crise_End\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime), length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime));
        end
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'Crise_Start__START__')
        fprintf('dir %d %s, marker error Crise_Start__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'Crise_End__START__')
        fprintf('dir %d %s, marker error Crise_End__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'SlowWave__START__')
        fprintf('dir %d %s, marker error SlowWave__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers,'SlowWave_EMG__START__')
        fprintf('dir %d %s, marker error SlowWave_EMG__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
end

end