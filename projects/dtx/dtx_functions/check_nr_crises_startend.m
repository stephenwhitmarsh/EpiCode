function check_nr_crises_startend(cfg, MuseStruct, ipart)

for idir = 1:length(MuseStruct{ipart})
    
    if isfield(MuseStruct{ipart}{idir}.markers.Crise_Start,'synctime') && isfield(MuseStruct{ipart}{idir}.markers.Crise_End,'synctime')
        if length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime) ~= length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime)
            fprintf('dir %d %s, %d Crise_Start and %d Crise_End\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime), length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime));
        end
    end
    
    if isfield(MuseStruct{ipart}{idir}.markers.SlowWave_begin,'synctime')
        try
            if isfield(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__,'synctime')
                if length(MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime) ~= length(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__.synctime)
                    fprintf('dir %d %s, %d SlowWave_begin and %d SlowWave_EMG__START__\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime), length(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__.synctime));
                end
            end
                catch
                    fprintf('dir %d %s, %d SlowWave_begin and 0 SlowWave_EMG__START__\n', idir, cfg.directorylist{ipart}{idir}, size(MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime,2));
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
    
    if isfield(MuseStruct{ipart}{idir}.markers,'SlowWave_begin__START__')
        fprintf('dir %d %s, marker error SlowWave_begin__START__\n', idir, cfg.directorylist{ipart}{idir});
    end
    
end

end