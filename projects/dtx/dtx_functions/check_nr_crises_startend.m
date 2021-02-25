function check_nr_crises_startend(cfg, MuseStruct, ipart)

for idir = 1:length(MuseStruct{ipart})
    
    %be sure all markers are in MuseStruct
    marker.synctime = [];
    MuseStruct{ipart}{idir}.markers.SlowWave              = ft_getopt(MuseStruct{ipart}{idir}.markers,'SlowWave',        marker);
    MuseStruct{ipart}{idir}.markers.Crise_NoSW            = ft_getopt(MuseStruct{ipart}{idir}.markers,'Crise_NoSW',        marker);
    MuseStruct{ipart}{idir}.markers.SlowWave_begin        = ft_getopt(MuseStruct{ipart}{idir}.markers,'SlowWave_begin',        marker);
    MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__ = ft_getopt(MuseStruct{ipart}{idir}.markers,'SlowWave_EMG__START__', marker);
    MuseStruct{ipart}{idir}.markers.SlowWave_EMG__END__   = ft_getopt(MuseStruct{ipart}{idir}.markers,'SlowWave_EMG__END__',   marker);
    MuseStruct{ipart}{idir}.markers.Crise_Start           = ft_getopt(MuseStruct{ipart}{idir}.markers,'Crise_Start',   marker);
    MuseStruct{ipart}{idir}.markers.Crise_End             = ft_getopt(MuseStruct{ipart}{idir}.markers,'Crise_End',   marker);
    %be sure synctime is in all markers
    MuseStruct{ipart}{idir}.markers.SlowWave.synctime              = ft_getopt(MuseStruct{ipart}{idir}.markers.SlowWave,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.Crise_NoSW.synctime            = ft_getopt(MuseStruct{ipart}{idir}.markers.Crise_NoSW,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime        = ft_getopt(MuseStruct{ipart}{idir}.markers.SlowWave_begin,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__.synctime = ft_getopt(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.SlowWave_EMG__END__.synctime   = ft_getopt(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__END__,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.Crise_Start.synctime           = ft_getopt(MuseStruct{ipart}{idir}.markers.Crise_Start,'synctime',[]);
    MuseStruct{ipart}{idir}.markers.Crise_End.synctime             = ft_getopt(MuseStruct{ipart}{idir}.markers.Crise_End,'synctime',[]);
    
    %same amount of crise start and end
    if length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime) ~= length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime)
        fprintf('dir %d %s, %d Crise_Start and %d Crise_End\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime), length(MuseStruct{ipart}{idir}.markers.Crise_End.synctime));
    end
    
    %same amount of omg and sw_begin
    if length(MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime) ~= length(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__.synctime)
        fprintf('dir %d %s, %d SlowWave_begin and %d SlowWave_EMG__START__\n', idir, cfg.directorylist{ipart}{idir}, length(MuseStruct{ipart}{idir}.markers.SlowWave_begin.synctime), length(MuseStruct{ipart}{idir}.markers.SlowWave_EMG__START__.synctime));
    end

    %même quantité de sw+crise_nosw que de crise_start
    nb_crises_1 = length(MuseStruct{ipart}{idir}.markers.SlowWave.synctime) + length(MuseStruct{ipart}{idir}.markers.Crise_NoSW.synctime);
    nb_crises_2 = length(MuseStruct{ipart}{idir}.markers.Crise_Start.synctime);
    if nb_crises_1 ~= nb_crises_2
        fprintf('dir %d %s, %d Crise_Start and %d SlowWave+Crise_NoSW\n', idir, cfg.directorylist{ipart}{idir}, nb_crises_2, nb_crises_1);
    end
    
    %erreur avec markers sur une plage au lieu de markers ponctuels
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