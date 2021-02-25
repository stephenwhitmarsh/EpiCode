function [MuseStruct_updated] = updateMarkers(cfg, MuseStruct_updated, labels)

% don't write original muse markers to file
cfg.muse.write     = false;
MuseStruct_new     = readMuseMarkers(cfg, true);

for ipart = 1 : size(MuseStruct_new, 2)
    for idir = 1 : size(MuseStruct_new{ipart}, 2)
        for label = labels
            try
                MuseStruct_updated{ipart}{idir}.markers.BAD__START__  = MuseStruct_new{ipart}{idir}.markers.BAD__START__;
                MuseStruct_updated{ipart}{idir}.markers.BAD__END__    = MuseStruct_new{ipart}{idir}.markers.BAD__END__;
            catch
            end
        end
    end
end