function [MuseStruct_updated] = updateMarkers(cfg, MuseStruct_updated, labels)

% don't write original muse markers to file
cfg.muse.write     = false;
MuseStruct_new     = readMuseMarkers(cfg, true);

for ipart = 1 : size(MuseStruct_new, 2)
    for idir = 1 : size(MuseStruct_new{ipart}, 2)
        for label = string(labels)
            if isfield(MuseStruct_new{ipart}{idir}.markers, label)
                MuseStruct_updated{ipart}{idir}.markers.(label)  = MuseStruct_new{ipart}{idir}.markers.(label);
                MuseStruct_updated{ipart}{idir}.markers.(label)  = MuseStruct_new{ipart}{idir}.markers.(label);
                fprintf('Replacing %s markers\n', label);
            else
                fprintf('No %s markers found\n', label);
            end
        end
    end
end