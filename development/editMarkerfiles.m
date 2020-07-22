
ipatient = 2;

% backup markerfiles
for ipart = 1 : 3
    for idir = 1 : size(config{ipatient}.directorylist{ipart},2)
        
        fname_mrk    = fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir},'Events.mrk');
        
        % backup markerfile
        if ~exist(config{ipatient}.muse.backupdir,'DIR')
            error('Backup directory does not exist');
        end
        [~, d] = fileparts(config{ipatient}.directorylist{ipart}{idir});
        if ~exist(fullfile(config{ipatient}.muse.backupdir, d), 'DIR')
            fprintf('Creating directory: %s\n', fullfile(config{ipatient}.muse.backupdir, d));
            eval(sprintf('!mkdir %s', fullfile(config{ipatient}.muse.backupdir, d)));
        end
        fname_backup = sprintf('Events_%s.mrk', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
        eval(sprintf('!cp %s %s', fname_mrk, fullfile(config{ipatient}.muse.backupdir, d, fname_backup)));
        fprintf('Succesfully backed up markerfile to %s\n',fullfile(config{ipatient}.muse.backupdir, d, fname_backup));
    end
end

[MuseStruct_orig] = readMuseMarkers(config{ipatient}, true);

for ipart = 1 : size(MuseStruct_orig, 2)
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        MuseStruct_new{ipart}{idir} = MuseStruct_orig{ipart}{idir};
        
        % rename marker
        if isfield(MuseStruct_orig{ipart}{idir}.markers,'SpikeHaT1_1')
            MuseStruct_new{ipart}{idir}.markers.Hspike = MuseStruct_orig{ipart}{idir}.markers.SpikeHaT1_1;
            MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,'SpikeHaT1_1');
        end
        
        % remove marker
        if isfield(MuseStruct_orig{ipart}{idir}.markers,'SpikeDetect')
            MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,'SpikeDetect');
        end
        fname_mrk = fullfile(config{ipatient}.rawdir, MuseStruct_new{ipart}{idir}.directory,'Events.mrk');
        writeMuseMarkerfile(MuseStruct_new{ipart}{idir}, fname_mrk);
    end
end