config = hspike_setparams;
old_name = 'Hspike2';
new_name = 'Hspike';

for ipatient = [1, 2]
    % backup markerfiles
    for ipart = 1 : 1
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
    
    %     for ipart = 2 : size(MuseStruct_orig, 2)
    ipart = 1;
    
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        MuseStruct_new{ipart}{idir} = MuseStruct_orig{ipart}{idir};
        
        % rename marker
        if isfield(MuseStruct_orig{ipart}{idir}.markers, old_name)
            MuseStruct_new{ipart}{idir}.markers.(new_name) = MuseStruct_orig{ipart}{idir}.markers.(old_name);
            MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers, old_name);
        end
        
        %             % remove marker
        %             if isfield(MuseStruct_orig{ipart}{idir}.markers,'SpikeDetect')
        %                 MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,'SpikeDetect');
        %             end
        
        % write markers
        fname_mrk = fullfile(config{ipatient}.rawdir, MuseStruct_new{ipart}{idir}.directory,'Events.mrk');
        writeMuseMarkerfile(MuseStruct_new{ipart}{idir}, fname_mrk);
        %         end
    end
end