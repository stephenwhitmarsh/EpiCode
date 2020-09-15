function edit_Markerfiles(cfg)

% Add, remove or rename markers in each Muse marker file of the patient
% 
% ATTENTION : edit_Markerfiles does not work for brainvision data, because there is an other format of Muse marker file
% 
% cfg fields : 
% * Optionals :
% cfg.editmarkerfile.torename = {'Baseline_Start', 'Analysis_Start'}; %{old1, new1; old2, new2; old3, new3};
% cfg.editmarkerfile.toremove = {}; %{remove1, remove2, remove3};
% cfg.editmarkerfile.toadd    = {'SlowWave_begin','SlowWave_R_begin','SlowWave_L_begin'}; %{add1, add2, add3};
% cfg.editmarkerfile.toadd_color = {'#00ff00','#00ff00','#00ff00'};%color of the new marker on muse
% * Mandatory : 
% cfg.rawdir
% cfg.directorylist
% cfg.muse.backupdir
% cfg.prefix
% 
% Dependencies on other EpiCode functions :
% - readMuseMarkers.m
% - writeMuseMarkerfile.m

%get default parameters : 
cfg.editmarkerfile              = ft_getopt(cfg, 'editmarkerfile', []);
cfg.editmarkerfile.torename     = ft_getopt(cfg.editmarkerfile, 'torename', {});
cfg.editmarkerfile.toremove     = ft_getopt(cfg.editmarkerfile, 'toremove', {});
cfg.editmarkerfile.toadd        = ft_getopt(cfg.editmarkerfile, 'toadd', {});
cfg.editmarkerfile.toadd_color  = ft_getopt(cfg.editmarkerfile, 'toadd_color', {});

if isempty(cfg.editmarkerfile.toadd_color) 
    for imarker = 1 : size(cfg.editmarkerfile.toadd,2)
        cfg.editmarkerfile.toadd_color{imarker} = '00ff00'; %default Muse marker color
    end
end
    
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

if isBrainvision
    error(sprintf('This script does not support Brainvision''s Muse marker file format. \nTo use it with this format, you need to adapt ''writeMuseMarkerfile.m'' for the ''.vmrk'' format and remove this error message.'));
end

% backup markerfiles
for ipart = 1 : size(cfg.directorylist,2)
    for idir = 1 : size(cfg.directorylist{ipart},2)
        
        %search marker file
        if isNeuralynx
            fname_mrk_temp    = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir},'Events*.mrk');
        elseif isMicromed %mrk file has the same name as the eeg file
            fname_mrk_temp    = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir},'.mrk']);
        end
        temp                  = dir(fname_mrk_temp);
        if size(temp,1) == 1
            fname_mrk{ipart}{idir}  = fullfile(temp.folder, temp.name);
        else
            error('Error when searching marker file %s', fname_mrk_temp);
        end
        
        % backup markerfile
        if ~exist(cfg.muse.backupdir,'DIR')
            error('Backup directory does not exist');
        end
        [~, d] = fileparts(cfg.directorylist{ipart}{idir});
        dir_backup = fullfile(cfg.muse.backupdir,cfg.prefix(1:end-1), d);
        
        if ~exist(dir_backup, 'DIR')
            fprintf('Creating directory: %s\n', dir_backup);
            %             eval(sprintf('!mkdir %s', fullfile(cfg.muse.backupdir, d)));
            mkdir(dir_backup);
        end
        fname_backup = sprintf('Events_%s.mrk', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
        %         eval(sprintf('!cp %s %s', fname_mrk, fullfile(cfg.muse.backupdir, d, fname_backup)));
        copyfile(fname_mrk{ipart}{idir},fullfile(dir_backup, fname_backup));
        fprintf('Succesfully backed up markerfile to %s\n',fullfile(dir_backup, fname_backup));
    end
end

[MuseStruct_orig] = readMuseMarkers(cfg, true);

for ipart = 1 : size(MuseStruct_orig, 2)
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        MuseStruct_new{ipart}{idir} = MuseStruct_orig{ipart}{idir};
        
        % rename marker
        for imarker = 1 : size(cfg.editmarkerfile.torename,1)
            if isfield(MuseStruct_orig{ipart}{idir}.markers,cfg.editmarkerfile.torename{imarker,1})
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,2}) = MuseStruct_orig{ipart}{idir}.markers.(cfg.editmarkerfile.torename{imarker,1});
                MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,cfg.editmarkerfile.torename{imarker,1});
            end
        end
        
        % remove marker
        for imarker = 1 : size(cfg.editmarkerfile.toremove,2)
            if isfield(MuseStruct_orig{ipart}{idir}.markers,cfg.editmarkerfile.toremove{imarker})
                MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,cfg.editmarkerfile.toremove{imarker});
            end
        end
        
        % add marker
        for imarker = 1 : size(cfg.editmarkerfile.toadd,2)
            if ~isfield(MuseStruct_orig{ipart}{idir}.markers,cfg.editmarkerfile.toadd{imarker})
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).events     = 0;
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).comment    = 'created with edit_Markerfiles.m';
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).color      = cfg.editmarkerfile.toadd_color{imarker};
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).editable   = 'Yes';
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).classid    = sprintf('+%d',numel(fieldnames(MuseStruct_new{ipart}{idir}.markers)));
                MuseStruct_new{ipart}{idir}.markers.(cfg.editmarkerfile.toadd{imarker}).classgroupid = '+3'; %FIXME : I don't know what is this field for
            end
        end
        
        %write new marker file
        %fname_mrk = fullfile(cfg.rawdir, MuseStruct_new{ipart}{idir}.directory,'Events.mrk');
        writeMuseMarkerfile(MuseStruct_new{ipart}{idir}, fname_mrk{ipart}{idir});
        
    end
end
