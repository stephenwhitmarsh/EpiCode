% edit_Markerfiles.m

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setting parameters : 
config   = dtx_setparams_patients_lgi1;
pat_list = 1:length(config);
torename = {'Baseline_Start', 'Analysis_Start'}; %{old1, new1; old2, new2; old3, new3};
toremove = {}; %{remove1, remove2, remove3};
toadd    = {'SlowWave_R_EMG','SlowWave_L_EMG'};%{'SlowWave_begin','SlowWave_R_begin','SlowWave_L_begin'}; %{add1, add2, add3};
toadd_color = {'#00ff00','#00ff00','#00ff00'};%color of the new marker on muse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ipatient = 12%1:size(config,2)
    
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(config{ipatient});
    
    % backup markerfiles
    for ipart = 1 : size(config{ipatient}.directorylist,2)
        for idir = 1 : size(config{ipatient}.directorylist{ipart},2)
            
            %search marker file
            if isNeuralynx
                fname_mrk_temp    = fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir},'Events*.mrk');
            elseif isMicromed %mrk file has the same name as the eeg file
                fname_mrk_temp    = fullfile(config{ipatient}.rawdir,[config{ipatient}.directorylist{ipart}{idir},'.mrk']);
			elseif isBrainvision %mrk file has the same name as the eeg file
                fname_mrk_temp    = fullfile(config{ipatient}.rawdir,[config{ipatient}.directorylist{ipart}{idir},'.vmrk']);
            end
            temp                  = dir(fname_mrk_temp);
            if size(temp,1) == 1
                fname_mrk{ipart}{idir}  = fullfile(temp.folder, temp.name);
            else
                error('Error when searching marker file %s', fname_mrk_temp);
            end
            
            % backup markerfile
            if ~exist(config{ipatient}.muse.backupdir,'DIR')
                error('Backup directory does not exist');
            end
            [~, d] = fileparts(config{ipatient}.directorylist{ipart}{idir});
            dir_backup = fullfile(config{ipatient}.muse.backupdir,config{ipatient}.prefix(1:end-1), d);
            
            if ~exist(dir_backup, 'DIR')
                fprintf('Creating directory: %s\n', dir_backup);
                mkdir(dir_backup);
            end
            fname_backup = sprintf('Events_%s.mrk', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            copyfile(fname_mrk{ipart}{idir},fullfile(dir_backup, fname_backup));
            fprintf('Succesfully backed up markerfile to %s\n',fullfile(dir_backup, fname_backup));
        end
    end
    
    [MuseStruct_orig] = readMuseMarkers(config{ipatient}, true);
    
    for ipart = 1 : size(MuseStruct_orig, 2)
        for idir = 1 : size(MuseStruct_orig{ipart}, 2)
            MuseStruct_new{ipart}{idir} = MuseStruct_orig{ipart}{idir};
            
            % rename marker
            for imarker = 1 : size(torename,1)
                if isfield(MuseStruct_orig{ipart}{idir}.markers,torename{imarker,1})
                    MuseStruct_new{ipart}{idir}.markers.(torename{imarker,2}) = MuseStruct_orig{ipart}{idir}.markers.(torename{imarker,1});
                    MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,torename{imarker,1});
                end
            end
            
            % remove marker
            for imarker = 1 : size(toremove,2)
                if isfield(MuseStruct_orig{ipart}{idir}.markers,toremove{imarker})
                    MuseStruct_new{ipart}{idir}.markers = rmfield(MuseStruct_new{ipart}{idir}.markers,toremove{imarker});
                end
            end
            
            % add marker
            for imarker = 1 : size(toadd,2)
                if ~isfield(MuseStruct_orig{ipart}{idir}.markers,toadd{imarker})
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).events     = 0;
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).comment    = 'created with edit_Markerfiles.m';
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).color      = toadd_color{imarker};
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).editable   = 'Yes';
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).classid    = sprintf('+%d',numel(fieldnames(MuseStruct_new{ipart}{idir}.markers)));
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).classgroupid = '+3'; %FIXME : I don't know what is this field for
                end
            end
            
            %write new marker file
            %fname_mrk = fullfile(config{ipatient}.rawdir, MuseStruct_new{ipart}{idir}.directory,'Events.mrk');
			if isMicromed || isNeuralynx
				writeMuseMarkerfile(MuseStruct_new{ipart}{idir}, fname_mrk{ipart}{idir});
			elseif isBrainvision
				writeMuseMarkerfile_Brainvision(MuseStruct_new{ipart}{idir}, fname_mrk{ipart}{idir});
            end
        end
    end
end