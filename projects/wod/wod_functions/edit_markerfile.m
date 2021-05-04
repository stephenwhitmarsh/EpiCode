

% writemusemarkerfiles do not work for brainvision data (because there is an other format of Muse marker file (.vmrk))

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'wod_functions']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparams_ouabaine;

ipart= 1;




pat_list = 1:size(config,2);
torename = {};%{'Baseline_Start', 'Analysis_Start'}; %{old1, new1; old2, new2; old3, new3};
toremove = {}; %{remove1, remove2, remove3};
toadd    = {'Seizure_start','Seizure_end'}; %{add1, add2, add3};


for ipatient = pat_list
    
    if isempty(config{ipatient})
        continue
    end
    
    
  %find concatenated LFP (see wod_concatenateLFP.m)
        config{ipatient}.rawdir                = config{ipatient}.concatdata_path;
        config{ipatient}.directorylist{ipart}  = {config{ipatient}.prefix};
        
    
    % backup markerfiles
    for ipart = 1 : size(config{ipatient}.directorylist,2)
        for idir = 1 : size(config{ipatient}.directorylist{ipart},2)
            
            fname_mrk_temp    = fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir},'Events*.mrk');
            temp              = dir(fname_mrk_temp);
            if size(temp,1) == 1
                fname_mrk{ipart}{idir}  = fullfile(temp.folder, temp.name);
            else
                %if no 'event.mrk' file (neuralynx), search for micromed
                %.mrk or brainvision . vmrk file, which has the same name
                %as the data file :
                fname_mrk_temp    = fullfile(config{ipatient}.rawdir,[config{ipatient}.directorylist{ipart}{idir},'.*mrk']);
                temp              = dir(fname_mrk_temp);
                if size(temp,1) == 1
                    fname_mrk{ipart}{idir}  = fullfile(temp.folder, temp.name);
                else
                    error('error when searching marker file %s', fname_mrk_temp);
                end
            end
                
            % backup markerfile
            if ~exist(config{ipatient}.muse.backupdir,'DIR')
                error('Backup directory does not exist');
            end
            [~, d] = fileparts(config{ipatient}.directorylist{ipart}{idir});
            dir_backup = fullfile(config{ipatient}.muse.backupdir,config{ipatient}.prefix(1:end-1), d);
            
            if ~exist(dir_backup, 'DIR')
                fprintf('Creating directory: %s\n', dir_backup);
                %             eval(sprintf('!mkdir %s', fullfile(config{ipatient}.muse.backupdir, d)));
                mkdir(dir_backup);
            end
            fname_backup = sprintf('Events_%s.mrk', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
            %         eval(sprintf('!cp %s %s', fname_mrk, fullfile(config{ipatient}.muse.backupdir, d, fname_backup)));
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
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).comment    = 'created with editMarkerfiles.m';
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).color      = '#00ff00'; %default Muse color
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).editable   = 'Yes';
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).classid    = sprintf('+%d',numel(fieldnames(MuseStruct_new{ipart}{idir}.markers)));
                    MuseStruct_new{ipart}{idir}.markers.(toadd{imarker}).classgroupid = '+3'; %I don't know what is this field for
                end
            end
            
            %write new marker file
            %fname_mrk = fullfile(config{ipatient}.rawdir, MuseStruct_new{ipart}{idir}.directory,'Events.mrk');
            writeMuseMarkerfile(MuseStruct_new{ipart}{idir}, fname_mrk{ipart}{idir});
            
        end
    end
end