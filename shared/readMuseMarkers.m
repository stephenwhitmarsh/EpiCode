function [MuseStruct]  = readMuseMarkers(cfg, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readMuseMarkers
%
% Crawls and searches through patient- and recording-directories to
% extracts marker timings (time and samples) created by Muse. These can be
% used to create an overview of events (plotmarkers.m), segment data with
% FieldTrip (e.g. plotpeaks.m) and create artefact files for Spyking-Circus
% (writeSpykingCircusDeadFile.m). The resultant MuseStruct can be edited
% and written into a new markerfile for Muse with writeMuseMarkers.m
%
% The first argument is the configuration struct defined by a _setparams.m
% The second argument is whether or not to force (re)reading the marker
% file or to load from disk (i.e. as saved a previous time).
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
% with help from Jean-Didier Lemarechal & Craig Richter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = fullfile(cfg.datasavedir,sprintf('%sMuseStruct.mat',cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('************************************\n');
    fprintf('** loading precomputed MuseStruct **\n');
    fprintf('************************************\n\n');
    load(fname,'MuseStruct');
else
    
    if force == true
        fprintf('*******************************************\n');
        fprintf('** forced redoing of MuseStruct creation **\n');
        fprintf('*******************************************\n\n');
    else
        fprintf('*************************\n');
        fprintf('** creating MuseStruct **\n');
        fprintf('*************************\n\n');
    end
    
    %get format to adapt script for each format
    %specificities :
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            fprintf('Extracting artefact timings from  %s \n',cfg.directorylist{ipart}{idir});
            
            % find muse event file and read header if neaded
            if isNeuralynx
                temp        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.mrk'));
                name_mrk    = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
            elseif isMicromed
                name_mrk    = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.mrk']);
            elseif isBrainvision
                name_mrk    = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.vmrk']);
                load(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '_header']), 'header'); %header create by Paul's script during the conversion from Deltamed to Brainvision
            end
            
            
            if ~exist(name_mrk,'file')
                error('%s not found', name_mrk);
            else
                fprintf('Found markerfile: %s \n',name_mrk);
            end
            
            % read muse event file
            f = fopen(name_mrk, 'rt');
            markfile = {};
            while true
                l = fgetl(f); 
                if ~ischar(l)
                    break
                end
                markfile{end + 1} = l;
            end
            fclose(f);
            
            if ~isempty(markfile)
                            
            if isNeuralynx || isMicromed 

                nmarkers        = str2double(markfile{find(strcmp('NUMBER OF MARKERS:', markfile)) + 1});
                classgroupid    = markfile(find(strcmp('CLASSGROUPID:', markfile)) + 1);
                name            = markfile(find(strcmp('NAME:', markfile)) + 1);
                comment         = markfile(find(strcmp('COMMENT:', markfile)) + 1);
                color           = markfile(find(strcmp('COLOR:', markfile)) + 1);
                editable        = markfile(find(strcmp('EDITABLE:', markfile)) + 1);
                classid         = markfile(find(strcmp('CLASSID:', markfile)) + 1);
                nrEvents        = str2double(markfile(find(strcmp('NUMBER OF SAMPLES:', markfile)) + 1));
                
                for i = 1:length(nrEvents)
                    if nrEvents(i) > 0
                        fprintf('found %d occurances of %s \n', nrEvents(i), name{i});
                    end
                end
                
                % Get the events, time in seconds from onset of file
                j = find(strcmp('LIST OF SAMPLES:', markfile)) + 2;
                for imarker = 1 : nmarkers
                    marks{imarker} = str2num(char(markfile(j(imarker):j(imarker) + nrEvents(imarker) - 1)));
                    
                    % Convert from index origin 0 to 1
                    if nrEvents(imarker) ~= 0
                        marks{imarker}(:, 1) = marks{imarker}(:, 1) + 1;
                    end
                end
                
                
            elseif isBrainvision %other muse marker file format than neuralynx or micromed
               
                %Get marker names and infos
                name_temp            = markfile(find(strcmp('#TRIGGER COMMENTS', markfile)) + 1 : find(strcmp('Brain Vision Data Exchange Marker File, Version 1.0', markfile)) - 1);
                nmarkers             = length(name_temp);
                
                for imarker = 1 : nmarkers
                    name_temp{imarker}        = split(name_temp{imarker});
                    name{imarker}             = name_temp{imarker}{2};
                    classgroupid{imarker}     = [];
                    comment{imarker}          = [];
                    color{imarker}            = name_temp{imarker}{4}(8:end);
                    editable{imarker}         = [];
                    classid{imarker}          = [];
                end

                % Get the events, time in seconds from onset of file
                %index beginning of marker samples : 
                begin_mrk_list_index = find(strcmp('[Marker Infos]', markfile)) + 1;
                
                for i = 1:length(markfile)
                    markfile{i} = split(markfile{i},',');
                end
                
                for imarker = 1:nmarkers
                    sample_index = 0;
                    for isample = begin_mrk_list_index:length(markfile)
                        if strcmp(markfile{isample}{2},name{imarker})
                            sample_index = sample_index+1;
                            marks{imarker}(sample_index, 1) = 1;
                            marks{imarker}(sample_index, 2) = str2num(markfile{isample}{3})/header.Fs;
                        end
                    end
                    if sample_index == 0
                         marks{imarker} = [];
                    end
                end
                        
                
                for imarker = 1:nmarkers
                    nrEvents(imarker)        = size(marks{imarker},1);
                    if nrEvents(imarker)>0
                        fprintf('found %d occurances of %s \n', size(marks{imarker},1), name{imarker});
                    end
                end
                
                
                
            end
            
            
                
                                
                
                % recover "real time"
                
                if isNeuralynx
                    %from first Neurlynx .txt file
                    temp        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
                    [~, f, ~]   = fileparts(temp(1).name);
                    f           = fopen(fullfile(temp(1).folder,[f,'.txt']));
                    clear timestring
                    
                    % depending on amplifier, there are somewhat different
                    % formats of the txt file
                    ftype       = 'none';
                    while 1
                        tline = fgetl(f);
                        if ~ischar(tline), break, end
                        searchstring1 = '## Time Opened (m/d/y)';
                        searchstring2 = '-TimeCreated';
                        try
                            if length(tline) >= max(length(searchstring1))
                                if strcmp(tline(1:length(searchstring1)),searchstring1)
                                    timestring = tline;
                                    ftype = 1;
                                    disp('Great, found timestamp in header file - Type 1');
                                    break
                                end
                            end
                            if length(tline) >= max(length(searchstring2))
                                if strcmp(tline(1:length(searchstring2)),searchstring2)
                                    timestring = tline;
                                    ftype = 2;
                                    disp('Great, found timestamp in header file - Type 2');
                                    break
                                end
                            end
                        catch
                            disp('Warning: something weird happened reading the txt time');
                        end
                    end
                    fclose(f);
                    
                    % add real time of onset of file
                    timestring = strsplit(timestring);
                    switch ftype
                        case 1
                            headerdate = [cell2mat(timestring(5)) ' ' cell2mat(timestring(7))];
                            MuseStruct{ipart}{idir}.starttime  = datetime(headerdate,'Format','MM/dd/yy HH:mm:ss.SSS');
                            
                            
                        case 2
                            headerdate = [cell2mat(timestring(2)) ' ' cell2mat(timestring(3))];
                            MuseStruct{ipart}{idir}.starttime  = datetime(headerdate,'Format','yy/MM/dd HH:mm:ss.SSS');
                    end
                    
                elseif isMicromed
                    %read .bni text header. because time at
                    %beginning is not available on the header
                    %created by fieldtrip
                    
                    datafile = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
                    hdr = ft_read_header(datafile);
                    
                    f = fopen([datafile(1:end-3) 'bni']);
                    clear timestring
                    ftype = 'none';
                    search_date_time = 0;
                    while search_date_time<2
                        tline = fgetl(f);
                        if ~ischar(tline), break, end
                        searchstring1 = 'Date = ';
                        searchstring2 = 'Time = ';
                        try
                            if length(tline) >= max(length(searchstring1))
                                if strcmp(tline(1:length(searchstring1)),searchstring1)
                                    datestring = tline;
                                    ftype = 1;
                                    disp('Great, found datestamp in header file');
                                    search_date_time = search_date_time+1;
                                end
                            end
                            if length(tline) >= max(length(searchstring2))
                                if strcmp(tline(1:length(searchstring2)),searchstring2)
                                    timestring = tline;
                                    search_date_time = search_date_time+1;
                                    disp('Great, found timestamp in header file');
                                end
                            end
                        catch
                            disp('Warning: something weird happened reading the .bni time');
                        end
                    end
                    fclose(f);
                    
                    % add real time of onset of file
                    datestring = strsplit(datestring);
                    timestring = strsplit(timestring);
                    headerdate = [cell2mat(datestring(3)) ' ' cell2mat(timestring(3))];
                    MuseStruct{ipart}{idir}.starttime  = datetime(headerdate,'Format','MM/dd/yy HH:mm:ss');
                
                    
                elseif isBrainvision
                    MuseStruct{ipart}{idir}.starttime = header.date;
                    
                    
                end %end of recover real time. different for each format
                    
                
                MuseStruct{ipart}{idir}.directory  = cfg.directorylist{ipart}{idir};

                
                
                % create markers details in MuseStruct
                for imarker = 1 : nmarkers
                    name{imarker} = strrep(name{imarker},'-','_'); % cant make fieldnames with minusses
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).events         = []; 
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).comment        = comment{imarker};
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).color          = color{imarker};
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).editable       = editable{imarker};
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).classid        = classid{imarker};
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).classgroupid   = classgroupid{imarker};
                    for ievent = 1 : nrEvents(imarker)
                        MuseStruct{ipart}{idir}.markers.(name{imarker}).trialnum                 = marks{imarker}(ievent,1);
                        MuseStruct{ipart}{idir}.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent,2);
                        MuseStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent)            = seconds(marks{imarker}(ievent,2)) + MuseStruct{ipart}{idir}.starttime;                        
                    end
                end
            else
                fprintf('\n\n %s is empty!!! \n\n',name_mrk);
            end
        end
    end
    save(fname,'MuseStruct');
end