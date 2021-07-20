function [MuseStruct]  = readMuseMarkers_discontinuousMicromed(cfg, force)
%In addition to readMuseMarkers.m, this scripts reads the text file
%associated with Micromed data and reconstruct real time if the EEG file
%was composed of several fragments, from different times. 

fname = fullfile(cfg.datasavedir, sprintf('%sMuseStruct.mat', cfg.prefix));

write = ft_getopt(cfg.muse, 'write', true);

if exist(fname,'file') && force == false
    fprintf('************************************\n');
    fprintf('** loading precomputed MuseStruct **\n');
    fprintf('************************************\n\n');
    load(fname,'MuseStruct');
    return
else
    fprintf('*************************\n');
    fprintf('** creating MuseStruct **\n');
    fprintf('*************************\n\n');
end

%get format to adapt script for each format
%specificities :
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% Go through different parts
for ipart = 1 : size(cfg.directorylist, 2)

    % Go through directory list
    for idir = 1 : size(cfg.directorylist{ipart}, 2)

        fprintf('Extracting artefact timings from  %s \n',cfg.directorylist{ipart}{idir});

        % find muse event file and read header if neaded
        if isNeuralynx
            temp        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir},'*.mrk'));
            name_mrk    = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
        elseif isMicromed
            name_mrk    = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.mrk']);
        elseif isBrainvision
            name_mrk    = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.vmrk']);
            load(fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '_header.mat']), 'header'); %header create by Paul's script during the conversion from Deltamed to Brainvision
        end

        if ~exist(name_mrk,'file')
            error('%s not found', name_mrk);
        else
            fprintf('Found markerfile: %s \n',name_mrk);
        end

        % read muse event file
        clear f markfile marks
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

        if isempty(markfile)
            fprintf('\n\n %s is empty!!! \n\n', name_mrk);
            continue
        end

        if isNeuralynx || isMicromed

            nmarkers        = str2double(markfile{find(strcmp('NUMBER OF MARKERS:', markfile)) + 1});
            classgroupid    = markfile(find(strcmp('CLASSGROUPID:', markfile)) + 1);
            name            = markfile(find(strcmp('NAME:', markfile)) + 1);
            comment         = markfile(find(strcmp('COMMENT:', markfile)) + 1);
            color           = markfile(find(strcmp('COLOR:', markfile)) + 1);
            editable        = markfile(find(strcmp('EDITABLE:', markfile)) + 1);
            classid         = markfile(find(strcmp('CLASSID:', markfile)) + 1);
            nrEvents        = str2double(markfile(find(strcmp('NUMBER OF SAMPLES:', markfile)) + 1));
            %correct bug paul :
            if length(name)-length(classgroupid) == 1
                classgroupid{end+1} = '+3';
            end


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
            clear name_temp name_mrk name classgroupid comment color editable classid

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
            hdr         = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp(1).name));
            f           = fopen(fullfile(temp(1).folder,[f,'.txt']));
            clear timestring
            
            if f >= 0
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
            
            else %error while loading katia text file header
                warning('Clock time not found, they will be ignored');
                MuseStruct{ipart}{idir}.starttime  = datetime.empty;
            end

        elseif isMicromed
            %read .txt text header. because time at
            %beginning is not available on the header
            %created by fieldtrip
            %also, the .txt file contains info about segment timings of the
            %data, in case data is not continuous

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
            MuseStruct{ipart}{idir}.starttime   = header.date;
            datafile                            = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
            hdr                                 = ft_read_header(datafile);

        end %end of recover real time. different for each format
        
        %%%%%%%%%
        
        %search if micromed data is continuous or segmented, and find where
        %data are cut
        iscontinuous = true;
        if isMicromed
            clear f segment* temp
            f = fopen([datafile(1:end-3) 'txt']);
            headerfile = {};
            while true
                l = fgetl(f);
                if ~ischar(l)
                    break
                end
                headerfile{end + 1} = l;
            end
            fclose(f);
            
            segment_line = find(strcmp('[SEGMENTS]', headerfile)) + 1; 
            segments = {};
            i=0;
            while ~isempty(headerfile{segment_line})
                i=i+1;
                segments(i) = headerfile(segment_line);
                segment_line = segment_line +1;
            end
            if ~isempty(segments)
                iscontinuous = false;
                MuseStruct{ipart}{idir}.iscontinuous = 'no';
                for isegment = 1:length(segments)
                    temp                         = strsplit(segments{isegment},'\t');
                    segment_datastart(isegment)  = str2double(erase(temp{6}, 'Sample: ')); %start sample in the cut data
                    segment_realstart(isegment)  = str2double(erase(temp{3}, 'Time: ')); %start sample in the raw data before cut
                end
                timediff =  segment_realstart - segment_datastart; %time to add to the cut data, to be synchronized with the raw non-cut data
            end
        end
        %%%%%%%%%%

        MuseStruct{ipart}{idir}.directory  = cfg.directorylist{ipart}{idir};
        
        if iscontinuous
            MuseStruct{ipart}{idir}.endtime    = MuseStruct{ipart}{idir}.starttime + seconds(hdr.nSamples / hdr.Fs - hdr.nSamplesPre / hdr.Fs);
        else
            segment_dataend = [diff(segment_datastart), (hdr.nSamples - hdr.nSamplesPre) - segment_datastart(end)]; %length of data segments 
            MuseStruct{ipart}{idir}.segments          = table; 
            MuseStruct{ipart}{idir}.segments.begclock = MuseStruct{ipart}{idir}.starttime + seconds(segment_realstart'./hdr.Fs); 
            MuseStruct{ipart}{idir}.segments.endclock = MuseStruct{ipart}{idir}.starttime + seconds((segment_realstart'+segment_dataend')./hdr.Fs);%make a draw to understand it
            MuseStruct{ipart}{idir}.endtime           = MuseStruct{ipart}{idir}.starttime + seconds((hdr.nSamples - hdr.nSamplesPre + timediff(end)) / hdr.Fs );
        end
        
        % create markers details in MuseStruct
        for imarker = 1 : nmarkers
            name{imarker} = strrep(name{imarker},'-','_'); % cant make fieldnames with minusses
            MuseStruct{ipart}{idir}.markers.(name{imarker}).events         = nrEvents(imarker);
            MuseStruct{ipart}{idir}.markers.(name{imarker}).comment        = comment{imarker};
            MuseStruct{ipart}{idir}.markers.(name{imarker}).color          = color{imarker};
            MuseStruct{ipart}{idir}.markers.(name{imarker}).editable       = editable{imarker};
            MuseStruct{ipart}{idir}.markers.(name{imarker}).classid        = classid{imarker};
            MuseStruct{ipart}{idir}.markers.(name{imarker}).classgroupid   = classgroupid{imarker};
            for ievent = 1 : nrEvents(imarker)
                MuseStruct{ipart}{idir}.markers.(name{imarker}).trialnum                 = marks{imarker}(ievent,1);
                MuseStruct{ipart}{idir}.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent,2);
                if isempty(MuseStruct{ipart}{idir}.starttime)
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent) = NaN;
                else
                    MuseStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent) = seconds(marks{imarker}(ievent,2)) + MuseStruct{ipart}{idir}.starttime;
                    if ~iscontinuous
                        timediff_idx_event = find(MuseStruct{ipart}{idir}.markers.(name{imarker}).synctime(ievent) >= segment_datastart./hdr.Fs, 1, 'last');
                        MuseStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent) = MuseStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent) + seconds(timediff(timediff_idx_event)/hdr.Fs);
                    end
                end
            end
        end
    end
    

    % check if data directory exists, if not create it
    if ~isfolder(cfg.datasavedir)
        ft_notice('creating directory %s', cfg.datasavedir);
        mkdir(cfg.datasavedir);
    end

end

if write
    save(fname,'MuseStruct');
end
