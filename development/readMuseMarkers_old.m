function [MuseStruct_micro, MuseStruct_macro]  = readMuseMarkers(cfg,force)

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
% Code by Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
% with help from Jean-Didier Lemarechal & Craig Richter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all');

%% Micro electrodes, after repeat same for macro electrodes


for type = ["micro","macro"]
    
    fname = fullfile(cfg.datasavedir,sprintf('%sMuseStruct_%s.mat',cfg.prefix,type));
    
    if exist(fname,'file') && force == false
        fprintf('******************************************\n');
        fprintf('** loading precomputed MuseStruct %s **\n',type);
        fprintf('******************************************\n\n');
        load(fname,'MuseStruct');
    else
        
        if force == true
            fprintf('*************************************************\n');
            fprintf('** forced redoing of MuseStruct creation %s **\n',type);
            fprintf('*************************************************\n\n');
        else
            fprintf('*******************************\n');
            fprintf('** creating MuseStruct %s **\n',type);
            fprintf('*******************************\n\n');    
        end
        
        % Go through every data-directory in patient-directory
        ldir = dir(fullfile(cfg.rawdir,cfg.directory_searchstring));
        
        for idir = 1 : size(ldir,1)
            
            fprintf('Extracting artefact timings from  %s \n',ldir(idir).name);
            
            % check if directory has all requested files
            hasfiles = 0;
            channr = 1;
            for ichan = 1 : size(cfg.labels.(type),2)
                if ldir(idir).isdir
                    d = dir2(fullfile(cfg.rawdir,ldir(idir).name,['*',cfg.labels.(type){ichan},'*.ncs']));
                    if ~isempty(d)
                        
                        % add some bookkeeping for later
                        MuseStruct{idir}.directory = fullfile(ldir(idir).folder,ldir(idir).name);
                        MuseStruct{idir}.filenames{channr} = d.name;
                        channr = channr + 1;
                        hasfiles = 1;
                    else
                        fprintf('*** WARNING: cannot find channel %s data in %s\n',cfg.labels.(type){ichan},fullfile(cfg.rawdir,ldir(idir).name));
                    end
                end
            end
            
            if hasfiles
                
                fprintf('Extracting artefact timings based on header of %s \n',MuseStruct{idir}.filenames{1});
                temp = dir(fullfile(cfg.rawdir,ldir(idir).name,'*.mrk'));
                name_mrk = fullfile(cfg.rawdir,ldir(idir).name,temp.name);
                
                if ~exist(name_mrk)
                    error('%s not found', name_mrk);
                end
                
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
                    nmarkers        = str2num(markfile{strmatch('NUMBER OF MARKERS:', markfile, 'exact') + 1});
                    classgroupid    = markfile(strmatch('CLASSGROUPID:', markfile, 'exact') + 1);
                    name            = markfile(strmatch('NAME:', markfile, 'exact') + 1);
                    comment         = markfile(strmatch('COMMENT:', markfile, 'exact') + 1);
                    color           = markfile(strmatch('COLOR:', markfile, 'exact') + 1);
                    editable        = markfile(strmatch('EDITABLE:', markfile, 'exact') + 1);
                    classid         = markfile(strmatch('CLASSID:', markfile, 'exact') + 1);
                    nrEvents        = str2num(char(markfile(strmatch('NUMBER OF SAMPLES:', markfile, 'exact') + 1)));
                    
                    for i = 1:length(nrEvents)
                        if nrEvents(i) > 0
                            fprintf('found %d occurances of %s \n', nrEvents(i), name{i});
                        end
                    end
                    
                    % Get the events, time in seconds from onset of file
                    j = strmatch('LIST OF SAMPLES:', markfile, 'exact') + 2;
                    for i = 1 : nmarkers
                        marks{i} = str2num(char(markfile(j(i):j(i) + nrEvents(i) - 1)));
                        
                        % Convert from index origin 0 to 1
                        if nrEvents(i) ~= 0
                            marks{i}(:, 1) = marks{i}(:, 1) + 1;
                        end
                    end
                    
                    % take the Neurlynx data headerinfo from the first .txt file to
                    % recover real time
                    datafile = fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames{1});
                    hdr = ft_read_header(datafile);
                    
                    f = fopen([datafile(1:end-3) 'txt']);
                    clear timestring
                    ftype = 'none';
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
                            MuseStruct{idir}.starttime  = datetime(headerdate,'Format','MM/dd/yy HH:mm:ss.SSS');
                            
                        case 2
                            headerdate = [cell2mat(timestring(2)) ' ' cell2mat(timestring(3))];
                            MuseStruct{idir}.starttime  = datetime(headerdate,'Format','yy/MM/dd HH:mm:ss.SSS');    
                    end        
                    
                    % convert time to samples - FROM ONSET OF FILE
                    for imarker = 1 : nmarkers
                        name{imarker} = strrep(name{imarker},'-','_'); % cant make fieldnames with minusses
                        
                       
                        MuseStruct{idir}.markers.(name{imarker}).events         = [];
                        MuseStruct{idir}.markers.(name{imarker}).comment        = comment{imarker};
                        MuseStruct{idir}.markers.(name{imarker}).color          = color{imarker};
                        MuseStruct{idir}.markers.(name{imarker}).editable       = editable{imarker};
                        MuseStruct{idir}.markers.(name{imarker}).classid        = classid{imarker};
                        MuseStruct{idir}.markers.(name{imarker}).classgroupid   = classgroupid{imarker};
                        
                        for ievent = 1 : nrEvents(imarker)
                            
                            % determine the location of the marker, expressed in samples
                            MuseStruct{idir}.markers.(name{imarker}).trialnum                 = marks{imarker}(ievent,1);
                            MuseStruct{idir}.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent,2);
                            MuseStruct{idir}.markers.(name{imarker}).begsample                = (MuseStruct{idir}.markers.(name{imarker}).trialnum-1)*hdr.nSamples + 1;    % of the trial, relative to the start of the datafile
                            MuseStruct{idir}.markers.(name{imarker}).endsample                = (MuseStruct{idir}.markers.(name{imarker}).trialnum  )*hdr.nSamples;        % of the trial, relative to the start of the datafile
                            MuseStruct{idir}.markers.(name{imarker}).offset                   = round(MuseStruct{idir}.markers.(name{imarker}).synctime*hdr.Fs);           % this is the offset (in samples) relative to time t=0 for this trial
                            MuseStruct{idir}.markers.(name{imarker}).offset                   = MuseStruct{idir}.markers.(name{imarker}).offset + hdr.nSamplesPre;         % and time t=0 corrsponds with the nSamplesPre'th sample
                            MuseStruct{idir}.markers.(name{imarker}).clock(ievent)            = seconds(marks{imarker}(ievent,2)) + MuseStruct{idir}.starttime;
                            %                     MuseStruct{idir}.markers.(name{imarker}).syncsample(ievent)       = MuseStruct{idir}.markers.(name{imarker}).begsample + MuseStruct{idir}.markers.(name{imarker}).offset;
                            
                            % event structure
                            MuseStruct{idir}.markers.(name{imarker}).events(end+1).value      = [];
                            MuseStruct{idir}.markers.(name{imarker}).events(end ).sample      = MuseStruct{idir}.markers.(name{imarker}).begsample + MuseStruct{idir}.markers.(name{imarker}).offset(ievent);
                            MuseStruct{idir}.markers.(name{imarker}).events(end ).duration    = 0;
                            MuseStruct{idir}.markers.(name{imarker}).events(end ).offset      = MuseStruct{idir}.markers.(name{imarker}).offset(ievent);
                            MuseStruct{idir}.markers.(name{imarker}).events(end ).time        = marks{imarker}(ievent,2);
                        end
                    end
                else
                    fprintf('\n\n %s is empty!!! \n\n',name_mrk);
                end
            end
        end
        save(fname,'MuseStruct');
    end
    eval(sprintf('MuseStruct_%s = MuseStruct;',type));
end % type
