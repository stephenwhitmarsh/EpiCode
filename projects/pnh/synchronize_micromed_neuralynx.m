function [micromed_markers]  = synchronize_micromed_neuralynx(cfg)

fname_output = fullfile(cfg.datasavedir,'markers_sync_stim.mat');
if exist(fname_output,'file') && cfg.force == false
    load(fname_output,'micromed_markers');
else

    feature('DefaultCharacterSet', 'UTF8') %# for all Character support, or 'CP1252'
    template_width      = 27; % officially 25ms, so redo according to samplerate
    template_pad        = 15;
    
    % prepare MicroMed
    clear micromed*
    micromed_path           = fullfile(cfg.patientdir,'eegmicromed');
    micromed_datafilelist   = dir2(fullfile(micromed_path,'*.TRC'));
    micromed_txtfiles       = dir2(fullfile(micromed_path,'*.txt'));
    micromed_infofiles      = dir2(fullfile(micromed_path,'*.bni'));
    
    for ifile = 1 : size(micromed_datafilelist,1)
        micromed_datafile{ifile}                = micromed_datafilelist(ifile);
        micromed_hdr{ifile}                     = ft_read_header(fullfile(micromed_datafile{ifile}.folder,micromed_datafile{ifile}.name)); % read header from first ncs file
        fname                                   = fullfile(micromed_path,micromed_txtfiles(ifile).name);
        opts                                    = detectImportOptions(fname);
        opts.Delimiter                          = {'\t', ' '};
        opts.CommentStyle                       = '[';
        opts.VariableNames                      = {'Date','Time','x1','Sample','x2','Label','Current','Frequency','Duration','x3'};
        opts.VariableTypes                      = {'char','char','char','char','char','char','char','char','char','char'};
        T                                       = readtable(fname, opts,'FileEncoding','CP1252');
        %         T = removevars(T,{'x1','x2','x3'}); % only MATLAB 2018
        T(strcmp(T.Current,'Interrogatoire'),:) = [];
        T.Current                               = erase(T.Current,'mA');
        T.Current                               = str2double(T.Current);
        T.Frequency                             = erase(T.Frequency,'Hz');
        T.Frequency                             = str2double(T.Frequency);
        T.Duration                              = erase(T.Duration,'µsec');
        T.Duration                              = str2double(T.Duration);
        T.Sample                                = str2double(T.Sample);
        T.DataFileNr                            = repmat(ifile,height(T),1);
        T.DataFileName                          = repmat(micromed_datafilelist(ifile).name,height(T),1);
        T.DateFileFolder                        = repmat(micromed_datafilelist(ifile).folder,height(T),1);
        T.DateTime                              = datetime([cell2mat(T.Date), cell2mat(T.Time)],'Format','yyy/MM/ddHH:mm:ss.SSS');
        T.DateTime                              = datetime(T.DateTime,'Format','yyy/MM/dd HH:mm:ss.SSS');
        T.Date                                  = datetime(cell2mat(T.Date), 'Format','yyy/MM/dd');
        T.Time                                  = datetime(cell2mat(T.Time),'Format','HH:mm:ss.SSS');
        
        % remove all special characters
        T.Label                                 = regexprep(T.Label,'[^a-zA-Z0-9-]','');
        
        % create labels consistent with channel names
%         temp = cell2mat(T.Label);
        for i = 1 : height(T)
            s = strsplit(T.Label{i,:},'-');
            T.contact1{i} = [s{1}(1:end-1) '_' s{1}(end)];
            T.contact2{i} = [s{2}(1:end-1) '_' s{2}(end)];
        end
        micromed_markers_file{ifile}            = T;
        
        fname                                   = fullfile(micromed_infofiles(ifile).folder,micromed_infofiles(ifile).name);
        opts                                    = detectImportOptions(fname,'FileType','text','Delimiter','=');
        opts.VariableNames                      = {'key','value'};
        opts.Delimiter                          = {'='};
        opts.VariableTypes                      = {'char','char'};
        T2                                      = readtable(fname,opts);
        micromed_date{ifile}                    = datetime(cell2mat(T2.value(find(strcmp(T2.key,'Date')))),'Format','MM/dd/yyyy');
        micromed_time{ifile}                    = datetime(cell2mat(T2.value(find(strcmp(T2.key,'Time')))),'Format','HH/mm/ss');
        micromed_datetime{ifile}                = datetime([cell2mat(T2.value(find(strcmp(T2.key,'Date')))), cell2mat(T2.value(find(strcmp(T2.key,'Time'))))],'Format','MM/dd/yyyyHH:mm:ss');
        micromed_datetime{ifile}                = datetime(micromed_datetime{ifile},'Format','yyy/MM/dd HH:mm:ss');
        
        % load sync channel from MicroMed
        hdr_MicroMed{ifile}                     = ft_read_header(fullfile(micromed_datafile{ifile}.folder,micromed_datafile{ifile}.name));
        sync_MicroMed{ifile}                    = ft_read_data(fullfile(micromed_datafile{ifile}.folder,micromed_datafile{ifile}.name),'chanindx',find(strcmp(hdr_MicroMed{ifile}.label,'SYNC1')));
        
    end
    
    % prepare NeuraLynx
    clear neuralynx*
    neuralynx_dirlist_temp = dir2(fullfile(cfg.patientdir,'eeg'));
    for idir = 1 : size(neuralynx_dirlist_temp,1)
        neuralynx_datafiles{idir}   = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,'*m*.ncs'));
        neuralynx_txtfiles{idir}    = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,'*.txt'));
        neuralynx_hdr{idir}         = ft_read_header(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,neuralynx_datafiles{idir}(1).name)); % read header from first ncs file
        neuralynx_syncfiles{idir}   = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,[cfg.timecode,'.ncs']));
        neuralynx_synchdr{idir}     = ft_read_header(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,neuralynx_syncfiles{idir}(1).name));
        
        f = fopen(fullfile(neuralynx_txtfiles{idir}(1).folder,neuralynx_txtfiles{idir}(1).name));
        clear timestring
        while 1
            tline = fgetl(f);
            if ~ischar(tline), break, end
            searchstring = '## Time Opened (m/d/y)';
            try
                if strcmp(tline(1:length(searchstring)),searchstring)
                    neuralynx_timestring = tline;
                    disp('Great, found timestamp in header file');
                    break
                end
            catch
                disp('Warning: something weird happened reading the txt time');
            end
        end
        fclose(f);
        
        % add real time of onset of neuralynx file
        neuralynx_timestring         = strsplit(neuralynx_timestring);
        neuralynx_headerdate{idir}   = [cell2mat(neuralynx_timestring(5)) ' ' cell2mat(neuralynx_timestring(7))];
        neuralynx_starttime{idir}    = datetime(neuralynx_headerdate{idir},'Format','MM/dd/yy HH:mm:ss.SSS');
        neuralynx_endtime{idir}      = neuralynx_starttime{idir} + seconds(neuralynx_hdr{idir}.nSamples /neuralynx_hdr{idir}.Fs);
        
        % load sync channel from NeuraLynx
        sync_NeuraLynx{idir}        = ft_read_data(fullfile(neuralynx_syncfiles{idir}.folder,neuralynx_syncfiles{idir}.name));
        sync_NeuraLynx{idir}        = downsample(sync_NeuraLynx{idir},neuralynx_synchdr{idir}.Fs / micromed_hdr{ifile}.Fs);
        
    end
    
    for iMicroMed = 1 : size(micromed_markers_file,2)
        
        % debug recording in which timecode was only recorded midway
        if strcmp(micromed_datafilelist(iMicroMed).name,'EEG_338059.TRC')
            startindx = 529500;
        else
            startindx = 1;
        end

        % Extract MicroMed timecode
        [~, locs_MM]        = findpeaks(cfg.timecodeoriMM*diff(sync_MicroMed{iMicroMed}),'MinPeakProminence',50);
        locs_MM             = locs_MM(locs_MM > startindx);
        template_MM         = sync_MicroMed{iMicroMed}(locs_MM(1)-template_pad:locs_MM(1)+template_width+template_pad);
        template_MM_lags    = locs_MM(1)-template_pad:locs_MM(1)+template_width+template_pad;
        
        %                 [istart_out,istop_out,dist_out] = findsignal(sync_MicroMed,template_MM,'MaxDistance',+inf);
        %                 figure
        %                 subplot(4,1,2);
        %                 hist(dist_out,10000)
        
        [istart_out_MM,istop_out_MM,dist_out_MM] = findsignal(sync_MicroMed{iMicroMed},template_MM,'MaxDistance',2e5);
        istart_out_MM = istart_out_MM(istart_out_MM > startindx);
        markers_MM = zeros(1,length(sync_MicroMed));
        markers_MM(repmat(istart_out_MM+template_pad,template_width,1)+(1:template_width)') = 1;
        
        % plot overlap with timecourse
        fig = figure;
        subplot(2,1,1);
        plot(template_MM_lags,template_MM);
        subplot(2,1,2); hold;
        plot(sync_MicroMed{iMicroMed});
%                 plot(sync_MicroMed{iMicroMed}(529500:end));

        plot(markers_MM*max(abs(sync_MicroMed{iMicroMed})));
        
        % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,['timecode_template_',micromed_datafile{iMicroMed}.name(1:end-4),'.pdf']),'-r600');
        
        % sort timecodes
        istart_out_MM = sort(istart_out_MM);
        %   figure; hist(diff(istart_out_MM),1000)
        
        % those that are close together identifies them as timecodes
        istimestamp = diff(istart_out_MM) < 400;
        isfirsttimestamp = diff(istimestamp);
        isfirsttimestamp(isfirsttimestamp == -1) = 0;
        istart_out_MM = istart_out_MM(find(isfirsttimestamp)+1); % samplerate / 10
        
        % Search all neuralynx files
        for iNeuraLynx = 1 : size(neuralynx_datafiles,2)
            
            % look for any marker in the MicroMed marker file that corresponds
            % with the approximate time recorded in the NeuraLynx datafile
            hasmarker = 0;
            for iMarker = 1 : size(micromed_markers_file{iMicroMed},1)
                if isbetween(micromed_markers_file{iMicroMed}.DateTime(iMarker),neuralynx_starttime{iNeuraLynx},neuralynx_endtime{iNeuraLynx})
                    hasmarker = hasmarker + 1;
                end
            end
            
            if hasmarker > 0
                fprintf('Found %d markers in %s within datetime of %s! \nLoading data \n',hasmarker,micromed_markers_file{iMicroMed}.DataFileName(1,:), neuralynx_datafiles{idir}(1).name);
                
                % neuralynx
                
                [pks_NL, locs_NL]   = findpeaks(cfg.timecodeoriNL*diff(sync_NeuraLynx{iNeuraLynx}),'MinPeakProminence',50);
                template_NL         = sync_NeuraLynx{iNeuraLynx}(locs_NL(1)-template_pad:locs_NL(1)+template_width+template_pad);
                %                 [istart_out,istop_out,dist_out] = findsignal(sync_NeuraLynx,template_NL,'MaxDistance',+inf);
                %                 subplot(3,1,2);
                %                 hist(dist_out,10000)
                [istart_out_NL,istop_out_NL,dist_out_NL]    = findsignal(sync_NeuraLynx{iNeuraLynx},template_NL,'MaxDistance',2e5);
                markers_NL                                  = zeros(1,length(sync_NeuraLynx{iNeuraLynx}));
                markers_NL(repmat(istart_out_NL+template_pad,template_width,1)+(1:template_width)') = 1;
                
                figure;
                subplot(2,1,1);
                plot(template_NL);
                subplot(2,1,2); hold;
                plot(sync_NeuraLynx{iNeuraLynx}(1:200000));
                plot(markers_NL(1:200000)*max(abs(sync_NeuraLynx{iNeuraLynx}(1:200000))));
                
                clear istart_out_map istop_out_map dist_out_map dist_time_map istart_out_MM_clean
                icounter = 1;
                for itimecode = 1:length(istart_out_MM)
                    try
                    timecodetemplate = markers_MM(istart_out_MM(itimecode):istart_out_MM(itimecode)+450); % samplerate / 10 * 1.
                    [istart_out_map_all,istop_out_map_all,dist_out_map_all] = findsignal(markers_NL,timecodetemplate,'MaxNumSegments',4, 'MaxDistance',0);
                    if ~isempty(istart_out_map_all)
                        fprintf('%d/%d: Found %d marker(s) mapping between NeuraLynx and MicroMed, %d selected in total\n', itimecode, length(istart_out_MM), length(istart_out_map_all),icounter);
                        timediff = (neuralynx_starttime{iNeuraLynx} + seconds(istart_out_map_all / hdr_MicroMed{iMicroMed}.Fs)) - (micromed_datetime{iMicroMed} + seconds(istart_out_MM(itimecode) / hdr_MicroMed{iMicroMed}.Fs)) ; % use resample frequencies
                        [~, I]                          = min(abs(timediff));
                        istart_out_map(icounter)        = istart_out_map_all(I);
                        istart_out_MM_clean(icounter)   = istart_out_MM(itimecode);
                        istop_out_map(icounter)         = istop_out_map_all(I);
                        dist_out_map(icounter)          = dist_out_map_all(I);
                        dist_time_map(icounter)         = timediff(I);
                        icounter = icounter + 1;
                    else
                        fprintf('%d/%d: Found no marker(s) between NeuraLynx and MicroMed, %d selected in total\n', itimecode, length(istart_out_MM),icounter);
                    end
                    catch
                        fprintf('failed\n');
                    end
                end
                
                timethreshold       = abs(dist_time_map) < minutes(10);
%                 timethreshold       = abs(dist_time_map) < minutes(70);
    
                dist_time_map       = dist_time_map(timethreshold);
                istart_out_map      = istart_out_map(timethreshold);
                istart_out_MM_clean = istart_out_MM_clean(timethreshold);
                istop_out_map       = istop_out_map(timethreshold);
                dist_out_map        = dist_out_map(timethreshold);
                
                clear timecodetemplate_raw timecodetemplate timecodetemplate_found_raw timecodetemplate_found
                for itimecode = 1:length(istart_out_MM_clean)
                    timecodetemplate_raw(itimecode,:)         = sync_MicroMed{iMicroMed}(istart_out_MM_clean(itimecode):istart_out_MM_clean(itimecode)+450);
                    timecodetemplate(itimecode,:)             = markers_MM(istart_out_MM_clean(itimecode):istart_out_MM_clean(itimecode)+450);
                    timecodetemplate_found_raw(itimecode,:)   = sync_NeuraLynx{iNeuraLynx}(istart_out_map(itimecode):istart_out_map(itimecode)+450);
                    timecodetemplate_found(itimecode,:)       = markers_NL(istart_out_map(itimecode):istart_out_map(itimecode)+450);
                end
                
                fig = figure;
                subplot(1,4,1); imagesc(timecodetemplate_raw);         set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight; title('MicroMed raw signal');
                subplot(1,4,2); imagesc(timecodetemplate);             set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight; title('MicroMed binarized');
                subplot(1,4,3); imagesc(timecodetemplate_found_raw);   set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight; title('Neuralynx raw signal');
                subplot(1,4,4); imagesc(timecodetemplate_found);       set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight; title('Neuralynx binarized');
                
                % print to file
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,['timecode_binarized_',micromed_datafile{iMicroMed}.name(1:end-4),'-',neuralynx_datafiles{iNeuraLynx}(1).name(1:end-4),'.pdf']),'-r300');
                
                fig = figure;
                subplot(4,1,1); plot(istart_out_MM_clean); title('MicroMed samples of timecodes');
                subplot(4,1,2); plot(istart_out_map); title('Mapped to Neuralynx samples of timecodes');
                subplot(4,1,3); plot(istart_out_map-istart_out_MM_clean); title('Distance in samples');
                subplot(4,1,4); plot(dist_time_map) ; title('Approx. distance in time');
                
                % print to file
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,['timecode_mapping_',micromed_datafile{iMicroMed}.name(1:end-4),'-',neuralynx_datafiles{iNeuraLynx}(1).name(1:end-4),'.pdf']),'-r300');
                
                sampleshift = median(istart_out_map - istart_out_MM_clean); % calculate resampleing from the beginning
                samples_NL  = (micromed_markers_file{iMicroMed}.Sample + sampleshift) / micromed_hdr{iMicroMed}.Fs * neuralynx_hdr{idir}.Fs;
                seconds_NL  = (micromed_markers_file{iMicroMed}.Sample + sampleshift) / micromed_hdr{iMicroMed}.Fs;
               
                inNeuraLynxFile = samples_NL>0 & samples_NL < neuralynx_hdr{iNeuraLynx}.nSamples; % sort out with resample freq
                micromed_markers_file{iMicroMed}.Sample_NL(inNeuraLynxFile)             = samples_NL(inNeuraLynxFile);
                micromed_markers_file{iMicroMed}.Seconds_NL(inNeuraLynxFile)            = seconds_NL(inNeuraLynxFile);
                micromed_markers_file{iMicroMed}.samplehift(inNeuraLynxFile)           = repmat(sampleshift,sum(inNeuraLynxFile),1); 
                
                temp = repmat(neuralynx_datafiles{iNeuraLynx}(1).folder,height(micromed_markers_file{iMicroMed}),1);
                micromed_markers_file{iMicroMed}.DataFileFolder_NL(inNeuraLynxFile,:)   = temp(inNeuraLynxFile,:);     
            else
                fprintf('Found no marker in %s within datetime of %s \n',micromed_markers_file{iMicroMed}.DataFileName(1,:), neuralynx_datafiles{idir}(1).name);
            end
        end
    end 
    micromed_markers = MultiOuterJoin(micromed_markers_file{:});
    save(fname_output,'micromed_markers');
end
disp('done');
