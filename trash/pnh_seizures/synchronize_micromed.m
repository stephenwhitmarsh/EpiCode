function [micromed_markers]  = synchronize_micromed

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/

patientdir = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/';
patientdir = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/';


ft_defaults
% feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
feature('DefaultCharacterSet', 'UTF8') %# for all Character support


clear neuralynx*
neuralynx_dirlist_temp = dir2(fullfile(patientdir,'eeg'));
for idir = 1 : size(neuralynx_dirlist_temp,1)
    neurolynx_datafiles{idir}   = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,'*m*.ncs'));
    neurolynx_txtfiles{idir}    = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,'*.txt'));
    neurolynx_hdr{idir}         = ft_read_header(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,neurolynx_datafiles{idir}(1).name)); % read header from first ncs file
    neurolynx_syncfiles{idir}   = dir2(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,'*SYNC*.ncs'));
    neurolynx_synchdr{idir}     = ft_read_header(fullfile(neuralynx_dirlist_temp(idir).folder,neuralynx_dirlist_temp(idir).name,neurolynx_syncfiles{idir}(1).name));
end

clear micromed*
micromed_path           = fullfile(patientdir,'eegmicromed');
micromed_datafilelist   = dir2(fullfile(micromed_path,'*.TRC'));
micromed_txtfiles       = dir2(fullfile(micromed_path,'*.txt'));
micromed_infofiles      = dir2(fullfile(micromed_path,'*.bni'));

for ifile = 1 : size(micromed_datafilelist,1)
        micromed_datafile{ifile}                    = micromed_datafilelist(ifile);
        micromed_hdr{ifile}                         = ft_read_header(fullfile(micromed_datafile{ifile}.folder,micromed_datafile{ifile}.name)); % read header from first ncs file
        fname                                   = fullfile(micromed_path,micromed_txtfiles(ifile).name);
        opts                                    = detectImportOptions(fname);
        opts.Delimiter                          = {'\t', ' '};
        opts.CommentStyle                       = '[';
        opts.VariableNames                      = {'Date','Time','x1','Sample','x2','Label','Current','Frequency','Duration','x3'};
        opts.VariableTypes                      = {'char','char','char','char','char','char','char','char','char','char'};
        T                                       = readtable(fname, opts,'FileEncoding','CP1252');
        %         T                                       = removevars(T,{'x1','x2','x3'});
        T(strcmp(T.Current,'Interrogatoire'),:) = [];
        T.Label                                 = erase(T.Label,'ù');
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

end



                
%% serach for each neuralynx directory which markers occur in the micromed


% take the Neurlynx data headerinfo from the first .txt file to
% recover real time


for idir = 1 : size(neurolynx_datafiles,2)
    
    f = fopen(fullfile(neurolynx_txtfiles{idir}(1).folder,neurolynx_txtfiles{idir}(1).name));
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
    neuralynx_timestring                  = strsplit(neuralynx_timestring);
    neuralynx_headerdate                  = [cell2mat(neuralynx_timestring(5)) ' ' cell2mat(neuralynx_timestring(7))];
    neuralynx_starttime                   = datetime(neuralynx_headerdate,'Format','MM/dd/yy HH:mm:ss.SSS');
    neuralynx_endtime                     = neuralynx_starttime + seconds(neurolynx_hdr{idir}.nSamples /neurolynx_hdr{idir}.Fs);
    
    % loop over events from MicroMed and see if any fall within the
    % NeuraLynx time
    
    hasmarker = 0;
    for iMicroMed = 1 : size(micromed_markers_file,2)
        for iMarker = 1 : size(micromed_markers_file{iMicroMed},1)
            if isbetween(micromed_markers_file{iMicroMed}.DateTime(iMarker),neuralynx_starttime,neuralynx_endtime)
                hasmarker = 1;
            end
        end
        
        if hasmarker == 1
            fprintf('Found marker in %s within datetime of %s. Loading data \n',micromed_markers_file{iMicroMed}.DataFileName(1,:), neurolynx_datafiles{idir}(1).name);
            % now we have found a MicroMed file with relevant marker
            % for current NeuraLynx file!
            
            
            % load sync channel from NeuraLynx
            sync_NeuraLynx      = ft_read_data(fullfile(neurolynx_syncfiles{idir}.folder,neurolynx_syncfiles{idir}.name));
            sync_NeuraLynx      = downsample(sync_NeuraLynx,4);
            
            % load sync channel from MicroMed
            hdr_MicroMed        = ft_read_header(fullfile(micromed_datafile{iMicroMed}.folder,micromed_datafile{iMicroMed}.name));
            sync_MicroMed       = ft_read_data(fullfile(micromed_datafile{iMicroMed}.folder,micromed_datafile{iMicroMed}.name),'chanindx',find(strcmp(hdr_MicroMed.label,'SYNC1')));
            
            template_width      = 27;
            template_pad        = 15;
            
            % neuralynx
            [pks_NL, locs_NL]   = findpeaks(-diff(sync_NeuraLynx),'MinPeakProminence',50);
            template_NL         = sync_NeuraLynx(locs_NL(1)-template_pad:locs_NL(1)+template_width+template_pad);
            %                 [istart_out,istop_out,dist_out] = findsignal(sync_NeuraLynx,template_NL,'MaxDistance',+inf);
            %                 subplot(3,1,2);
            %                 hist(dist_out,10000)
            [istart_out_NL,istop_out_NL,dist_out_NL]    = findsignal(sync_NeuraLynx,template_NL,'MaxDistance',2e5);
            markers_NL                                  = zeros(1,length(sync_NeuraLynx));
            markers_NL(repmat(istart_out_NL+template_pad,template_width,1)+(1:template_width)') = 1;
            
            figure;
            subplot(2,1,1);
            plot(template_NL);
            
            subplot(2,1,2); hold;
            plot(sync_NeuraLynx);
            plot(markers_NL*max(abs(sync_NeuraLynx)));
                
            % MicroMed timecode
            [pks_MM, locs_MM]   = findpeaks(-diff(sync_MicroMed),'MinPeakProminence',50);
            template_MM         = sync_MicroMed(locs_MM(1)-template_pad:locs_MM(1)+template_width+template_pad);
            template_MM_lags    = locs_MM(1)-template_pad:locs_MM(1)+template_width+template_pad;
            
            %                 [istart_out,istop_out,dist_out] = findsignal(sync_MicroMed,template_MM,'MaxDistance',+inf);
            %                 figure
            %                 subplot(4,1,2);
            %                 hist(dist_out,10000)
            
            [istart_out_MM,istop_out_MM,dist_out_MM] = findsignal(sync_MicroMed,template_MM,'MaxDistance',2e5);
            markers_MM = zeros(1,length(sync_MicroMed));
            markers_MM(repmat(istart_out_MM+template_pad,template_width,1)+(1:template_width)') = 1;
            
            figure;
            subplot(3,1,1);
            plot(template_MM_lags,template_MM);
            
            subplot(3,1,2); hold;
            plot(sync_MicroMed(1:200000));
            plot(markers_MM(1:200000));
            
            % sort timecodes
            istart_out_MM = sort(istart_out_MM);
            figure; hist(diff(istart_out_MM),1000)
            
            % select only those that are close together, so are a timecode
            istimestamp = diff(istart_out_MM) < 400;
            isfirsttimestamp = diff(istimestamp);
            isfirsttimestamp(isfirsttimestamp == -1) = 0;
            istart_out_MM = istart_out_MM(find(isfirsttimestamp)+1); % samplerate / 10
            
            clear istart_out_map istop_out_map dist_out_map dist_time_map istart_out_MM_clean
            icounter = 1;
            for itimecode = 1:length(istart_out_MM)
                timecodetemplate = markers_MM(istart_out_MM(itimecode):istart_out_MM(itimecode)+2500); % samplerate / 10 * 1.
                [istart_out_map_all,istop_out_map_all,dist_out_map_all] = findsignal(markers_NL,timecodetemplate,'MaxNumSegments',4, 'MaxDistance',0);
                if ~isempty(istart_out_map_all)
                    fprintf('%d/%d: Found %d marker(s) mapping between NeuraLynx and MicroMed, %d selected in total\n', itimecode, length(istart_out_MM), length(istart_out_map_all),icounter);
                    timediff = (micromed_datetime{iMicroMed} + seconds(istart_out_MM(itimecode) / hdr_MicroMed.Fs)) - (neuralynx_starttime + seconds(istart_out_map_all / hdr_MicroMed.Fs)); % use resample frequencies
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
            end
            
            timethreshold       = abs(dist_time_map) < minutes(10);
            dist_time_map       = dist_time_map(timethreshold);
            istart_out_map      = istart_out_map(timethreshold);
            istart_out_MM_clean = istart_out_MM_clean(timethreshold);
            istop_out_map       = istop_out_map(timethreshold);
            dist_out_map        = dist_out_map(timethreshold);
            
            
            figure;
            ploti = 1;
            nrofplots = min(20,length(istart_out_MM_clean));
            for itimecode = 1:nrofplots
                
                timecodetemplate_raw        = sync_MicroMed(istart_out_MM_clean(itimecode):istart_out_MM_clean(itimecode)+450);
                timecodetemplate            = markers_MM(istart_out_MM_clean(itimecode):istart_out_MM_clean(itimecode)+450);
                timecodetemplate_found_raw  = sync_NeuraLynx(istart_out_map(itimecode):istart_out_map(itimecode)+450);
                timecodetemplate_found      = markers_NL(istart_out_map(itimecode):istart_out_map(itimecode)+450);
                
                subplot(nrofplots,4,ploti+ploti*3-3);
                plot(timecodetemplate_raw);         set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight
                
                subplot(nrofplots,4,ploti+ploti*3-2);
                plot(timecodetemplate);             set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight
                
                subplot(nrofplots,4,ploti+ploti*3-1);
                plot(timecodetemplate_found_raw);   set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight
                
                subplot(nrofplots,4,ploti+ploti*3);
                plot(timecodetemplate_found);       set(gca, 'XTickLabel', ''); set(gca, 'YTickLabel', ''); axis tight
                
                ploti = ploti + 1;
                
            end
            
            figure; 
            subplot(4,1,1);     plot(istart_out_MM_clean); title('MicroMed samples of timecodes');
            subplot(4,1,2);     plot(istart_out_map); title('Mapped to Neuralynx samples of timecodes');
            subplot(4,1,3);     plot(istart_out_MM_clean-istart_out_map); title('Distance in samples');
            subplot(4,1,3);     plot(dist_out_map); title('Approx. distance in time');

            sampleshift = median(istart_out_MM_clean-istart_out_map)*4; % calculate resampleing from the beginning
            
            micromed_markers_file{iMicroMed}.Sample_NL = micromed_markers_file{iMicroMed}.Sample + sampleshift;
            micromed_markers_file{iMicroMed}.DataFileFolder_NL = repmat(neurolynx_datafiles{idir}(1).folder,height(micromed_markers_file{iMicroMed}),1);

  
        else
            fprintf('Found no marker in %s within datetime of %s \n',micromed_markers_file{iMicroMed}.DataFileName(1,:), neurolynx_datafiles{idir}(1).name);
        end
    end
end

micromed_markers = outerjoin(micromed_markers_file{:},'MergeKeys', true);

disp('done');
