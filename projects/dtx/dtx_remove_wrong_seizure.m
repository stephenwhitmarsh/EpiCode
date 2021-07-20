function [MuseStruct_corrected] = dtx_remove_wrong_seizure(cfg, MuseStruct, remove_seizure_between_2_files)

cfg.type = ft_getopt(cfg, 'type');
if ~strcmp(cfg.type, 'dtx')
    MuseStruct_corrected = MuseStruct;
    return
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_corrected_dtx.mat']);
MuseStruct_corrected = MuseStruct;

fprintf('For %s\n',cfg.prefix(1:end-1));
for ipart = 1 : length(MuseStruct)
    for idir = 1:length(MuseStruct{ipart})
        if ~isfield(MuseStruct{ipart}{idir},'markers')
            continue
        end
        if ~isfield(MuseStruct{ipart}{idir}.markers,'SlowWave')
            continue
        end
        if ~isfield(MuseStruct{ipart}{idir}.markers.SlowWave,'clock')
            continue
        end
        if isempty(MuseStruct{ipart}{idir}.markers.SlowWave.clock) 
            continue
        end
        
        n_removeSeizure{idir} = 0;
        n_removeSlowWave{idir} = 0;
        removeSeizure_index = [];
        removeSlowWave_index = [];
        n_last_Crise_Start_removed{idir} = 0; 
        n_last_Crise_End_removed{idir} = 0;
        
        SlowWave_orig = MuseStruct{ipart}{idir}.markers.SlowWave.clock;
        Crise_Start_orig = MuseStruct{ipart}{idir}.markers.Crise_Start.clock;
        Crise_End_orig = MuseStruct{ipart}{idir}.markers.Crise_End.clock;
        
        %remove the first crise_end if it is the end of a
        %seizure begining in the previous dir
        if Crise_End_orig(1)-Crise_Start_orig(1)<0
            if remove_seizure_between_2_files
                fprintf('part %d dir %d : remove first Crise_End\n', ipart, idir);
                cfgtemp                     = [];
                cfgtemp.bad.part_list       = ipart;
                cfgtemp.bad.dir_list        = idir;
                cfgtemp.bad.markerStart     = "Crise_End";
                cfgtemp.bad.markerEnd       = "Crise_End";
                cfgtemp.bad.sample_list     = 1;
                cfgtemp.bad.time_from_begin = -1;
                cfgtemp.bad.time_from_end   = 1;
                [MuseStruct_corrected]      = addMuseBAD(cfgtemp,MuseStruct_corrected);
            end
            remove_1st_Crise_End{idir}=1; %correct index for other seizures
            
        else
            remove_1st_Crise_End{idir}=0;
        end

        %go throught all markers
        iSlowWave = 1;
        iSeizure = 1;
        iresult = 0;
        
        while iSeizure <= length(Crise_Start_orig)  && iSlowWave <= length(SlowWave_orig)
            removeseizure = 0;
            
            %SlowWave without seizure
            if seconds(Crise_Start_orig(iSeizure)-SlowWave_orig(iSlowWave)) > 10
                fprintf('part %d dir %d : remove %dth SlowWave, because of no seizure\n', ipart, idir, iSlowWave);
                cfgtemp                     = [];
                cfgtemp.bad.part_list       = ipart;
                cfgtemp.bad.dir_list        = idir;
                cfgtemp.bad.markerStart     = "SlowWave";
                cfgtemp.bad.markerEnd       = "SlowWave";
                cfgtemp.bad.sample_list     = iSlowWave;
                cfgtemp.bad.time_from_begin = -1;
                cfgtemp.bad.time_from_end   = 1;
                [MuseStruct_corrected]      = addMuseBAD(cfgtemp,MuseStruct_corrected);
                n_removeSlowWave{idir} = n_removeSlowWave{idir}+1;
                removeSlowWave_index(n_removeSlowWave{idir}) = iSlowWave;
                iSlowWave = iSlowWave+1; %go to the next SlowWave
                removeseizure = 1;
            end
            
            %Seizure without SlowWave
            if removeseizure == 0 && seconds(Crise_Start_orig(iSeizure)-SlowWave_orig(iSlowWave)) <0
                fprintf('part %d dir %d : remove %dth seizure, because of no SlowWave\n', ipart, idir, iSeizure);
                cfgtemp                     = [];
                cfgtemp.bad.part_list       = ipart;
                cfgtemp.bad.dir_list        = idir;
                cfgtemp.bad.markerStart     = "Crise_Start";
                cfgtemp.bad.markerEnd       = "Crise_End";
                cfgtemp.bad.sample_list     = iSeizure;
                cfgtemp.bad.time_from_begin = -1;
                cfgtemp.bad.time_from_end   = 1;
                [MuseStruct_corrected]      = addMuseBAD(cfgtemp,MuseStruct_corrected);
                n_removeSeizure{idir} = n_removeSeizure{idir}+1;
                removeSeizure_index(n_removeSeizure{idir}) = iSeizure;
                iSeizure = iSeizure+1; %go to the next seizure
                removeseizure = 1;
            end
            
            %remove the last seizure of the file if seizure ends in the next file
            if length(Crise_Start_orig)-length(Crise_End_orig)==1 && iSeizure==length(Crise_Start_orig)
                if remove_seizure_between_2_files
                    fprintf('part %d dir %d : remove last seizure, because cut into 2 files\n', ipart, idir);
                    cfgtemp                     = [];
                    cfgtemp.bad.part_list       = ipart;
                    cfgtemp.bad.dir_list        = idir;
                    cfgtemp.bad.markerStart     = "Crise_Start";
                    cfgtemp.bad.markerEnd       = "Crise_Start";
                    cfgtemp.bad.sample_list     = iSeizure;
                    cfgtemp.bad.time_from_begin = -1;
                    cfgtemp.bad.time_from_end   = 1;
                    [MuseStruct_corrected]      = addMuseBAD(cfgtemp,MuseStruct_corrected);
                    remove_cut_Seizure{idir}=1;
                    break
                else %add CriseStart and SlowWave without Crise_End
                    remove_cut_Seizure{idir}=0;
                    iresult = iresult+1;
                    marker = ["SlowWave", "Crise_Start"];
                    eventindex = {iSlowWave, iSeizure};
                    for i=1:2
                        imarker = marker(i);
                        i_eventindex = eventindex{i};
                        for ifield = ["synctime", "clock"]
                            MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(iresult) = ...
                                MuseStruct{ipart}{idir}.markers.(imarker).(ifield)(i_eventindex);
                        end
                    end
                    
                end
                break %break because it is the last seizure
            else
                remove_cut_Seizure{idir}=0;
            end
            
            %Write markers to MuseStruct_corrected if seizure is not to remove
            if removeseizure == 0
                iresult = iresult+1;
                marker = ["SlowWave", "Crise_Start", "Crise_End"];
                eventindex = {iSlowWave, iSeizure, iSeizure+remove_1st_Crise_End{idir}};
                for i=1:3
                    imarker = marker(i);
                    i_eventindex = eventindex{i};
                    for ifield = ["synctime", "clock"]
                        MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(iresult) = ...
                            MuseStruct{ipart}{idir}.markers.(imarker).(ifield)(i_eventindex);
                    end
                end
                %go to the next SlowWave and the next
                %Seizure
                iSlowWave = iSlowWave+1;
                iSeizure = iSeizure+1;
            end
        end %while
        
        %Add index values to MuseStruct
        marker = ["SlowWave", "Crise_Start", "Crise_End"];
        if remove_1st_Crise_End{idir}
            removeindex = {removeSlowWave_index, removeSeizure_index, [remove_1st_Crise_End{idir}, removeSeizure_index]};
        else
            removeindex = {removeSlowWave_index, removeSeizure_index, removeSeizure_index};
        end
        
        for i=1:3
            imarker = marker(i);
            i_removeindex = removeindex{i};
            MuseStruct_corrected{ipart}{idir}.markers.(imarker).originalindex_removed = i_removeindex;
        end
        if remove_cut_Seizure{idir}
            MuseStruct_corrected{ipart}{idir}.markers.SlowWave.originalindex_removed(end+1) = length(SlowWave_orig);
            MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.originalindex_removed(end+1) = length(Crise_Start_orig);
        end
        
        %Remove wrong marker if there are any, after the last good seizure (so not removed by the while loop):
        
        %Find for SlowWaves
        if iresult + n_removeSlowWave{idir} + remove_cut_Seizure{idir} < length(SlowWave_orig)
            n_lastSlowWaves_removed = length(SlowWave_orig)-iresult;
            n_removeSlowWave{idir} = n_removeSlowWave{idir}+n_lastSlowWaves_removed;
            index_temp = iresult+1 : length(SlowWave_orig);
            cfgtemp                     = [];
            cfgtemp.bad.part_list       = ipart;
            cfgtemp.bad.dir_list        = idir;
            cfgtemp.bad.markerStart     = "SlowWave";
            cfgtemp.bad.markerEnd       = "SlowWave";
            cfgtemp.bad.sample_list     = index_temp;
            cfgtemp.bad.time_from_begin = -1;
            cfgtemp.bad.time_from_end   = 1;
            [MuseStruct_corrected]      = addMuseBAD(cfgtemp,MuseStruct_corrected);
            MuseStruct_corrected{ipart}{idir}.markers.SlowWave.originalindex_removed(end+1:end+length(index_temp)) = index_temp;
            fprintf('part %d dir %d : remove %d last SlowWaves, because of no seizure\n', ipart, idir,n_lastSlowWaves_removed);
        end
        
        %Find for Seizures
        if iresult + n_removeSeizure{idir} + remove_cut_Seizure{idir} < length(Crise_Start_orig)
            index_temp_start = iresult+1 : length(Crise_Start_orig);
            index_temp_end = iresult+1 : length(Crise_End_orig);
            n_lastSeizures_removed = length(Crise_Start_orig)-iresult;
            n_last_Crise_Start_removed{idir} = length(index_temp_start);
            n_last_Crise_End_removed{idir} = length(index_temp_end);
            cfgtemp                     = [];
            cfgtemp.bad.part_list       = ipart;
            cfgtemp.bad.dir_list        = idir;
            cfgtemp.bad.markerStart     = "Crise_Start";
            cfgtemp.bad.markerEnd       = "Crise_Start";
            cfgtemp.bad.sample_list     = index_temp_start;
            cfgtemp.bad.time_from_begin = -1;
            cfgtemp.bad.time_from_end   = 1;
            [MuseStruct_corrected] = addMuseBAD(cfgtemp,MuseStruct_corrected);
            cfgtemp.bad.markerStart     = "Crise_End";
            cfgtemp.bad.markerEnd       = "Crise_End";
            cfgtemp.bad.sample_list     = index_temp_end;
            [MuseStruct_corrected] = addMuseBAD(cfgtemp,MuseStruct_corrected);
            MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.originalindex_removed(end+1:end+length(index_temp_start)) = index_temp_start;
            MuseStruct_corrected{ipart}{idir}.markers.Crise_End.originalindex_removed(end+1:end+length(index_temp_end)) = index_temp_end;
            fprintf('part %d dir %d : remove %d last seizures, because of no SlowWave\n', ipart, idir, n_lastSeizures_removed);
        end
        
        %Remove events
        if ~remove_seizure_between_2_files && length(MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.synctime) - length(MuseStruct_corrected{ipart}{idir}.markers.Crise_End.synctime) == 1
            %if last cut seizure is not
            %removed, there is one Crise_End
            %missing
            imarker = 'Crise_End';
            for ifield = ["synctime", "clock"]
                MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield) = ...
                    MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(1:iresult-1);
            end
            marker_list = [];
            marker_list = ["SlowWave", "Crise_Start"];
        elseif ~remove_seizure_between_2_files && length(MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.synctime) - length(MuseStruct_corrected{ipart}{idir}.markers.Crise_End.synctime) == -1
            % if begins with a Crise_End alone, leave it
            imarker = 'Crise_End';
            for ifield = ["synctime", "clock"]
                MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield) = ...
                    MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(1:iresult+1);
            end
            marker_list = [];
            marker_list = ["SlowWave", "Crise_Start"];
        else
            marker_list = [];
            marker_list = ["SlowWave", "Crise_Start", "Crise_End"];
        end
        for imarker = marker_list
            for ifield = ["synctime", "clock"]
                MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield) = ...
                    MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(1:iresult);
            end
        end
        
        %add infos of nr of removed seizure per dir in MuseStruct
        MuseStruct_corrected{ipart}{idir}.removemarkers_Info.removeLastSeizure = remove_cut_Seizure{idir};
        MuseStruct_corrected{ipart}{idir}.removemarkers_Info.removeFirstCrise_End = remove_1st_Crise_End{idir};
        MuseStruct_corrected{ipart}{idir}.removemarkers_Info.n_removed_SlowWave = n_removeSlowWave{idir}+remove_cut_Seizure{idir};
        MuseStruct_corrected{ipart}{idir}.removemarkers_Info.n_removed_Crise_Start = n_removeSeizure{idir}+remove_cut_Seizure{idir} + n_last_Crise_Start_removed{idir};
        MuseStruct_corrected{ipart}{idir}.removemarkers_Info.n_removed_Crise_End = n_removeSeizure{idir}+remove_1st_Crise_End{idir} + n_last_Crise_End_removed{idir};
        
        %safety check
        if remove_seizure_between_2_files
            if ~(length(MuseStruct_corrected{ipart}{idir}.markers.SlowWave.synctime) == length(MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.synctime) &&...
                    length(MuseStruct_corrected{ipart}{idir}.markers.SlowWave.synctime) == length(MuseStruct_corrected{ipart}{idir}.markers.Crise_End.synctime))
                error('something wrong with the removal of markers. %s part %d dir %d', cfg.prefix(1:end-1), ipart, idir);
            end
        end

    end %idir
end %ipart