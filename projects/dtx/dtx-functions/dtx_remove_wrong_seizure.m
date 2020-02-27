function [MuseStruct_corrected] = dtx_remove_wrong_seizure(cfg, MuseStruct, force)
%check if there is slowave without seizure or seizure without slowwave and
%remove the correspondant marker.
%A seizure beginning in one file and ending in the next file is ignored
%Add indexes and number of removed markers in MuseStruct
%Paul Baudin


fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct.mat']);
if exist(fname,'file') && force == false
    fprintf('******************************\n');
    fprintf('****** Loading results  ******\n');
    fprintf('******************************\n\n');
    load(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_corrected.mat']),'MuseStruct_corrected');
else
    
    if exist(fname,'file') && force == true
        fprintf('********************************\n');
        fprintf('******** Forced redoing ********\n');
        fprintf('********************************\n\n');
    else
        fprintf('*************************\n');
        fprintf('*** Removing seizures ***\n');
        fprintf('*************************\n\n');
    end
    
    MuseStruct_corrected = MuseStruct;
    
    for ipart = 1 : length(MuseStruct)
        
        for idir = 1:length(MuseStruct{ipart})
            
            if isfield(MuseStruct{ipart}{idir},'markers')
                if isfield(MuseStruct{ipart}{idir}.markers,'SlowWave')
                    if isfield(MuseStruct{ipart}{idir}.markers.SlowWave,'clock')
                        if ~isempty(MuseStruct{ipart}{idir}.markers.SlowWave.clock) %If no clock : don't change the structure
                            
                            n_ignoreSeizure{idir} = 0;
                            n_ignoreSlowWave{idir} = 0;
                            n_ignoreLastSeizure{idir} = 0;
                            
                            ignoreSeizure_index = [];
                            ignoreSlowWave_index = [];
                            ignoreLastSeizure_index = [];
                            
                            SlowWave_orig = MuseStruct{ipart}{idir}.markers.SlowWave.clock;
                            Crise_Start_orig = MuseStruct{ipart}{idir}.markers.Crise_Start.clock;
                            Crise_End_orig = MuseStruct{ipart}{idir}.markers.Crise_End.clock;
                            
                            %ignore the first crise_end if it is the end of a
                            %seizure begining in the previous dir
                            if Crise_End_orig(1)-Crise_Start_orig(1)<0
                                ignore_1st_Crise_End{idir}=1;
                            else
                                ignore_1st_Crise_End{idir}=0;
                            end
                            
                            %go throught all markers
                            iSlowWave = 1;
                            iSeizure = 1;
                            iresult = 0;
                            
                            while iSeizure <= length(Crise_Start_orig)  && iSlowWave <= length(SlowWave_orig)
                                iresult = iresult+1;
                                ignoreseizure = 0;
                                
                                %SlowWave without seizure
                                if seconds(Crise_Start_orig(iSeizure)-SlowWave_orig(iSlowWave)) > 10
                                    n_ignoreSlowWave{idir} = n_ignoreSlowWave{idir}+1;
                                    ignoreSlowWave_index(n_ignoreSlowWave{idir}) = iSlowWave;
                                    iSlowWave = iSlowWave+1;
                                    ignoreseizure = 1;
                                end
                                
                                %Seizure without SlowWave
                                if ignoreseizure == 0 && seconds(Crise_Start_orig(iSeizure)-SlowWave_orig(iSlowWave)) <0
                                    n_ignoreSeizure{idir} = n_ignoreSeizure{idir}+1;
                                    ignoreSeizure_index(n_ignoreSeizure{idir}) = iSeizure;
                                    iSeizure = iSeizure+1;
                                    ignoreseizure = 1;
                                end
                                
                                %ignore the last seizure of the file if seizure ends in the next file
                                if length(Crise_Start_orig)-length(Crise_End_orig)==1 && iSeizure==length(Crise_Start_orig)
                                    n_ignoreLastSeizure{idir} = n_ignoreLastSeizure{idir}+1;
                                    ignoreLastSeizure_index(n_ignoreLastSeizure{idir}) = idir;
                                    iresult = iresult-1;
                                    ignore_last_Seizure{idir}=1;
                                    break
                                else
                                    ignore_last_Seizure{idir}=0;
                                end
                                
                                %Write markers to MuseStruct_corrected if seizure is not to ignore
                                if ignoreseizure == 0
                                    
                                    marker = ["SlowWave", "Crise_Start", "Crise_End"];
                                    eventindex = {iSlowWave, iSeizure, iSeizure+ignore_1st_Crise_End{idir}};
                                    
                                    for i=1:3
                                        imarker = marker(i);
                                        i_eventindex = eventindex{i};
                                        for ifield = ["synctime", "clock"]
                                            MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(iresult) = ...
                                                MuseStruct{ipart}{idir}.markers.(imarker).(ifield)(i_eventindex);
                                        end
                                    end
                                    
                                    iSlowWave = iSlowWave+1;
                                    iSeizure = iSeizure+1;
                                end
                                
                                
                            end %while

                            %Remove values after the last good seizure
                            for imarker = ["SlowWave", "Crise_Start", "Crise_End"]
                                for ifield = ["synctime", "clock"]
                                    MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield) = ...
                                        MuseStruct_corrected{ipart}{idir}.markers.(imarker).(ifield)(1:iresult);
                                end
                            end
                            
                            %Add index values
                            marker = ["SlowWave", "Crise_Start", "Crise_End"];
                            if ignore_1st_Crise_End{idir}
                                ignoreindex = {ignoreSlowWave_index, ignoreSeizure_index, [ignore_1st_Crise_End{idir}, ignoreSeizure_index]};
                            else
                                ignoreindex = {ignoreSlowWave_index, ignoreSeizure_index, ignoreSeizure_index};
                            end
                            for i=1:3
                                imarker = marker(i);
                                i_ignoreindex = ignoreindex{i};
                                MuseStruct_corrected{ipart}{idir}.markers.(imarker).originalindex_ignored = i_ignoreindex;
                            end
                            if ignore_last_Seizure{idir}
                                MuseStruct_corrected{ipart}{idir}.markers.SlowWave.originalindex_ignored(end+1) = length(SlowWave_orig);
                                MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.originalindex_ignored(end+1) = length(Crise_Start_orig);
                            end
                            
                            
                            %add infos of nr of removed seizure per dir
                            MuseStruct_corrected{ipart}{idir}.ignoremarkers_Info.ignoreLastSeizure = ignore_last_Seizure{idir};
                            MuseStruct_corrected{ipart}{idir}.ignoremarkers_Info.ignoreFirstCrise_End = ignore_1st_Crise_End{idir};
                            MuseStruct_corrected{ipart}{idir}.ignoremarkers_Info.n_ignored_SlowWave = n_ignoreSlowWave{idir}+ignore_last_Seizure{idir};
                            MuseStruct_corrected{ipart}{idir}.ignoremarkers_Info.n_ignored_Crise_Start = n_ignoreSeizure{idir}+ignore_last_Seizure{idir};
                            MuseStruct_corrected{ipart}{idir}.ignoremarkers_Info.n_ignored_Crise_End = n_ignoreSeizure{idir}+ignore_1st_Crise_End{idir};
                            
                            %safety check
                            if ~(length(MuseStruct_corrected{ipart}{idir}.markers.SlowWave.synctime) == length(MuseStruct_corrected{ipart}{idir}.markers.Crise_Start.synctime) &&...
                                    length(MuseStruct_corrected{ipart}{idir}.markers.SlowWave.synctime) == length(MuseStruct_corrected{ipart}{idir}.markers.Crise_End.synctime))
                                error('something wrong with the removal of markers. %s part %d dir %d', cfg.prefix(1:end-1), ipart, idir);
                            end
                            
                        end
                    end
                end
            end
            
        end %idir
        
    end %ipart
    
end

% save data
save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_corrected.mat']),'MuseStruct_corrected');

end









