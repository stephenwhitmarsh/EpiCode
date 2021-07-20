function seizure_infos = dtx_stats_seizure_timings(cfg, MuseStruct_concat,ipart)
% computes some descriptive stats for data with or without EMG, with
% seizure start/end (rats) or slowwave (patients). For rats, time is 
% relative to the injection of DTX.

% ## INPUT : 
% - cfg :
%mandatory : 
% cfg.seizuretimings.marker_start
% cfg.seizuretimings.winsize        
% cfg.seizuretimings.winstep        
% 
% Optional :
% cfg.seizuretimings.marker_end
% cfg.seizuretimings.analysis_start
% cfg.seizuretimings.analysis_end
% cfg.seizuretimings.injection_clock
% 
% - MuseStruct_concat : output from concatenateMuseMarkers
% - ipart : nr. of the analysis part concerned (see cfg.directorylist in the setparams script) 

%% get marker timings
time_start = MuseStruct_concat{ipart}.markers.(cfg.seizuretimings.marker_start);
if isfield(cfg.seizuretimings, 'marker_end')
    time_end   = MuseStruct_concat{ipart}.markers.(cfg.seizuretimings.marker_end);
else
    time_end = time_start;
end

%% get start and end times of analysis
if isfield(cfg.seizuretimings, 'analysis_start')
    analysis_start = MuseStruct_concat{ipart}.markers.(cfg.seizuretimings.analysis_start);
else %begining of the data
    analysis_start.synctime = 0;
    analysis_start.clock    = MuseStruct_concat{ipart}.starttime;
    analysis_start.dir      = 1;
end

if size(analysis_start.synctime,2) >1
    analysis_start.synctime = analysis_start.synctime(1);
    analysis_start.clock    = analysis_start.clock(1);
    analysis_start.dir      = analysis_start.dir(1);
end

if isfield(cfg.seizuretimings, 'analysis_end')
    analysis_end = MuseStruct_concat{ipart}.markers.(cfg.seizuretimings.analysis_end);
else 
    analysis_end.synctime = [];
    analysis_end.clock    = MuseStruct_concat{ipart}.endtime;
    analysis_end.dir      = size(MuseStruct_concat{ipart}.directory,2);
end

%safety check
if size(time_start.synctime,2) ~= size(time_end.synctime,2)
    error('%s : there are not the same number of marker_start than marker_end', cfg.prefix(1:end-1));
end
if isempty(analysis_start.clock) || isempty(analysis_end.clock)
    error('%s : miss analysis_start or analysis_end', cfg.prefix(1:end-1));
end

%% output infos from cfg
seizure_infos                   = cfg.seizuretimings;
seizure_infos.time_start        = time_start;
seizure_infos.time_end          = time_end;
seizure_infos.analysis_start    = analysis_start;
seizure_infos.analysis_end      = analysis_end;

%% baseline duration before injection
if isfield(cfg.seizuretimings, 'injection_clock') 
    seizure_infos.baselineduration = cfg.seizuretimings.injection_clock - analysis_start.clock;
end

%% record duration
seizure_infos.recordduration = analysis_end.clock - analysis_start.clock;

%% post injection duration
if isfield(cfg.seizuretimings, 'injection_clock') 
    seizure_infos.postinjduration = analysis_end.clock - cfg.seizuretimings.injection_clock;
end

%% nb of seizures
seizure_infos.nrseizures = size(time_start.synctime,2);

%% time between 2 begins of seizures
%remove inter-seizure interval between 2 dirs are there may have missing
%data
seizure_infos.timebetween2seizures      = diff(time_start.clock);
seizure_infos.x_timebetween2seizures    = time_start.clock(2:end); 
changedir                               = logical(diff(time_start.dir));
seizure_infos.timebetween2seizures      = seizure_infos.timebetween2seizures(~changedir); %time from the previous seizure, removing the change of file
seizure_infos.x_timebetween2seizures    = seizure_infos.x_timebetween2seizures(~changedir); %time at the seizure

%% seizures duration
if isfield(cfg.seizuretimings, 'marker_end')
    seizure_infos.seizureduration = time_end.clock - time_start.clock;
    seizure_infos.x_seizureduration = time_start.clock;
end

%% CV2 time between 2 seizures 
datatemp = seizure_infos.timebetween2seizures;
cv2_data = NaN;
if length(datatemp)>2
    for i = 1:length(datatemp)-1
        cv2_data(i) = 2*abs(datatemp(i)-datatemp(i+1))/(datatemp(i)+datatemp(i+1));
    end
end
clear datatemp
seizure_infos.timebetween2seizures_cv2 = cv2_data;

%% CV2 seizure duration
if isfield(cfg.seizuretimings, 'marker_end')
    datatemp = seizure_infos.seizureduration;
    cv2_data = NaN;
    if length(datatemp)>2
        for i = 1:length(datatemp)-1
            cv2_data(i) = 2*abs(datatemp(i)-datatemp(i+1))/(datatemp(i)+datatemp(i+1));
        end
    end
    clear datatemp
    seizure_infos.seizureduration_cv2 = cv2_data;
end

%% same values but over a sliding time window 
if ~isempty(ft_getopt(cfg.seizuretimings, 'winsize'))
    seizure_infos.statsovertime.winsize = cfg.seizuretimings.winsize;
    seizure_infos.statsovertime.winstep = cfg.seizuretimings.winstep;
    seizure_infos.statsovertime.time_unit = 'seconds';
    
    %first window
    twin_start = analysis_start.clock;
    twin_end   = analysis_start.clock + seconds(cfg.seizuretimings.winsize);
    i_window   = 0;
    
    %go trough each window
    while twin_end < analysis_end.clock
        
        i_window   = i_window+1;
        
        %keep window time
        if isfield(cfg.seizuretimings, 'injection_clock')
            %time relative to injection
            seizure_infos.statsovertime.starttime_orig(i_window) = twin_start;
            seizure_infos.statsovertime.endtime_orig(i_window)   = twin_end;
            seizure_infos.statsovertime.starttime(i_window)      = twin_start - cfg.seizuretimings.injection_clock;
            seizure_infos.statsovertime.endtime(i_window)        = twin_end - cfg.seizuretimings.injection_clock;
        else
            seizure_infos.statsovertime.starttime(i_window)      = twin_start;
            seizure_infos.statsovertime.endtime(i_window)        = twin_end;
        end
        
        seizure_infos.statsovertime.nb_seizures(i_window) = sum(time_start.clock > twin_start & time_start.clock < twin_end);
        
        for iparam= ["timebetween2seizures", "seizureduration"]
            %mean and std of param
            idx = seizure_infos.(sprintf('x_%s',iparam)) > twin_start & seizure_infos.(sprintf('x_%s',iparam)) < twin_end;
            seizure_infos.statsovertime.(sprintf('%s_nb_seizures',iparam))(i_window) = sum(idx);
            seizure_infos.statsovertime.(sprintf('%s_mean',iparam))(i_window) = nanmean(seconds(seizure_infos.(iparam)(idx)));
            seizure_infos.statsovertime.(sprintf('%s_std',iparam))(i_window)  = nanstd(seconds(seizure_infos.(iparam)(idx)));
            %cv2 of param
            idx = idx(1:end-1); %one less value for cv2
            seizure_infos.statsovertime.(sprintf('%s_cv2',iparam))(i_window) = nanmean(seizure_infos.(sprintf('%s_cv2',iparam))(idx));
        end
        
        % remove values if this is part of a data hole (it is an exception 
        % for a few data, see dtx_eegvideo_setparams.m), otherwise it will be artefactely
        % zero
        keepwindow = true;
        cfg.missingdata = ft_getopt(cfg, 'missingdata');
        if ~isempty(cfg.missingdata)
            if all(twin_start > cfg.missingdata)
                keepwindow = true;
            elseif all(twin_end < cfg.missingdata)
                keepwindow = true;
            else
                keepwindow = false;
            end
        end
        if ~keepwindow
            seizure_infos.statsovertime.nb_seizures(i_window) = nan;
            for iparam= ["timebetween2seizures", "seizureduration"]
                %mean and std of param
                seizure_infos.statsovertime.(sprintf('%s_nb_seizures',iparam))(i_window) = nan;
                seizure_infos.statsovertime.(sprintf('%s_mean',iparam))(i_window) = nan;
                seizure_infos.statsovertime.(sprintf('%s_std',iparam))(i_window)  = nan;
                %cv2 of param
                seizure_infos.statsovertime.(sprintf('%s_cv2',iparam))(i_window) = nan;
            end
        end
        
        %go to the next window
        twin_start = twin_start + seconds(cfg.seizuretimings.winstep);
        twin_end   = twin_end   + seconds(cfg.seizuretimings.winstep);
        
    end %while
end

%% save to disk
seizure_infos = orderfields(seizure_infos);
try, seizure_infos.statsovertime = orderfields(seizure_infos.statsovertime); end
save(fullfile(cfg.datasavedir,[cfg.prefix,'seizures_infos.mat']), 'seizure_infos');

end

