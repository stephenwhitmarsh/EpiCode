function emg_infos = dtx_stats_emg_timings(cfg, MuseStruct_concat,ipart)
% computes some descriptive stats on EMG data

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

%% get emg timings
marker.synctime = [];
marker.clock = datetime.empty;
marker.dir = [];
emg_start = ft_getopt(MuseStruct_concat{ipart}.markers,cfg.emgtimings.emg_start,marker);%concatenateMuseMarker(cfg,MuseStruct_concat,ipart, cfg.emgtimings.emg_start);
emg_end   = ft_getopt(MuseStruct_concat{ipart}.markers,cfg.emgtimings.emg_end,marker);%;%concatenateMuseMarker(cfg,MuseStruct_concat,ipart, cfg.emgtimings.emg_end);
eeg_start = ft_getopt(MuseStruct_concat{ipart}.markers,cfg.emgtimings.eeg_start,marker);

%% get start and end times of analysis
if isfield(cfg.emgtimings, 'analysis_start')
    analysis_start = MuseStruct_concat{ipart}.markers.(cfg.emgtimings.analysis_start);
else %begining of the data
    analysis_start.synctime = 0;
    analysis_start.clock    = MuseStruct_concat{ipart}.starttime;
    analysis_start.dir      = 1;
end
if isfield(cfg.emgtimings, 'analysis_end')
    analysis_end = MuseStruct_concat{ipart}.markers.(cfg.emgtimings.analysis_end);
else 
    analysis_end.synctime = [];%n_samples/hdr{1}.Fs;
    analysis_end.clock    = MuseStruct_concat{ipart}.endtime;
    analysis_end.dir      = size(MuseStruct_concat{ipart}.directory,2);
end

%safety check
if size(emg_start.synctime,2) ~= size(emg_end.synctime,2)
    error('%s : there are not the same number of emg_start than emg_end', cfg.prefix(1:end-1));
end
if isempty(analysis_start.clock) || isempty(analysis_end.clock)
    error('%s : miss analysis_start or analysis_end', cfg.prefix(1:end-1));
end

%% output cfg
emg_infos                   = cfg.emgtimings;
emg_infos.time_start        = emg_start;
emg_infos.time_end          = emg_end;
emg_infos.analysis_start    = analysis_start;
emg_infos.analysis_end      = analysis_end;

%% baseline duration before injection
if isfield(cfg.emgtimings, 'injection_clock') 
    emg_infos.baselineduration = cfg.emgtimings.injection_clock - analysis_start.clock;
end

%% record duration
emg_infos.recordduration = analysis_end.clock - analysis_start.clock;

%% post injection duration
if isfield(cfg.emgtimings, 'injection_clock') 
    emg_infos.postinjduration = analysis_end.clock - cfg.emgtimings.injection_clock;
end

%% nb of emg
emg_infos.nr_emg = size(emg_start.synctime,2);

%% emg duration
if ~isempty(emg_end.clock)
    emg_infos.emg_duration = seconds(emg_end.clock - emg_start.clock);
else
    emg_infos.emg_duration = [];
end

%% CV2 emg duration
arraytemp = emg_infos.emg_duration;
cv2_data = NaN;
if length(arraytemp)>2
    for i = 1:length(arraytemp)-1
        cv2_data(i) = 2*abs(arraytemp(i)-arraytemp(i+1))/(arraytemp(i)+arraytemp(i+1));
    end
end
clear arraytemp
emg_infos.emg_duration_cv2 = cv2_data;

%% eeg emg delay
if ~isempty(emg_start.clock)
    emg_infos.eeg_emg_delay = seconds(emg_start.clock - eeg_start.clock);
else
    emg_infos.eeg_emg_delay = [];
end

%% CV2 eeg emg delay
arraytemp = emg_infos.eeg_emg_delay;
cv2_data = NaN;
if length(arraytemp)>2
    for i = 1:length(arraytemp)-1
        cv2_data(i) = 2*abs(arraytemp(i)-arraytemp(i+1))/(arraytemp(i)+arraytemp(i+1));
    end
end
clear arraytemp
emg_infos.eeg_emg_delay_cv2 = cv2_data;

%% same values over a sliding time window : temps entre 2
emg_infos.statsovertime.winsize = cfg.emgtimings.winsize;
emg_infos.statsovertime.winstep = cfg.emgtimings.winstep;
emg_infos.statsovertime.time_unit = 'seconds';

%first window
twin_start = analysis_start.clock;
twin_end   = analysis_start.clock + seconds(cfg.emgtimings.winsize);
i_window   = 0;

%go trough each window
while twin_end < analysis_end.clock
    
    i_window   = i_window+1;
    
    %keep window time
    if isfield(cfg.emgtimings, 'injection_clock')
        %time relative to injection
        emg_infos.statsovertime.starttime_orig(i_window) = twin_start;
        emg_infos.statsovertime.endtime_orig(i_window)   = twin_end;
        emg_infos.statsovertime.starttime(i_window)      = twin_start - cfg.emgtimings.injection_clock;
        emg_infos.statsovertime.endtime(i_window)        = twin_end - cfg.emgtimings.injection_clock;
    else
        emg_infos.statsovertime.starttime(i_window)      = twin_start;
        emg_infos.statsovertime.endtime(i_window)        = twin_end;
    end
       
    emg_infos.statsovertime.nb_emg(i_window) = sum(emg_start.clock > twin_start & emg_start.clock < twin_end);
    
    for iparam= ["emg_duration", "eeg_emg_delay"]
        %mean and std of param
        idx = emg_infos.time_start.clock > twin_start & emg_infos.time_end.clock < twin_end;
        emg_infos.statsovertime.(sprintf('%s_nb_emg',iparam))(i_window) = sum(idx);
        emg_infos.statsovertime.(sprintf('%s_mean',iparam))(i_window) = nanmean(emg_infos.(iparam)(idx));
        emg_infos.statsovertime.(sprintf('%s_std',iparam))(i_window)  = nanstd(emg_infos.(iparam)(idx));
        %cv2 of param
        idx = idx(1:end-1); %one less value for cv2
        emg_infos.statsovertime.(sprintf('%s_cv2',iparam))(i_window) = nanmean(emg_infos.(sprintf('%s_cv2',iparam))(idx));
    end

    %remove values if this is part of missing file (see dtx_eegvideo_setparams.m)
    keepwindow = true;
    if ~isempty(cfg.missingdata)
        if twin_start > cfg.missingdata(1) && twin_start > cfg.missingdata(2)
            keepwindow = true;
        elseif twin_end < cfg.missingdata(1) && twin_end < cfg.missingdata(2)
            keepwindow = true;
        else
            keepwindow = false;
        end
    end
    if ~keepwindow
        emg_infos.statsovertime.nb_seizures(i_window) = nan;
        for iparam= ["emg_duration", "eeg_emg_delay"]
            %mean and std of param
            emg_infos.statsovertime.(sprintf('%s_nb_seizures',iparam))(i_window) = nan;
            emg_infos.statsovertime.(sprintf('%s_mean',iparam))(i_window) = nan;
            emg_infos.statsovertime.(sprintf('%s_std',iparam))(i_window)  = nan;
            %cv2 of param
            emg_infos.statsovertime.(sprintf('cv2_%s_mean',iparam))(i_window) = nan;
            emg_infos.statsovertime.(sprintf('cv2_%s_std',iparam))(i_window)  = nan;
        end
    end
    
    %go to the next window
    twin_start = twin_start + seconds(cfg.emgtimings.winstep);
    twin_end   = twin_end   + seconds(cfg.emgtimings.winstep);
    
end %while

%% save to disk
emg_infos = orderfields(emg_infos);
emg_infos.statsovertime = orderfields(emg_infos.statsovertime);
save(fullfile(cfg.datasavedir,[cfg.prefix,'emg_infos.mat']), 'emg_infos');

end

