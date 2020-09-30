function seizure_infos = dtx_stats_seizure_timings(cfg, MuseStruct, ipart)

%mandatory : 
% cfg.seizuretimings.marker_start
% Optional :
% cfg.seizuretimings.marker_end
% cfg.seizuretimings.analysis_start
% cfg.seizuretimings.analysis_end
% cfg.seizuretimings.injection_clock

%get all marker timings
time_start = concatenateMuseMarker(cfg,MuseStruct,ipart, cfg.seizuretimings.marker_start);
if isfield(cfg.seizuretimings, 'marker_end')
    time_end   = concatenateMuseMarker(cfg,MuseStruct,ipart, cfg.seizuretimings.marker_end);
else
    time_end = time_start;
end

%load header
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
for idir = 1:length(MuseStruct{ipart})
    if isNeuralynx
        filelist = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(1).name));
    elseif isMicromed
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']));
    elseif isBrainvision
        hdr{idir} = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']));
    end
end

%get start and end times of analysis
if isfield(cfg.seizuretimings, 'analysis_start')
    analysis_start = concatenateMuseMarker(cfg,MuseStruct,ipart, cfg.seizuretimings.analysis_start);
else %begining of the data
    analysis_start.synctime = 0;
    analysis_start.clock    = MuseStruct{ipart}{1}.starttime;
    analysis_start.dir      = 1;
end
if isfield(cfg.seizuretimings, 'analysis_end')
    analysis_end = concatenateMuseMarker(cfg,MuseStruct,ipart, cfg.seizuretimings.analysis_end);
else %end of the data
    n_samples = 0;
    for idir = 1:size(hdr,2)
        n_samples = n_samples + hdr{idir}.nSamples;
    end
    analysis_end.synctime = n_samples/hdr{1}.Fs;%all data must have the same Fs of course
    analysis_end.clock    = MuseStruct{ipart}{end}.endtime;%all data must have the same Fs of course
    analysis_end.dir      = size(MuseStruct{ipart},2);%all data must have the same Fs of course
end

%safety check
if size(time_start.synctime,2) ~= size(time_end.synctime,2)
    error('%s : there are not the same number of marker_start than marker_end', cfg.prefix(1:end-1));
end
if isempty(analysis_start.synctime) || isempty(analysis_end.synctime)
    error('%s : miss analysis_start or analysis_end', cfg.prefix(1:end-1));
end

%output cfg
seizure_infos                   = cfg.seizuretimings;
seizure_infos.time_start        = time_start;
seizure_infos.time_end          = time_end;
seizure_infos.analysis_start    = analysis_start;
seizure_infos.analysis_end      = analysis_end;

%baseline duration before injection
if isfield(cfg.seizuretimings, 'injection_clock') 
    seizure_infos.baselineduration = cfg.seizuretimings.injection_clock - analysis_start.clock;
end

%record duration
seizure_infos.recordduration = analysis_end.clock - analysis_start.clock;

%post injection duration
if isfield(cfg.seizuretimings, 'injection_clock') 
    seizure_infos.postinjduration = analysis_end.clock - cfg.seizuretimings.injection_clock;
end

%nb of seizures
seizure_infos.nrseizures = size(time_start.synctime,2);

%time between 2 begin of seizures
%remove inter-seizure interval between 2 dirs are they more often occur on
%cut data
seizure_infos.timebetween2seizures      = diff(time_start.clock);
seizure_infos.x_timebetween2seizures    = time_start.clock(2:end); 
changedir                               = logical(diff(time_start.dir));
seizure_infos.timebetween2seizures      = seizure_infos.timebetween2seizures(~changedir); %time from the previous seizure, removing the change of file
seizure_infos.x_timebetween2seizures    = seizure_infos.x_timebetween2seizures(~changedir); %time at the seizure

%seizures duration
if isfield(cfg.seizuretimings, 'marker_end')
    seizure_infos.seizureduration = time_end.clock - time_start.clock;
end

end

