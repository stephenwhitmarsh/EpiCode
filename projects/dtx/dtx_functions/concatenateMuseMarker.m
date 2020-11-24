function marker = concatenateMuseMarker(cfg, MuseStruct, ipart, markername)
% Concatenate markers over all the dirs of one MuseStruct part.
% First line of the output array : clock or synctime
% Second line of the output array : nr of the dir
% clock time : real time, in datetime format
% synctime : time in seconds since the begining of the first file
% dir : integer number of the data directiry (see cfg.directorylist)

% concatenate clock and synctime with different methods. Clock is supposed
% to always be related to the time of the day, so no correction to do.
% synctime is the time since the begining of the file

[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

marker.clock    = datetime.empty;
marker.synctime = [];
marker.dir      = [];
length_previous = 0;

%find hdr for each dir
ft_progress('init','text');
for idir = 1:size(MuseStruct{ipart},2)
    ft_progress(idir/size(MuseStruct{ipart},2),'reading header for dir %d from %d',idir, size(MuseStruct{ipart},2));
    if isNeuralynx
        temp  	 = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{1}, '.ncs']));
        fname 	 = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
    elseif isMicromed
        fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
    elseif isBrainvision
        fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
    end
    hdr{ipart}{idir}     = ft_read_header(fname);
end
ft_progress('close');

%concatenate clock times
for idir = 1:size(MuseStruct{ipart},2) 
    if isfield(MuseStruct{ipart}{idir}.markers, markername)
        if isfield(MuseStruct{ipart}{idir}.markers.(markername), 'clock')
            
            clock_temp = MuseStruct{ipart}{idir}.markers.(markername).clock;
            
            for isample = 1:size(clock_temp,2)
                marker.clock(isample+length_previous) = clock_temp(isample); 
                marker.dir(isample+length_previous) = idir;
            end
            
            length_previous = length_previous + size(clock_temp,2); 
            
        else
            %fprintf('Marker %s exists but is empty in %s\n', markername, MuseStruct{ipart}{idir}.directory);
        end
        
        
    else
        %fprintf('Marker %s does not exist in %s\n',markername, MuseStruct{ipart}{idir}.directory);
    end
    
end

% concatenate synctime.
dirOnset = 0;
for idir = 1:size(MuseStruct{ipart},2)
    if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
        dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
        continue
    end
    if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'synctime')
        dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
        continue
    end
    
    synctime_temp = MuseStruct{ipart}{idir}.markers.(markername).synctime;
    
    for isample = 1:size(synctime_temp,2)
        marker.synctime(end+1) = synctime_temp(isample)+dirOnset/hdr{ipart}{idir}.Fs;
    end
    
    dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
    
end

fprintf('%d times found for %s\n',size(marker.clock,2),markername);



end


