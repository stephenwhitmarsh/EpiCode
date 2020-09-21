function MuseStruct = read_nev_Muse(cfg,MuseStruct)

for ipart = 1:size(MuseStruct,2)
    for idir = 1:size(MuseStruct{ipart},2)
        
        %read nev events
        temppath = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.nev');
        temp     = dir(temppath);
        if size(temp,1) == 1 
            eventname = fullfile(temp(1).folder, temp(1).name);
        else
            error('found 0 or several .nev files');
        end
        
        nev      = read_neuralynx_nev(eventname,'eventformat','neuralynx_nev');
        
        if isempty(nev)
            continue
        end
        
        %read header to convert TimeStamp to seconds
        temp  	 = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.LFP.channel{1}, '.ncs']));
        fname 	 = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);        
        hdr      = ft_read_header(fname);
        
        %remove spaces and minuses
        for ievent = 1:size(nev,2)
            nev(ievent).EventString = strrep(nev(ievent).EventString, ' ', '_');
            nev(ievent).EventString = strrep(nev(ievent).EventString, '-', '_');
        end
      
        %add events to MuseStruct
        cfgtemp                      = [];
        cfgtemp.editmarkerfile.toadd = unique({nev.EventString});
        MuseStruct                   = editMuseMarkers(cfgtemp, MuseStruct);
        
        for markername = string(unique(cfgtemp.editmarkerfile.toadd))
           for ievent = 1:size(nev,2)
               if strcmp(nev(ievent).EventString, markername)
                   
                   MuseStruct{ipart}{idir}.markers.(markername).synctime = ft_getopt(MuseStruct{ipart}{idir}.markers.(markername),'synctime', []);
                   MuseStruct{ipart}{idir}.markers.(markername).clock    = ft_getopt(MuseStruct{ipart}{idir}.markers.(markername),'clock', datetime.empty);
                   
                   marker_time = (nev(ievent).TimeStamp - hdr.FirstTimeStamp) / hdr.TimeStampPerSample / hdr.Fs;
                   MuseStruct{ipart}{idir}.markers.(markername).synctime(end+1) = marker_time;
                   
                   if isempty(MuseStruct{ipart}{idir}.starttime)
                       MuseStruct{ipart}{idir}.markers.(markername).clock(end+1) = NaN;
                   else
                       MuseStruct{ipart}{idir}.markers.(markername).clock(end+1) = MuseStruct{ipart}{idir}.starttime + seconds(marker_time);
                   end
                   
               end
           end
        end
    end
end
