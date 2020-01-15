function Muse2txt(MuseStruct)

%% Text files for overview markers

startend = {'BAD__START__','BAD__END__';...
    'VF1__START__','VF1__END__';...
    'VF3__START__','VF3__END__';...
    'VF4__START__','VF4__END__';...
    'VF5__START__','VF5__END__';...
    'P','P';...
    'PP__START__','PP__END__';...
    'RR__START__','RR__END__'};

for imarker = 1 : length(startend)
    Starttime       = [];
    Endtime         = [];
    Startsample     = [];
    Endsample       = [];
    Duration        = [];
    Directory       = [];
    Event_start     = [];
    Event_end       = [];
    for idir = 1 : size(MuseStruct,2)
        try % some might be empty
            for ievent = 1 : size(MuseStruct{idir}.markers.(startend{imarker,1}).events,2)
                Event_start = [Event_start; startend{imarker,1}];
                Event_end   = [Event_end; startend{imarker,2}];
                Directory   = [Directory;   MuseStruct{idir}.directory];
                Starttime   = [Starttime;   MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).time];
                Endtime     = [Endtime;     MuseStruct{idir}.markers.(startend{imarker,2}).events(ievent).time];
                Duration    = [Duration;    MuseStruct{idir}.markers.(startend{imarker,2}).events(ievent).time - MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).time];
                Startsample = [Startsample; MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).sample];
                Endsample   = [Endsample;   MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).sample];
            end
        catch
        end
    end
    overview = table(Event_start,Event_end,Directory,Starttime,Endtime,Duration,Startsample,Endsample);
    writetable(overview,['/Users/stephen.whitmarsh/Documents/VF/',startend{imarker,1} '-to-' startend{imarker,2} '.txt']);
end


