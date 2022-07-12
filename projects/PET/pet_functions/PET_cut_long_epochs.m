function [MuseStruct_new] = PET_cut_long_epochs(cfg, MuseStruct_orig, markername, max_epoch_length, new_epoch_length)

% Cut epochs into sub-epochs if they are too long.

MuseStruct_new = MuseStruct_orig;
markerstart    = cfg.muse.startmarker.(markername);
markerend      = cfg.muse.endmarker.(markername);

for ipart = 1 : size(MuseStruct_orig, 2)
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        if ~isfield(MuseStruct_orig{ipart}{idir}, 'markers')
            continue
        end
        if ~isfield(MuseStruct_orig{ipart}{idir}.markers, markerstart)
            continue
        end
        if ~isfield(MuseStruct_orig{ipart}{idir}.markers.(markerstart), 'synctime')
            continue
        end
        if length(MuseStruct_orig{ipart}{idir}.markers.(markerstart).synctime) <= 1
            % need at least 2 events to define the period between 2 events
            continue
        end
          
        % Cut too long epochs if required
        trial_length = MuseStruct_new{ipart}{idir}.markers.(markerend).synctime - MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime;
        to_cut = find(trial_length > max_epoch_length);
        if isempty(to_cut)
            continue
        end
        
        toadd.synctime.start = [];
        toadd.synctime.end   = [];
        toadd.clock.start    = datetime.empty;
        toadd.clock.end      = datetime.empty;
        
        for itrial = to_cut
            n_trials = floor(trial_length(itrial)/new_epoch_length);
            for i_addtrial = 0:n_trials-1
                toadd.synctime.start(end+1) = MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(itrial) +...
                    new_epoch_length * i_addtrial;
                toadd.synctime.end(end+1)   = MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(itrial) +...
                    new_epoch_length * (i_addtrial+1);
                toadd.clock.start(end+1)    = MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(itrial) +...
                    seconds(new_epoch_length * i_addtrial);
                toadd.clock.end(end+1)      = MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(itrial) +...
                    seconds(new_epoch_length * (i_addtrial+1));
            end
        end
        
        MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(to_cut) = [];
        MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(end+1 : end+size(toadd.synctime.start,2)) = toadd.synctime.start;
        MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime = sort(MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime);
               
        MuseStruct_new{ipart}{idir}.markers.(markerend).synctime(to_cut) = [];
        MuseStruct_new{ipart}{idir}.markers.(markerend).synctime(end+1 : end+size(toadd.synctime.start,2)) = toadd.synctime.end;
        MuseStruct_new{ipart}{idir}.markers.(markerend).synctime = sort(MuseStruct_new{ipart}{idir}.markers.(markerend).synctime);
               
        MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(to_cut) = [];
        MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(end+1 : end+size(toadd.synctime.start,2)) = toadd.clock.start;
        MuseStruct_new{ipart}{idir}.markers.(markerstart).clock = sort(MuseStruct_new{ipart}{idir}.markers.(markerstart).clock);
               
        MuseStruct_new{ipart}{idir}.markers.(markerend).clock(to_cut) = [];
        MuseStruct_new{ipart}{idir}.markers.(markerend).clock(end+1 : end+size(toadd.synctime.start,2)) = toadd.clock.end;
        MuseStruct_new{ipart}{idir}.markers.(markerend).clock = sort(MuseStruct_new{ipart}{idir}.markers.(markerend).clock);
       
    end
end
