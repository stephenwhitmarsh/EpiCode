function MuseStruct = padHypnogram(MuseStruct)

for ipart = 1 : 3
    
    for idir = 1 : size(MuseStruct{ipart}, 2)
        if isfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__START__')
            
            if MuseStruct{ipart}{idir}.markers.AWAKE__END__.clock(1) ~= MuseStruct{ipart}{idir}.endtime
                MuseStruct{ipart}{idir}.markers.PRESLEEP__START__.synctime(1)     = 0;
                MuseStruct{ipart}{idir}.markers.PRESLEEP__START__.clock(1)        = MuseStruct{ipart}{idir}.starttime;
                
                MuseStruct{ipart}{idir}.markers.PRESLEEP__END__.synctime(1)       = MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(1);
                MuseStruct{ipart}{idir}.markers.PRESLEEP__END__.clock(1)          = MuseStruct{ipart}{idir}.markers.AWAKE__END__.clock(1);
                
                MuseStruct{ipart}{idir-1}.markers.PRESLEEP__START__.synctime(1)   = 7200 - MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(1);
                MuseStruct{ipart}{idir-1}.markers.PRESLEEP__START__.clock(1)      = MuseStruct{ipart}{idir}.starttime - seconds(7200 - MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(1));
                
                MuseStruct{ipart}{idir-1}.markers.PRESLEEP__END__.synctime(1)     = 7200;
                MuseStruct{ipart}{idir-1}.markers.PRESLEEP__END__.clock(1)        = MuseStruct{ipart}{idir-1}.endtime;
                
                fprintf('Added Presleep period to part %d, dir %d and %d\n', ipart, idir, idir-1);
                
                % replace
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(1)      = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(1)    = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.clock(1)         = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.clock(1)       = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.events           = MuseStruct{ipart}{idir}.markers.AWAKE__END__.events - 1;
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.events         = MuseStruct{ipart}{idir}.markers.AWAKE__START__.events - 1;
                if MuseStruct{ipart}{idir}.markers.AWAKE__END__.events == 0
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__END__');
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__START__');
                end
                
            else
                
                MuseStruct{ipart}{idir+1}.markers.PRESLEEP__START__.synctime(1)     = 0;
                MuseStruct{ipart}{idir+1}.markers.PRESLEEP__START__.clock(1)        = MuseStruct{ipart}{idir+1}.starttime;
                
                MuseStruct{ipart}{idir+1}.markers.PRESLEEP__END__.synctime(1)       = MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.synctime(1);
                MuseStruct{ipart}{idir+1}.markers.PRESLEEP__END__.clock(1)          = MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.clock(1);
                
                total = seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
                
                MuseStruct{ipart}{idir}.markers.PRESLEEP__START__.synctime(1)   = total - 7200 + MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.synctime(1);
                MuseStruct{ipart}{idir}.markers.PRESLEEP__START__.clock(1)      = MuseStruct{ipart}{idir}.endtime - seconds(7200) + seconds(MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.synctime(1));
                
                MuseStruct{ipart}{idir}.markers.PRESLEEP__END__.synctime(1)     = total;
                MuseStruct{ipart}{idir}.markers.PRESLEEP__END__.clock(1)        = MuseStruct{ipart}{idir}.endtime;
                
                fprintf('Added Presleep period to part %d, dir %d and %d\n', ipart, idir, idir+1);
                
                % replace
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(1)      = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(1)    = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.clock(1)         = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.clock(1)       = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.events           = MuseStruct{ipart}{idir}.markers.AWAKE__END__.events - 1;
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.events         = MuseStruct{ipart}{idir}.markers.AWAKE__START__.events - 1;
                if MuseStruct{ipart}{idir}.markers.AWAKE__END__.events == 0
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__END__');
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__START__');
                end
                MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.synctime(1)      = [];
                MuseStruct{ipart}{idir+1}.markers.AWAKE__START__.synctime(1)    = [];
                MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.clock(1)         = [];
                MuseStruct{ipart}{idir+1}.markers.AWAKE__START__.clock(1)       = [];
                MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.events           = MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.events - 1;
                MuseStruct{ipart}{idir+1}.markers.AWAKE__START__.events         = MuseStruct{ipart}{idir+1}.markers.AWAKE__START__.events - 1;
                if MuseStruct{ipart}{idir+1}.markers.AWAKE__END__.events == 0
                    MuseStruct{ipart}{idir+1}.markers = rmfield(MuseStruct{ipart}{idir+1}.markers, 'AWAKE__END__');
                    MuseStruct{ipart}{idir+1}.markers = rmfield(MuseStruct{ipart}{idir+1}.markers, 'AWAKE__START__');
                end
                
            end
            
            break
        end
    end
    
    remaining = 7200;
    
    for idir = size(MuseStruct{ipart}, 2) : -1 : 1
        if isfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__START__')
            if isfield(MuseStruct{ipart}{idir}.markers.AWAKE__START__, 'synctime')
                total = seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime);
                
                MuseStruct{ipart}{idir}.markers.POSTSLEEP__START__.synctime(1)     = MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(end);
                MuseStruct{ipart}{idir}.markers.POSTSLEEP__START__.clock(1)        = MuseStruct{ipart}{idir}.markers.AWAKE__START__.clock(end);
                
                MuseStruct{ipart}{idir}.markers.POSTSLEEP__END__.synctime(1)       = min(total, remaining);
                MuseStruct{ipart}{idir}.markers.POSTSLEEP__END__.clock(1)          = MuseStruct{ipart}{idir}.endtime;
                
                fprintf('Added Postsleep period (%0.f seconds) to part %d, dir %d\n', (total - MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(end)), ipart, idir);
                
                remaining = remaining - (total - MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(end));
                
                % replace
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.synctime(end)      = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.synctime(end)    = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.clock(end)         = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.clock(end)       = [];
                MuseStruct{ipart}{idir}.markers.AWAKE__END__.events             = MuseStruct{ipart}{idir}.markers.AWAKE__END__.events - 1;
                MuseStruct{ipart}{idir}.markers.AWAKE__START__.events           = MuseStruct{ipart}{idir}.markers.AWAKE__START__.events - 1;
                if MuseStruct{ipart}{idir}.markers.AWAKE__END__.events == 0
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__END__');
                    MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers, 'AWAKE__START__');
                end
                
                % if note yet required duration, add in subsequent dirs
                i = 1;
                while remaining > 0
                    total = seconds(MuseStruct{ipart}{idir+i}.endtime - MuseStruct{ipart}{idir+i}.starttime);
                    
                    MuseStruct{ipart}{idir+i}.markers.POSTSLEEP__START__.synctime(1)   = 0;
                    MuseStruct{ipart}{idir+i}.markers.POSTSLEEP__START__.clock(1)      = MuseStruct{ipart}{idir+i}.starttime;
                    
                    MuseStruct{ipart}{idir+i}.markers.POSTSLEEP__END__.synctime(1)     = min(total, remaining);
                    MuseStruct{ipart}{idir+i}.markers.POSTSLEEP__END__.clock(1)        = MuseStruct{ipart}{idir+i}.starttime + seconds(min(total, remaining));
                    
                    fprintf('%0.f seconds were remaining, I removed %0.f in part %d, dir %d\n', remaining, min(total, remaining), ipart, idir+i);
                    
                    remaining = remaining - min(total, remaining);
                    
                    i = i + 1;
                    
                end
            end
            
            break
        end
    end
    
    % remove "NO_SCORE"
    for idir = 1 : size(MuseStruct{ipart}, 2)
        if isfield(MuseStruct{ipart}{idir}.markers, 'NO_SCORE__START__')
            MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers,'NO_SCORE__START__');
            MuseStruct{ipart}{idir}.markers = rmfield(MuseStruct{ipart}{idir}.markers,'NO_SCORE__END__');
        end
    end
    
end
