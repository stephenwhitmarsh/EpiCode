function [cfg, MuseStruct_new] = PET_add_interIED_epochs(cfg, MuseStruct_orig, markername, toi, include_startend_of_file)

% Define periods between two successive markers
MuseStruct_new = MuseStruct_orig;
cfg.interIEDtrials_maxlength = ft_getopt(cfg, 'interIEDtrials_maxlength', Inf);

for ipart = 1 : size(MuseStruct_orig, 2)
    for idir = 1 : size(MuseStruct_orig{ipart}, 2)
        if ~isfield(MuseStruct_orig{ipart}{idir}, 'markers')
            continue
        end
        if ~isfield(MuseStruct_orig{ipart}{idir}.markers, markername)
            continue
        end
        if ~isfield(MuseStruct_orig{ipart}{idir}.markers.(markername), 'synctime')
            continue
        end
        if length(MuseStruct_orig{ipart}{idir}.markers.(markername).synctime) <= 1
            % need at least 2 events to define the period between 2 events
            continue
        end
           
        markerstart = sprintf('inter_%s__START__', markername);
        MuseStruct_new{ipart}{idir}.markers.(markerstart).clock    = MuseStruct_orig{ipart}{idir}.markers.(markername).clock(1:end-1);
        MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime = MuseStruct_orig{ipart}{idir}.markers.(markername).synctime(1:end-1);
        
        markerend = sprintf('inter_%s__END__', markername);
        MuseStruct_new{ipart}{idir}.markers.(markerend).clock      = MuseStruct_orig{ipart}{idir}.markers.(markername).clock(2:end);
        MuseStruct_new{ipart}{idir}.markers.(markerend).synctime   = MuseStruct_orig{ipart}{idir}.markers.(markername).synctime(2:end);
        
        if include_startend_of_file
            %include end of the file
            MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(end+1)    = MuseStruct_new{ipart}{idir}.markers.(markerend).clock(end);
            MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(end+1) = MuseStruct_new{ipart}{idir}.markers.(markerend).synctime(end);
            MuseStruct_new{ipart}{idir}.markers.(markerend).clock(end+1)    = MuseStruct_new{ipart}{idir}.endtime;
            MuseStruct_new{ipart}{idir}.markers.(markerend).synctime(end+1) = seconds(MuseStruct_new{ipart}{idir}.endtime - MuseStruct_new{ipart}{idir}.starttime);

            %include start of the file
            MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(end+1)    = MuseStruct_new{ipart}{idir}.starttime;
            MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(end+1) = 0;
            MuseStruct_new{ipart}{idir}.markers.(markerend).clock(end+1)    = MuseStruct_new{ipart}{idir}.markers.(markerstart).clock(1);
            MuseStruct_new{ipart}{idir}.markers.(markerend).synctime(end+1) = MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime(1);

            %sort the new events
            MuseStruct_new{ipart}{idir}.markers.(markerstart).clock    = sort(MuseStruct_new{ipart}{idir}.markers.(markerstart).clock);
            MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime    = sort(MuseStruct_new{ipart}{idir}.markers.(markerstart).synctime);
            MuseStruct_new{ipart}{idir}.markers.(markerend).clock    = sort(MuseStruct_new{ipart}{idir}.markers.(markerend).clock);
            MuseStruct_new{ipart}{idir}.markers.(markerend).synctime    = sort(MuseStruct_new{ipart}{idir}.markers.(markerend).synctime);
        end
    end
end

% add configuration to settings 

cfg.muse.startmarker.interIED  = markerstart;
cfg.muse.endmarker.interIED    = markerend;

cfg.epoch.toi.interIED         = toi;
cfg.epoch.pad.interIED         = 0;

if isfield(cfg, 'LFP')
    if isfield(cfg.LFP, 'name')
        cfg.LFP.toi.interIED           = toi;
        cfg.LFP.pad.interIED           = 0;      
    end
end

if isfield(cfg, 'FFT')
    if isfield(cfg.FFT, 'name')
        cfg.FFT.toi.interIED           = toi;
        cfg.FFT.pad.interIED           = 0;
    end
end

if isfield(cfg, 'TFR')
    if isfield(cfg.TFR, 'name')
        cfg.TFR.toi.interIED           = toi;
        cfg.TFR.pad.interIED           = 0;
    end
end

if isfield(cfg, 'spike')
    if isfield(cfg.spike, 'name')
        cfg.spike.toi.interIED         = toi;
        cfg.spike.pad.interIED         = 0;
    end
end