function [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg, MuseStruct, SpikeRaw, force)

% [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg, MuseStruct, SpikeRaw, force)
% Make a Fieldtrip spike trials structure based on a Fieldtrip raw spike
% structure and on timings defined by Muse Markers.
% Trials separated between 2 files (ie begin on one file and end on the
% following file) are not created.
%
% ### Necessary input:
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.channel    = micro electrode names
% MuseStruct            = info (e.g. events, files) of original data,
%                         used to segment the spikes into trials
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% ### Optional cfg fields :
% cfg.circus.postfix     = string postfix appended to spike data results.
%                         Default = [].
% cfg.circus.part_list  = list of parts to analyse. Can be an array of
%                         integers, or 'all'. Default = 'all'.
%
% ### Output:
% SpikeTrials           = spike data epoched in FieldTrip trial data structure

% get the default cfg options
cfg.spike.postfix           = ft_getopt(cfg.spike,  'postfix', []);
cfg.spike.overlap           = ft_getopt(cfg.spike,  'overlap', []);
cfg.circus.part_list        = ft_getopt(cfg.circus, 'part_list', 'all');
cfg.circus.channelname      = ft_getopt(cfg.circus, 'channelname', []);

fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, "SpikeTrials_Timelocked", cfg.spike.postfix, ".mat"));

if nargin == 1
    if exist(fname, 'file')
        fprintf('Reading %s\n', fname);
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname, 'SpikeTrials');
            catch ME
                err_count = err_count + 1;
                disp('Something went wrong loading the file. Trying again...')                
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        SpikeTrials = {};
        return
    end
end

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    
    % repeat to deal with load errors
    count = 0;
    err_count = 0;
    while count == err_count
        try
            load(fname, 'SpikeTrials');
        catch ME
            err_count = err_count + 1;
            disp('Something went wrong loading the file. Trying again...')                
        end
        count = count + 1;
    end
    return
else
    fprintf('(re-)computing SpikeTrials_MuseMarkers for %s\n', cfg.prefix);
end

hyplabels = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"];

% start loop
for ipart = cfg.circus.part_list
    
    if ipart > size(SpikeRaw, 2)
        SpikeTrials{ipart} = [];
        continue
    end

    if isempty(SpikeRaw{ipart})
        SpikeTrials{ipart} = [];
        continue
    end

    if isfield(SpikeRaw{ipart}, 'hdr')
        hdr = SpikeRaw{ipart}.hdr;
    else
        disp('Can''t find header information, recovering it now');
        [~, ~, ~, hdr_in] = writeSpykingCircusFileList(cfg, false);
        hdr = hdr_in{ipart};
    end

    % create trials
    clear Trials
    for markername = string(cfg.spike.name)

        SpikeTrials{ipart}.(markername) = []; %initialize 

        hyplabels_trl   = [];
        Startsample     = int64([]);
        Endsample       = int64([]);
        Startsample_dir = int64([]);
        Endsample_dir   = int64([]);
        Startsec_dir    = [];
        Endsec_dir      = [];
        Offset          = int64([]);
        Trialnr         = [];
        Trialnr_dir     = [];
        Filenr          = [];
        FileOffset      = [];
        Starttime       = [];
        Endtime         = [];
        Directory       = [];
        TrialDirOnset   = int64([]);
        t0              = [];
        t1              = [];
        % first sample starts at 0, add cumulative directory along the way        
        dirOnset        = int64(0);
        trialcount      = 1;
        StartTrialnr    = [];
        EndTrialnr      = [];
        overlapcount    = [];
        
        for othermarkername = string(cfg.spike.overlap)
            overlapping.(othermarkername) = [];
        end
        
        % loop over directories
        for idir = 1 : size(MuseStruct{ipart}, 2)
            
            % make sure to do this before any 'continue'
            if idir > 1
                if isfield(cfg.circus, 'correct_chunk')
                    if cfg.circus.correct_chunk{ipart} %% FIXME: WORKAROUND FOR BUG IN SPYKING CIRCUS
                        dirOnset = dirOnset + hdr{idir-1}.nSamples - 512;
                    else
                        dirOnset = dirOnset + hdr{idir-1}.nSamples;
                    end
                else
                    dirOnset = dirOnset + hdr{idir-1}.nSamples;
                end
            end
            
            % Check whether data has events
            if ~isfield(cfg.muse.startmarker, char(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers, char(cfg.muse.startmarker.(markername)))
                fprintf('No events starting with %s found in filenr %d\n', cfg.muse.startmarker.(markername), idir);
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)), 'synctime')
                continue
            end
            if isempty(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2))
                continue
            end

            % starting trialcount fresh at every directory
            trialcount_dir = 1;
            
            % add trialnr before sorting for reference & checks
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr     = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr       = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime, 2);

            % sort times in case markers have been combined
            [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr     = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr(sidx);
            
            [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime      = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock         = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr       = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr(sidx);      

            % loop over events 
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2)

                ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) * hdr{idir}.Fs);
                es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(ievent) * hdr{idir}.Fs);
                
                % idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime * hdr{idir}.Fs) >= ss, 1, 'first');
                % idx = ievent;                
%                 if isempty(es); continue; end
                
                Startsample     = [Startsample;  int64(ss + (cfg.spike.toi.(markername)(1) - cfg.spike.pad.(markername)) * hdr{idir}.Fs) + dirOnset];
                Endsample       = [Endsample;    int64(es + (cfg.spike.toi.(markername)(2) + cfg.spike.pad.(markername)) * hdr{idir}.Fs) + dirOnset];
                Offset          = [Offset;       int64(     (cfg.spike.toi.(markername)(1) - cfg.spike.pad.(markername)) * hdr{idir}.Fs)];

                % extra info for trialinfo
                Startsample_dir = [Startsample_dir;  int64(ss + (cfg.spike.toi.(markername)(1) - cfg.spike.pad.(markername)) * hdr{idir}.Fs)];
                Endsample_dir   = [Endsample_dir;    int64(es + (cfg.spike.toi.(markername)(2) + cfg.spike.pad.(markername)) * hdr{idir}.Fs)];  
                Startsec_dir    = [Startsec_dir;  MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent)];
                Endsec_dir      = [Endsec_dir;    MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(ievent)];        
                t0              = [t0; ss];
                t1              = [t1; es];
                TrialDirOnset   = [TrialDirOnset; dirOnset];
                StartTrialnr    = [StartTrialnr; MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr(ievent)];
                EndTrialnr      = [EndTrialnr;   MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr(ievent)];
                Starttime       = [Starttime;    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent) + seconds(cfg.spike.toi.(markername)(1))];
                Endtime         = [Endtime;      MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(ievent)   + seconds(cfg.spike.toi.(markername)(2))];
                Trialnr_dir     = [Trialnr_dir;  trialcount_dir];
                Trialnr         = [Trialnr;      trialcount];
                Filenr          = [Filenr;       idir];
                FileOffset      = [FileOffset;   dirOnset];
                Directory       = [Directory;    cfg.directorylist{ipart}{idir}];           
                trialcount      = trialcount + 1;
                trialcount_dir  = trialcount_dir + 1;

                % will be used to find overlap between events
                trlstart        = Startsec_dir(ievent);
                trlend          = Endsec_dir(ievent);
                
                % find overlap with hypnogram markers
                overlap         = zeros(size(hyplabels));
                for ihyplabel = 1 : size(hyplabels, 2)

                    % if no hypnogram label is present continue
                    if ~isfield(MuseStruct{ipart}{idir}.markers, strcat(hyplabels{ihyplabel}, '__START__'))
                        continue
                    end
                    
                    for ihyp = 1 : size(MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).synctime, 2)
                        
                        hypstart = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).synctime(ihyp);
                        hypend   = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__END__')).synctime(ihyp);
                        
                        % end of trial overlaps with beginning of sleepstate
                        if hypstart > trlstart && hypstart < trlend
                            overlap(ievent, ihyplabel) = (trlend - hypstart);
                        end
                        
                        % sleepstage falls fully within trial
                        if hypstart > trlstart && hypend < trlend
                            overlap(ievent, ihyplabel) = (hypend - hypstart);
                        end
                        
                        % beginning of trial overlaps with end of sleepstate
                        if hypstart < trlstart && hypend > trlstart
                            overlap(ievent, ihyplabel) = (hypend - trlstart);
                        end
                        
                        % trial falls fully within sleepstate
                        if hypstart < trlstart && hypend > trlend && ~(hypend > trlstart)
                            overlap(ievent, ihyplabel) = (trlend - trlstart);
                        end
                    end % ihyp
                end % ihyplabel

                % add sleep stage to trialinfo
                [val, indx] = max(overlap, [], 2);
                if val > 0
                    hyplabels_trl = [hyplabels_trl, hyplabels(indx)];
                else
                    hyplabels_trl = [hyplabels_trl, "NO_SCORE"];
                end
                
                % overlap with other events
                for othermarkername = string(cfg.spike.overlap)
                    
                    overlap = 0;

                    % dont check overlap with itself
                    if strcmp(othermarkername, markername)
                        overlapping.(othermarkername) = [overlapping.(othermarkername), 0];
                        continue
                    end
                    
                    % if it doesn't exist
                    if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(othermarkername))
                        overlapping.(othermarkername) = [overlapping.(othermarkername), 0];                        
                        continue
                    end
                    
                    % if it doesn't have any events
                    if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(othermarkername)), 'synctime')
                        overlapping.(othermarkername) = [overlapping.(othermarkername), 0];                        
                        continue
                    end
                    
                    for iIED = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(othermarkername)).synctime, 2)
                        
                        IEDstart = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(othermarkername)).synctime(iIED); %% TODO: add toi 
                        IEDend   = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(othermarkername)).synctime(iIED);
                        
                        % other either before or after marker period
                        if ~(IEDend < trlstart || IEDstart > trlend)
                            overlap = overlap + 1;
                            fprintf('*** Found overlap with %s and %s, dir %d, trial %d \n', markername, othermarkername, idir, ievent);
                        end
       
                    end % iIED
                    
                    overlapping.(othermarkername) = [overlapping.(othermarkername), overlap];
                end % othermarkername
            end % ievent
        end %idir
        
        full_trial = Startsample > 0 & Endsample < SpikeRaw{ipart}.trialinfo.endsample;
        cfgtemp     = [];
        cfgtemp.trl = [Startsample, Endsample, Offset];
        cfgtemp.trl = cfgtemp.trl(full_trial, :); % so not to read before BOF or after EOF
        
        if isempty(cfgtemp.trl); continue; end

        % create spiketrials timelocked to events
        cfgtemp.trlunit                 = 'samples';
        cfgtemp.hdr                     = rmfield(hdr{1}, {'nSamples', 'orig'});
        cfgtemp.hdr.FirstTimeStamp      = 0;
        cfgtemp.hdr.TimeStampPerSample  = 1; 
        SpikeTrials{ipart}.(markername) = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
        
%         cfg_raster                  = [];
%         % cfgtemp.latency       = [cfg.stats.toi.(markername)(1), cfg.stats.toi.(markername)(2)];
%         cfg_raster.trialborders     = 'yes';
%         cfg_raster.linewidth        = 1;
%         cfg_raster.topplotsize      = 0.5;
%         cfg_raster.trialborders     = 'no';
%         ft_spike_plot_raster(cfg_raster, SpikeTrials{ipart}.(markername));
  
        if isfield(SpikeRaw{ipart},'clusternames'); SpikeTrials{ipart}.clusternames = SpikeRaw{ipart}.clusternames; end
        SpikeTrials{ipart}.(markername).hdr                     = cfgtemp.hdr;
        SpikeTrials{ipart}.(markername).trialinfo               = table;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_dir   = Trialnr_dir(full_trial); % trialnr which restarts for each dir
        SpikeTrials{ipart}.(markername).trialinfo.begsample     = Startsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endsample     = Endsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.begsample_dir = Startsample_dir(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endsample_dir = Endsample_dir(full_trial);  
        SpikeTrials{ipart}.(markername).trialinfo.offset        = Offset(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.trialduration = Endsample(full_trial)-Startsample(full_trial)+1;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr       = Trialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.idir          = Filenr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.fileoffset    = FileOffset(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.hyplabel      = hyplabels_trl(full_trial)';
        SpikeTrials{ipart}.(markername).trialinfo.starttime     = Starttime(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endtime       = Endtime(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.directory     = Directory(full_trial, :);
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_start = StartTrialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_end   = EndTrialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.t0            = t0(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.t1            = t1(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.overlapcount  = zeros(height(SpikeTrials{ipart}.(markername).trialinfo), 1);
        
        for othermarkername = string(cfg.spike.overlap)
            SpikeTrials{ipart}.(markername).trialinfo.(othermarkername) = overlapping.(othermarkername)(full_trial)';
            SpikeTrials{ipart}.(markername).trialinfo.overlapcount = SpikeTrials{ipart}.(markername).trialinfo.overlapcount + overlapping.(othermarkername)(full_trial)';
        end
        
        %% Detect artefacts
        artefact        = false(size(SpikeTrials{ipart}.(markername).trialinfo, 1), 1);
        artefact_length = zeros(size(SpikeTrials{ipart}.(markername).trialinfo, 1), 1);
        ft_progress('init','text')
        
        for ievent = 1 : size(SpikeTrials{ipart}.(markername).trialinfo, 1)
            ft_progress(ievent/size(SpikeTrials{ipart}.(markername).trialinfo, 1), 'Looking for overlap with artefacts in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.(markername).trialinfo, 1))
            trlstart = SpikeTrials{ipart}.(markername).trialinfo.begsample_dir(ievent);
            trlend   = SpikeTrials{ipart}.(markername).trialinfo.endsample_dir(ievent);

            for idir = 1 : size(MuseStruct{ipart}, 2)

                if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'synctime')
                    continue
                end

                for iart = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime, 2)
                    artstart = MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(iart);
                    artend   = MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(iart);

                    % full trial is before or after artefact
                    if trlend < artstart || trlstart > artend
                        continue
                    else
                        artefact(ievent) = true;
                        artefact_length(ievent) = (artend - artstart) + artefact_length(ievent); %in case several artefacts intersect this trial
                    end

                end % iart
            end % idir
        end % ievent
        
        ft_progress('close');

        % add artefact to trialinfo
        SpikeTrials{ipart}.(markername).trialinfo.artefact = artefact;
        SpikeTrials{ipart}.(markername).trialinfo.artefact_length = artefact_length;

    end % markername
end % ipart

save(fname, 'SpikeTrials', '-v7.3');
