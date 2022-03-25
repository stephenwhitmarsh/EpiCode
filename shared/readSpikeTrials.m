function [SpikeTrials] = readSpikeTrials(cfg, MuseStruct, SpikeRaw, force)

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
write                           = ft_getopt(cfg.spike,  'write', true);
cfg.spike.postfix               = ft_getopt(cfg.spike,  'postfix', []);
cfg.spike.overlap               = ft_getopt(cfg.spike,  'overlap', []);
cfg.circus.part_list            = ft_getopt(cfg.circus, 'part_list', 'all');
cfg.circus.channelname          = ft_getopt(cfg.circus, 'channelname', []);

% add markers to always look for overlap for
cfg.spike.overlap               = unique([cfg.spike.overlap, "BAD", "PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "PRESLEEP", "POSTSLEEP"], 'stable');
cfg.muse.startmarker.BAD        = 'BAD__START__';
cfg.muse.endmarker.BAD          = 'BAD__END__';
cfg.muse.startmarker.PHASE_1    = 'PHASE_1__START__';
cfg.muse.endmarker.PHASE_1      = 'PHASE_1__END__';
cfg.muse.startmarker.PHASE_2    = 'PHASE_2__START__';
cfg.muse.endmarker.PHASE_2      = 'PHASE_2__END__';
cfg.muse.startmarker.PHASE_3    = 'PHASE_3__START__';
cfg.muse.endmarker.PHASE_3      = 'PHASE_3__END__';
cfg.muse.startmarker.REM        = 'REM__START__';
cfg.muse.endmarker.REM          = 'REM__END__';
cfg.muse.startmarker.AWAKE      = 'AWAKE__START__';
cfg.muse.endmarker.AWAKE        = 'AWAKE__END__';
cfg.muse.startmarker.PRESLEEP   = 'PRESLEEP__START__';
cfg.muse.endmarker.PRESLEEP     = 'PRESLEEP__END__';
cfg.muse.startmarker.POSTSLEEP  = 'POSTSLEEP__START__';
cfg.muse.endmarker.POSTSLEEP    = 'POSTSLEEP__END__';

% used later to determine sleepstage
hyplabels                       = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "PRESLEEP", "POSTSLEEP"];
hypindex                        = false(length(cfg.spike.overlap), 1);
for i = 1 : length(cfg.spike.overlap)
    if any(strcmp(cfg.spike.overlap{i}, hyplabels))
        hypindex(i) = true;
    end
end

if nargin == 1
    for markername = string(cfg.spike.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'SpikeTrials_', markername, cfg.spike.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                    for ipart = 1 : size(temp.SpikeTrials, 2)
                        SpikeTrials{ipart}.(markername) = temp.SpikeTrials{ipart}.(markername);
                    end
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
        else
            fprintf('Will be (re-) computing SpikeTrial data for %s\n', markername);
        end
    end
    return
elseif ~force
    missing = [];
    for markername = string(cfg.spike.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'SpikeTrials_', markername, cfg.spike.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                    for ipart = 1 : size(temp.SpikeTrials, 2)
                        SpikeTrials{ipart}.(markername) = temp.SpikeTrials{ipart}.(markername);
                    end
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
        else
            fprintf('Will be (re-) computing SpikeTrial data for %s\n', markername);
            missing = [missing, markername];
        end
    end
    cfg.spike.name = missing;
end

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

% loop over markers
for markername = string(cfg.spike.name)
    
    % start loop
    for ipart = cfg.circus.part_list
        
        SpikeTrials{ipart}.(markername) = []; % initialize

        if ipart > size(SpikeRaw, 2)
            continue
        end
        
        if isempty(SpikeRaw{ipart})
            continue
        end
        
        if isfield(SpikeRaw{ipart}, 'hdr')
            hdr = SpikeRaw{ipart}.hdr;
        else
            disp('Can''t find header information, recovering it now');
            [~, ~, ~, hdr_in] = writeSpykingCircusFileList(cfg, false);
            hdr = hdr_in{ipart};
        end

        itrial          = 1;  
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
        
        overlap_sec     = [];
        overlap_cnt     = [];
        
        % loop over directories
        for idir = 1 : size(MuseStruct{ipart}, 2)
            
            fprintf('Processing directory %d from %d\n', idir, size(MuseStruct{ipart}, 2))
            
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
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime)
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
            
            overlap_sec = [overlap_sec; zeros(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2), size(cfg.spike.overlap, 2))];
            overlap_cnt = [overlap_cnt; zeros(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2), size(cfg.spike.overlap, 2))];
            
            ft_progress('init', 'text', 'Please wait...')
            
            % loop over events
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2)
                
                ft_progress(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2) / ievent, ...
                    'Processing event %d from %d', ievent, size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2))
                
                ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) * hdr{idir}.Fs);
                es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(ievent) * hdr{idir}.Fs);
                
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
                EndTrialnr      = [EndTrialnr;   MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr(ievent)];
                Starttime       = [Starttime;    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent) + seconds(cfg.spike.toi.(markername)(1))];
                Endtime         = [Endtime;      MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(ievent)   + seconds(cfg.spike.toi.(markername)(2))];
                Trialnr_dir     = [Trialnr_dir;  trialcount_dir];
                Trialnr         = [Trialnr;      trialcount];
                Filenr          = [Filenr;       idir];
                Directory       = [Directory;    cfg.directorylist{ipart}{idir}];
                trialcount      = trialcount + 1;
                trialcount_dir  = trialcount_dir + 1;
                
                % will be used to find overlap between events
                trlstart                = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) + cfg.epoch.toi.(markername)(1);
                trlend                  = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(ievent) + cfg.epoch.toi.(markername)(2);      

                % find overlap
                for iother = 1 : size(cfg.spike.overlap, 2)
                    
                    % if the event is not present continue
                    if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(cfg.spike.overlap{iother}))
                        continue
                    end
                    if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.spike.overlap{iother})), 'synctime')
                        continue
                    end
                    if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.spike.overlap{iother})).synctime)
                        continue
                    end
                    
                    % loop over instances of overlap-event
                    for i = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.spike.overlap{iother})).synctime, 2)
                        
                        other_start = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.spike.overlap{iother})).synctime(i);
                        other_end   = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(cfg.spike.overlap{iother})).synctime(i);
                        
                        % trial falls fully within other
                        if other_start < trlstart && other_end > trlend
                            overlap_sec(itrial, iother) = overlap_sec(itrial, iother) + (trlend - trlstart);
                            overlap_cnt(itrial, iother) = overlap_cnt(itrial, iother) + 1;
                            % other falls fully within trial
                        elseif other_start > trlstart && other_end < trlend
                            overlap_sec(itrial, iother) = overlap_sec(itrial, iother) + (other_end - other_start);
                            overlap_cnt(itrial, iother) = overlap_cnt(itrial, iother) + 1;
                            % trial overlaps with beginning of other
                        elseif other_start > trlstart && other_start < trlend
                            overlap_sec(itrial, iother) = overlap_sec(itrial, iother) + (trlend - other_start);
                            overlap_cnt(itrial, iother) = overlap_cnt(itrial, iother) + 1;
                            % trial overlaps with end of other
                        elseif other_start < trlstart && other_end > trlstart
                            overlap_sec(itrial, iother) = overlap_sec(itrial, iother) + (other_end - trlstart);
                            overlap_cnt(itrial, iother) = overlap_cnt(itrial, iother) + 1;
                        end
                        
                    end % i
                end % iother
                
                % add sleep stage to trialinfo
                sel = overlap_sec(itrial, hypindex);
                [val, indx] = max(sel);
                if val > 0
                    hyplabels_trl(itrial) = hyplabels(indx);
                else
                    hyplabels_trl(itrial) = "NO_SCORE";
                end
                
                itrial = itrial + 1;
                
            end % ievent
            
            ft_progress('close');
            
        end %idir
        
        % create Fieldtrip trl
        full_trial  = Startsample > 0 & Endsample < SpikeRaw{ipart}.trialinfo.endsample;
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
        SpikeTrials{ipart}.(markername).trialinfo.hyplabel      = hyplabels_trl(full_trial)';
        SpikeTrials{ipart}.(markername).trialinfo.starttime     = Starttime(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endtime       = Endtime(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.directory     = Directory(full_trial, :);
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_start = StartTrialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_end   = EndTrialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.t0            = t0(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.t1            = t1(full_trial);
        
        for ioverlap = 1 : size(cfg.spike.overlap, 2)
            othermarkername = cfg.spike.overlap{ioverlap};
            SpikeTrials{ipart}.(markername).trialinfo.([othermarkername '_sec']) = overlap_sec(full_trial, ioverlap);
            SpikeTrials{ipart}.(markername).trialinfo.([othermarkername '_cnt']) = overlap_cnt(full_trial, ioverlap);
        end
        
    end % ipart
    
    if write
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'SpikeTrials_', markername, cfg.spike.postfix, '.mat'));
        fprintf('Saving SpikeTrial data for %s\n', markername);
        saveMarker_SpikeTrials(SpikeTrials, markername, fname)
    end  
    
end % markername 