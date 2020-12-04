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
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');
cfg.circus.channelname   = ft_getopt(cfg.circus, 'channelname', []);

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'SpikeTrials_Timelocked.mat']);
if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname, 'SpikeTrials');
    return;
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

    % to deal with multichannel data
    if isfield(SpikeRaw{ipart}, 'hdr')
        hdr = SpikeRaw{ipart}.hdr;
    else
        if isempty(cfg.circus.channelname)
            temp        = dir(fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], [cfg.prefix, 'p', num2str(ipart), '-multifile-*.ncs']));
        else
            temp        = dir(fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], cfg.circus.channelname{1}, [cfg.prefix, 'p', num2str(ipart), '-multifile-*.ncs']));
        end

        hdr_fname   = fullfile(temp(1).folder, temp(1).name);
        hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    end

    % create trials
    clear Trials
    for markername = string(cfg.name)

        SpikeTrials{ipart}.(markername) = []; %initialize to avoid bug in case no marker are found and need to output empty structure

        % clock time of each event
        clocktimes = [];
        for ifile = 1 : size(MuseStruct{ipart}, 2)
            if ~isfield(cfg.muse.startmarker, char(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{ifile}.markers, char(cfg.muse.startmarker.(markername)))
                continue
            end
            if ~isfield(MuseStruct{ipart}{ifile}.markers.(cfg.muse.startmarker.(markername)), 'clock')
                continue
            end
            clocktimes = [clocktimes, MuseStruct{ipart}{ifile}.markers.(cfg.muse.startmarker.(markername)).clock];
        end

        % first sample starts at 0, add cumulative directory along the way
        hyplabels_trl   = [];
        Startsample     = [];
        Endsample       = [];
        Offset          = [];
        Trialnr         = [];
        Trialnr_dir     = [];
        Filenr          = [];
        FileOffset      = [];
        Starttime       = [];
        Endtime         = [];
        dirOnset(1)     = 0;
        trialcount      = 1;

        % loop over directories
        for idir = 1 : size(MuseStruct{ipart}, 2)

            % calculate duration of single file/directory
            temp                = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{1}(1:end-2), '*.ncs']));
            hdrtemp             = ft_read_header(fullfile(temp(1).folder, temp(1).name));
            dirOnset(idir+1)    = dirOnset(idir) + hdrtemp.nSamples;

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
            trialcount_dir      = 1;

            % sort times in case markers have been e.g. combined
            [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(sidx);

            [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime      = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(sidx);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock         = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(sidx);

            % loop over events
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2)

                ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) * hdr.Fs);
                idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime * hdr.Fs) >= ss, 1, 'first');
                es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(idx) * hdr.Fs);

                if isempty(es)
                    continue
                end

                Startsample     = [Startsample; ss + cfg.spike.toi.(markername)(1) * hdr.Fs + dirOnset(idir)];
                Endsample       = [Endsample;   es + cfg.spike.toi.(markername)(2) * hdr.Fs + dirOnset(idir)];
                Starttime       = [Starttime;   MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent) + seconds(cfg.spike.toi.(markername)(1))];
                Endtime         = [Endtime;     MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(idx) + seconds(cfg.spike.toi.(markername)(2))];
                Offset          = [Offset;      cfg.spike.toi.(markername)(1) * hdr.Fs];
                Trialnr_dir     = [Trialnr_dir; trialcount_dir];
                Trialnr         = [Trialnr;     trialcount];
                Filenr          = [Filenr;      idir];
                FileOffset      = [FileOffset;  dirOnset(idir)];

                trialcount      = trialcount + 1;
                trialcount_dir  = trialcount_dir + 1;

                % find overlap with hypnogram markers
                trlstart        = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent);
                trlend          = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(ievent);
                overlap         = zeros(size(hyplabels));

                for ihyplabel = 1 : size(hyplabels, 2)

                    if ~isfield(MuseStruct{ipart}{idir}.markers, strcat(hyplabels{ihyplabel}, '__START__'))
                        continue
                    end

                    for ihyp = 1 : size(MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).synctime, 2)

                        hypstart = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).clock(ihyp);
                        hypend   = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__END__')).clock(ihyp);

                        % end of trial overlaps with beginning of sleepstate
                        if hypstart > trlstart && hypstart < trlend
                            overlap(ihyplabel) = seconds(trlend - hypstart);
                        end

                        % sleepstage falls fully within trial
                        if hypstart > trlstart && hypend < trlend
                            overlap(ihyplabel) = seconds(hypend - hypstart);
                        end

                        % beginning of trial overlaps with end of sleepstate
                        if hypstart < trlstart && hypend > trlstart
                            overlap(ihyplabel) = seconds(hypend - trlstart);
                        end

                        % trial falls fully within sleepstate
                        if hypstart < trlstart && hypend > trlend && ~(hypend > trlstart)
                            overlap(ihyplabel) = seconds(trlend - trlstart);
                        end
                    end
                end % ihyplabel

                % add sleep stage to trialinfo
                [val, indx] = max(overlap, [], 2);
                if val > 0
                    hyplabels_trl = [hyplabels_trl, hyplabels(indx)];
                else
                    hyplabels_trl = [hyplabels_trl, "NO_SCORE"];
                end

            end % ievent
        end %idir

        full_trial = Startsample > 0 & Endsample < hdr.nSamples;% so not to read before BOF or after EOFs
        cfgtemp                         = [];
        cfgtemp.trl                     = [Startsample, Endsample, Offset];
        cfgtemp.trl                     = cfgtemp.trl(full_trial, :); % so not to read before BOF or after EOFs

        if isempty(cfgtemp.trl)
            continue
        end

        % create spiketrials timelocked to events
        cfgtemp.trlunit                                         = 'samples';
        cfgtemp.hdr                                             = hdr;
        SpikeTrials{ipart}.(markername)                         = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
        if isfield(SpikeRaw{ipart},'clusternames'); SpikeTrials{ipart}.clusternames = SpikeRaw{ipart}.clusternames; end
        SpikeTrials{ipart}.(markername).trialinfo.clocktime     = clocktimes;
        SpikeTrials{ipart}.(markername).hdr                     = hdr;
        SpikeTrials{ipart}.(markername).trialinfo               = table;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_dir   = Trialnr_dir(full_trial); % trialnr which restarts for each dir
        SpikeTrials{ipart}.(markername).trialinfo.begsample     = Startsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endsample     = Endsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.offset        = Offset(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.trialduration = Endsample(full_trial)-Startsample(full_trial)+1;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr       = Trialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.idir          = Filenr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.fileoffset    = FileOffset(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.hyplabel      = hyplabels_trl(full_trial)';
        SpikeTrials{ipart}.(markername).trialinfo.starttime     = Starttime(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endtime       = Endtime(full_trial);

        %% Detect artefacts
        artefact = false(size(SpikeTrials{ipart}.(markername).trialinfo, 1), 1);
        artefact_length = zeros(size(SpikeTrials{ipart}.(markername).trialinfo, 1), 1);
        ft_progress('init','text')
        for ievent = 1 : size(SpikeTrials{ipart}.(markername).trialinfo, 1)
            ft_progress(ievent/size(SpikeTrials{ipart}.(markername).trialinfo, 1), 'Looking for overlap with artefacts in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.(markername).trialinfo, 1))
            trlstart = SpikeTrials{ipart}.(markername).trialinfo.starttime(ievent);
            trlend   = SpikeTrials{ipart}.(markername).trialinfo.endtime(ievent);

            for idir = 1 : size(MuseStruct{ipart}, 2)

                if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'synctime')
                    continue
                end

                for iart = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.clock, 2)
                    artstart = MuseStruct{ipart}{idir}.markers.BAD__START__.clock(iart);
                    artend   = MuseStruct{ipart}{idir}.markers.BAD__END__.clock(iart);

                    % full trial is before artefact
                    if trlstart < artstart && trlend < artstart
                        continue
                    % full trial is after artefact
                    elseif trlstart > artend && trlend > artend
                        continue
                    else
                        artefact(ievent) = true;
                        artefact_length(ievent) = seconds(artend - artstart) + artefact_length(ievent); %in case several artefacts intersect this trial
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
