function [SpikeTrials] = readSpikeTrials_windowed(cfg, MuseStruct, SpikeRaw, force)

% [SpikeTrials] = readSpikeTrials_windowed(cfg, MuseStruct, SpikeRaw, force)
% Make a Fieldtrip spike trials structure based on equal overlapping
% windows.
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
%

% FIXME Remove the need for MuseStruct

% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');
cfg.circus.channelname   = ft_getopt(cfg.circus, 'channelname', []);

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'SpikeTrials_Windowed.mat']);
if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname, 'SpikeTrials');
    return;
else
    fprintf('(re-)computing SpikeTrials_MuseMarkers for %s\n', cfg.prefix);
end

for ipart = cfg.circus.part_list

    if ipart > size(SpikeRaw, 2)
        SpikeTrials{ipart}.window = [];
        continue
    end

    if isempty(SpikeRaw{ipart})
        SpikeTrials{ipart}.window = [];
        continue
    end

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

    % create trials of x seconds
    cfgtemp                                         = [];
    cfgtemp.trl(:, 1)                               = 0 : hdr.Fs * (cfg.spikewin.windowsize - cfg.spikewin.windowsize * cfg.spikewin.windowoverlap) : hdr.nSamples - (hdr.Fs * cfg.spikewin.windowsize);
    cfgtemp.trl(:, 2)                               = cfgtemp.trl(:, 1) + hdr.Fs * cfg.spikewin.windowsize - 1; %-1 because starts at zero
    cfgtemp.trl(:, 3)                               = zeros(size(cfgtemp.trl, 1), 1);
    cfgtemp.trlunit                                 = 'samples';
    cfgtemp.hdr                                     = hdr;
    SpikeTrials{ipart}.window                       = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
    if isfield(SpikeRaw{ipart},'clustername'); SpikeTrials{ipart}.clustername = SpikeRaw{ipart}.clustername; end
    SpikeTrials{ipart}.window.trialinfo             = table;
    SpikeTrials{ipart}.window.trialinfo.starttime   = (MuseStruct{ipart}{1}.starttime : seconds(cfg.spikewin.windowsize - cfg.spikewin.windowsize * cfg.spikewin.windowoverlap) : MuseStruct{ipart}{end}.endtime - seconds(cfg.spikewin.windowsize))';
    SpikeTrials{ipart}.window.trialinfo.endtime     = (MuseStruct{ipart}{1}.starttime + seconds(cfg.spikewin.windowsize -1/hdr.Fs) : seconds(cfg.spikewin.windowsize - cfg.spikewin.windowsize * cfg.spikewin.windowoverlap) : MuseStruct{ipart}{end}.endtime)';

    if height(SpikeTrials{ipart}.window.trialinfo) ~= size(cfgtemp.trl, 1)
        error('Difference in total duration of TRL and MuseMarker list - have you trimmed your directorylist?');
    end

    % find overlap with sleepstages
    overlap     = zeros(size(SpikeTrials{ipart}.window.trialinfo, 1), 6);
    hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"];

    ft_progress('init','text')
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        ft_progress(ievent/size(SpikeTrials{ipart}.window.trialinfo, 1), 'Looking for overlap with hypnogram in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.window.trialinfo, 1))
        trlstart = SpikeTrials{ipart}.window.trialinfo.starttime(ievent);
        trlend   = SpikeTrials{ipart}.window.trialinfo.endtime(ievent);

        for idir = 1 : size(MuseStruct{ipart}, 2)

            for ihyplabel = 1 : size(hyplabels, 2)

                if ~isfield(MuseStruct{ipart}{idir}.markers, strcat(hyplabels{ihyplabel}, '__START__'))
                    continue
                end

                for ihyp = 1 : size(MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).synctime, 2)

                    hypstart = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__START__')).clock(ihyp);
                    hypend   = MuseStruct{ipart}{idir}.markers.(strcat(hyplabels{ihyplabel}, '__END__')).clock(ihyp);

                    % end of trial overlaps with beginning of sleepstate
                    if hypstart > trlstart && hypstart < trlend
                        overlap(ievent, ihyplabel) = seconds(trlend - hypstart);
                    end

                    % sleepstage falls fully within trial
                    if hypstart > trlstart && hypend < trlend
                        overlap(ievent, ihyplabel) = seconds(hypend - hypstart);
                    end

                    % beginning of trial overlaps with end of sleepstate
                    if hypstart < trlstart && hypend > trlstart
                        overlap(ievent, ihyplabel) = seconds(hypend - trlstart);
                    end

                    % trial falls fully within sleepstate
                    if hypstart < trlstart && hypend > trlend && ~(hypend > trlstart)
                        overlap(ievent, ihyplabel) = seconds(trlend - trlstart);
                    end
                end % ihyp
            end % ihyplabel
        end % idir
    end % ievent
    ft_progress('close');

    % add sleep stage to trialinfo
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        [val, indx] = max(overlap(ievent, :));
        if val > 0
            SpikeTrials{ipart}.window.trialinfo.hyplabel(ievent) = string(hyplabels{indx});
        else
            SpikeTrials{ipart}.window.trialinfo.hyplabel(ievent) = "NO_SCORE";
        end
    end

    % find overlap with IEDs
    for markername = string(cfg.name)
        SpikeTrials{ipart}.window.trialinfo.(char(markername)) = zeros(size(SpikeTrials{ipart}.window.trialinfo, 1), 1);
    end

    ft_progress('init','text');
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        trlstart = SpikeTrials{ipart}.window.trialinfo.starttime(ievent);
        trlend   = SpikeTrials{ipart}.window.trialinfo.endtime(ievent);

        for idir = 1 : size(MuseStruct{ipart}, 2)

            for markername = string(cfg.name)

                if ~isfield(MuseStruct{ipart}{idir}.markers, char(markername))
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers.(char(markername)), 'synctime')
                    continue
                end

                for iIED = 1 : size(MuseStruct{ipart}{idir}.markers.(markername).synctime, 2)

                    IEDstart = MuseStruct{ipart}{idir}.markers.(markername).clock(iIED);
                    IEDend   = MuseStruct{ipart}{idir}.markers.(markername).clock(iIED);

%                     % end of trial only overlaps with beginning of sleepstate
%                     if IEDstart > trlstart && IEDstart < trlend
%                         SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) = SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) + 1;
%                     end
%
%                     % sleepstage falls fully within trial
%                     if IEDstart > trlstart && IEDend < trlend
%                         SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) = SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) + 1;
%                     end
%
%                     % beginning of trial only overlaps with end of sleepstate
%                     if IEDstart < trlstart && IEDend > trlstart
%                         SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) = SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) + 1;
%                     end
%
                    % ied falls fully within trial
                    if IEDstart > trlstart && IEDend < trlend
                        SpikeTrials{ipart}.window.trialinfo.(char(markername))(ievent) = SpikeTrials{ipart}.window.trialinfo.(char(markername))(ievent) + 1;
                        ft_progress(0,'Found: %s, dir %d, trial %d \n', markername, idir, ievent);
                    end
                end % ihyp
            end % ihyplabel
        end % idir
    end % ievent
    ft_progress('close');

    artefact        = false(size(SpikeTrials{ipart}.window.trialinfo, 1), 1);
    artefact_length = zeros(size(SpikeTrials{ipart}.window.trialinfo, 1), 1);
    ft_progress('init','text')
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        ft_progress(ievent/size(SpikeTrials{ipart}.window.trialinfo, 1), 'Looking for overlap with artefacts in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.window.trialinfo, 1))
        trlstart = SpikeTrials{ipart}.window.trialinfo.starttime(ievent);
        trlend   = SpikeTrials{ipart}.window.trialinfo.endtime(ievent);

        for idir = 1 : size(MuseStruct{ipart}, 2)

                if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'synctime')
                    continue
                end

                for iart = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime, 2)

                    artstart = MuseStruct{ipart}{idir}.markers.BAD__START__.clock(iart);
                    artend   = MuseStruct{ipart}{idir}.markers.BAD__END__.clock(iart);

                    %full trial is before artefact
                    if trlstart < artstart && trlend < artstart
                        continue
                        %full trial is after artefact
                    elseif trlstart > artend && trlend > artend
                        continue
                    else
                        artefact(ievent) = true;
                        artefact_length(ievent) = seconds(artend - artstart) + artefact_length(ievent);
                    end


                end % ihyp
        end % idir
    end % ievent
    ft_progress('close');

    % add artefact to trialinfo
    SpikeTrials{ipart}.window.trialinfo.artefact        = artefact;
    SpikeTrials{ipart}.window.trialinfo.artefact_length = artefact_length;

end % ipart

save(fname, 'SpikeTrials', '-v7.3');
