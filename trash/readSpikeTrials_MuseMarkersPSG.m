function [SpikeTrials] = readSpikeTrials_welch(cfg, MuseStruct, SpikeRaw, force)

% [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg, MuseStruct, SpikeRaw, force)
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

% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'SpikeTrials_MuseMarkers.mat']);
if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname, 'SpikeTrials');
    return;
else
    fprintf('(re-)computing SpikeTrials_MuseMarkers for %s\n', cfg.prefix);
end

for ipart = cfg.circus.part_list
    
    if ipart > size(SpikeRaw, 2)
        SpikeTrials{ipart}.welch = [];
        continue
    end
    
    if isempty(SpikeRaw{ipart})
        SpikeTrials{ipart}.welch = [];
        continue
    end
    
    temp        = dir(fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], [cfg.prefix, 'p', num2str(ipart), '-multifile-', cfg.circus.channel{1}, '*.ncs']));
    hdr_fname   = fullfile(temp(1).folder, temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data

    % create trials of x seconds
    cfgtemp                                         = [];
    cfgtemp.trl(:, 1)                               = 0 : hdr.Fs * (cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : hdr.nSamples - hdr.Fs * cfg.hyp.spikewindow;
    cfgtemp.trl(:, 2)                               = cfgtemp.trl(:, 1) + hdr.Fs * cfg.hyp.spikewindow;
    cfgtemp.trl(:, 3)                               = zeros(size(cfgtemp.trl, 1), 1);
    cfgtemp.trlunit                                 = 'samples';
    cfgtemp.hdr                                     = hdr;
    SpikeTrials{ipart}.welch                        = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
    SpikeTrials{ipart}.welch.trialinfo              = table;
    SpikeTrials{ipart}.welch.trialinfo.starttime    = (MuseStruct{ipart}{1}.starttime : seconds(cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : MuseStruct{ipart}{end}.endtime - seconds(cfg.hyp.spikewindow))';
    SpikeTrials{ipart}.welch.trialinfo.endtime      = (MuseStruct{ipart}{1}.starttime + seconds(cfg.hyp.spikewindow) : seconds(cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : MuseStruct{ipart}{end}.endtime)';
    
    if height(SpikeTrials{ipart}.welch.trialinfo) ~= size(cfgtemp.trl, 1)
        error('Difference in total duration of TRL and MuseMarker list - have you trimmed your directorylist?');
    end
    
    % find overlap with sleepstages
    overlap = zeros(size(SpikeTrials{ipart}.welch.trialinfo, 1), 6);
    for ievent = 1 : size(SpikeTrials{ipart}.welch.trialinfo, 1)
        trlstart = SpikeTrials{ipart}.welch.trialinfo.starttime(ievent);
        trlend   = SpikeTrials{ipart}.welch.trialinfo.endtime(ievent);
        for idir = 1 : size(MuseStruct{ipart}, 2)
            ihyplabel = 1;
            for hyplabel = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"]
                if ~isfield(MuseStruct{ipart}{idir}.markers, strcat(hyplabel, '__START__'))
                    continue
                end
                fprintf('Checking for overlap in part %d of %d, directory %d of %d, trial %d of %d with %s in\n', ipart, length(cfg.circus.part_list), idir, size(MuseStruct{ipart}, 2), ievent, size(SpikeTrials{ipart}.welch.trialinfo, 1), hyplabel);
                for ihyp = 1 : size(MuseStruct{ipart}{idir}.markers.(strcat(hyplabel, '__START__')).synctime, 2)
                    
                    hypstart = MuseStruct{ipart}{idir}.markers.(strcat(hyplabel, '__START__')).clock(ihyp);
                    hypend   = MuseStruct{ipart}{idir}.markers.(strcat(hyplabel, '__END__')).clock(ihyp);
                    
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
                end
                ihyplabel = ihyplabel + 1;
            end
        end
    end
    
    % add sleep stage to trialinfo    
    [~, indx] = max(overlap, [], 2);
    indx(indx == -1 | indx == 0 | indx == 6) = 0;
    SpikeTrials{ipart}.welch.trialinfo.stage  = indx;

end % ipart

save(fname, 'SpikeTrials');


