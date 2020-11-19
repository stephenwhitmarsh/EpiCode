function [SpikeTrials] = readSpikeTrials_windowed(cfg, MuseStruct, SpikeRaw, force)

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

% FIXME Remove the need for MuseStruct

% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');

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
    
    temp        = dir(fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], [cfg.prefix, 'p', num2str(ipart), '-multifile-*.ncs']));
    hdr_fname   = fullfile(temp(1).folder, temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % create trials of x seconds
    cfgtemp                                         = [];
    cfgtemp.trl(:, 1)                               = 0 : hdr.Fs * (cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : hdr.nSamples - (hdr.Fs * cfg.hyp.spikewindow);
    cfgtemp.trl(:, 2)                               = cfgtemp.trl(:, 1) + hdr.Fs * cfg.hyp.spikewindow;
    cfgtemp.trl(:, 3)                               = zeros(size(cfgtemp.trl, 1), 1);
    cfgtemp.trlunit                                 = 'samples';
    cfgtemp.hdr                                     = hdr;
    SpikeTrials{ipart}.window                       = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
    SpikeTrials{ipart}.window.trialinfo             = table;
    SpikeTrials{ipart}.window.trialinfo.starttime   = (MuseStruct{ipart}{1}.starttime : seconds(cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : MuseStruct{ipart}{end}.endtime - seconds(cfg.hyp.spikewindow))';
    SpikeTrials{ipart}.window.trialinfo.endtime     = (MuseStruct{ipart}{1}.starttime + seconds(cfg.hyp.spikewindow) : seconds(cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : MuseStruct{ipart}{end}.endtime)';
    
    if height(SpikeTrials{ipart}.window.trialinfo) ~= size(cfgtemp.trl, 1)
        error('Difference in total duration of TRL and MuseMarker list - have you trimmed your directorylist?');
    end
    
    % find overlap with sleepstages
    overlap     = zeros(size(SpikeTrials{ipart}.window.trialinfo, 1), 6);
    hyplabels   = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"];
    
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        fprintf('Looking for overlap with hypnogram in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.window.trialinfo, 1))
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
        SpikeTrials{ipart}.window.trialinfo.(markername) = zeros(size(SpikeTrials{ipart}.window.trialinfo, 1), 1);
    end
    
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        trlstart = SpikeTrials{ipart}.window.trialinfo.starttime(ievent);
        trlend   = SpikeTrials{ipart}.window.trialinfo.endtime(ievent);
        
        for idir = 1 : size(MuseStruct{ipart}, 2)

            for markername = string(cfg.name)
                
                if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
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
                    % trial falls fully within sleepstate
                    if IEDstart > trlstart && IEDend < trlend
                        SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) = SpikeTrials{ipart}.window.trialinfo.(markername)(ievent) + 1;
                        fprintf('Found: %s, dir %d, trial %d \n', markername, idir, ievent);
                    end
                end % ihyp
            end % ihyplabel
        end % idir 
    end % ievent 
    
    artefact = false(size(SpikeTrials{ipart}.window.trialinfo, 1), 1);
    for ievent = 1 : size(SpikeTrials{ipart}.window.trialinfo, 1)
        fprintf('Looking for overlap with artefacts in trial %d of %d \n', ievent, size(SpikeTrials{ipart}.window.trialinfo, 1))
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
                    
                    % end of trial overlaps with beginning of sleepstate
                    if artstart > trlstart && artstart < trlend
                        artefact(ievent) = true;
                    end
                    
                    % sleepstage falls fully within trial
                    if artstart > trlstart && artend < trlend
                        artefact(ievent) = true;
                    end
                    
                    % beginning of trial overlaps with end of sleepstate
                    if artstart < trlstart && artend > trlstart
                        artefact(ievent) = true;
                    end
                    
                    % trial falls fully within sleepstate
                    if artstart < trlstart && artend > trlend && ~(hypend > trlstart)
                        artefact(ievent) = true;
                    end
                end % ihyp
        end % idir 
    end % ievent    
    
    % add artefact to trialinfo
    SpikeTrials{ipart}.window.trialinfo.artefact = artefact;

    
    
    
end % ipart

save(fname, 'SpikeTrials', '-v7.3');


