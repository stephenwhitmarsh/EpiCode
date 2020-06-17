function [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg,MuseStruct,SpikeRaw,force,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg,MuseStruct,SpikeRaw,force,varargin)
% Make a Fieldtrip spike trials structure based on a Fieldtrip raw spike 
% structure and on timings defined by Muse Markers.
% Trials separated between 2 files (ie begin on one file and end on the
% following file) are not created.
%
% ### Necessary input:
%
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files
% cfg.circus.channel    = micro electrode names
%
% MuseStruct{ipart}     = info (e.g. events, files) of original data,
%                         used to segment the spikes into trials
%
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% ### Output:
%
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% SpikeTrials           = spike data epoched in FieldTrip trial data structure
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
% Modified by Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) || strcmp(varargin{1},'all')
    parts_to_read = 1:size(cfg.directorylist,2);
else
    parts_to_read = varargin{1};
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeTrials_MuseMarkers.mat']);
if exist(fname,'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname,'SpikeTrials');
    return;
else
    fprintf('(re-)computing SpikeTrials_MuseMarkers for %s\n', cfg.prefix);
end

for ipart = parts_to_read
    
    temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
    hdr_fname   = fullfile(temp(1).folder,temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % create trials
    clear Trials
    for ilabel = 1 : size(cfg.name,2)
        
        % clock time of each event
        clocktimes = [];
        for ifile = 1 : size(MuseStruct{ipart},2)
            if isfield(MuseStruct{ipart}{ifile}.markers,cfg.muse.startend{ilabel})
                if isfield(MuseStruct{ipart}{ifile}.markers.(cfg.muse.startend{ilabel}),'clock')
                    clocktimes = [clocktimes, MuseStruct{ipart}{ifile}.markers.(cfg.muse.startend{ilabel}).clock];
                end
            end
        end
        
        % create Fieldtrip trl based on concatinated files by adding nr of
        % samples of each file
        Startsample = [];
        Endsample   = [];
        Offset      = [];
        Trialnr     = [];
        Trialdir    = [];
        Filenr      = [];
        FileOffset  = [];
        
        dirOnset(1) = 0;
        trialcount = 1;
        for idir = 1 : size(MuseStruct{ipart},2)
            %compute dir onset of the next dir based of the length
            %of this dir. Do it before searching for events to avoid
            %error if one dir has no event
            temp        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*', cfg.circus.channel{1}(1:end-2),'*.ncs']));
            hdrtemp     = ft_read_header(fullfile(temp(1).folder,temp(1).name));
            dirOnset(idir+1)    = dirOnset(idir) + hdrtemp.nSamples; % assuming all channels have same sampleinfo
            
            if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{ilabel})
                if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}),'synctime')
                    if ~isempty(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}).synctime,2))
                        trialcount_dir = 1;
                        
                        for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}).synctime,2)
                            %end of trial : take the following marker, no need to have the same index has begining marker
                            ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime(ievent) * hdr.Fs);
                            idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime * hdr.Fs) >= ss,1,'first');
                            es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime(idx) * hdr.Fs);
                            
                            if ~isempty(es)
                                Startsample  = [Startsample; ss + cfg.epoch.toi{ilabel}(1) * hdr.Fs + dirOnset(idir)];
                                Endsample    = [Endsample;   es + cfg.epoch.toi{ilabel}(2) * hdr.Fs + dirOnset(idir)];
                                Offset       = [Offset; cfg.epoch.toi{ilabel}(1) * hdr.Fs];
                                Trialdir     = [Trialdir; trialcount_dir];
                                Trialnr      = [Trialnr; trialcount];
                                Filenr       = [Filenr; idir];
                                FileOffset   = [FileOffset; dirOnset(idir)];
                                trialcount   = trialcount + 1;
                                trialcount_dir = trialcount_dir+1;
                            end
                        end
                    else
                        fprintf('No events starting with %s found in filenr %d\n',cfg.muse.startend{ilabel},idir);
                    end
                end
            end
            
        end
        
        cfgtemp                         = [];
        cfgtemp.trl                     = [Startsample, Endsample, Offset];
        cfgtemp.trl(:,4)                = ones(size(cfgtemp.trl,1),1) * idir;
        cfgtemp.trl(:,5)                = Trialdir;                % trialnr. restart for each dir
        cfgtemp.trl(:,6)                = Startsample;                          % startsample
        cfgtemp.trl(:,7)                = Endsample;                            % endsample
        cfgtemp.trl(:,8)                = Offset;                               % offset
        cfgtemp.trl(:,9)                = Endsample-Startsample+1;              % duration in samples
        cfgtemp.trl(:,10)               = Trialnr;                              % 
        cfgtemp.trl(:,11)               = Filenr;                               %
        cfgtemp.trl(:,12)               = FileOffset;                           % 
        
        cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs
        
        % create spiketrials timelocked to events
        cfgtemp.trlunit                          = 'samples';
        cfgtemp.hdr                              = hdr;
        SpikeTrials{ipart}{ilabel}               = ft_spike_maketrials(cfgtemp,SpikeRaw{ipart});
        SpikeTrials{ipart}{ilabel}.clocktimes    = clocktimes;
        SpikeTrials{ipart}{ilabel}.hdr           = hdr;
        SpikeTrials{ipart}{ilabel}.analysis_name = cfg.name{ilabel};
        
    end % ilabel
    
end % ipart

save(fname,'SpikeTrials');

end

