function [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg,MuseStruct,SpikeRaw,force)

% [SpikeTrials] = readSpikeTrials_MuseMarkers(cfg,MuseStruct,SpikeRaw,force)
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
% 


% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');


if strcmp(cfg.circus.part_list,'all')
    cfg.circus.part_list = 1:size(cfg.directorylist,2);
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeTrials_MuseMarkers.mat']);
if exist(fname,'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname,'SpikeTrials');
    return;
else
    fprintf('(re-)computing SpikeTrials_MuseMarkers for %s\n', cfg.prefix);
end

for ipart = cfg.circus.part_list
    
    temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
    hdr_fname   = fullfile(temp(1).folder,temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % create trials
    clear Trials
    for markername = string(cfg.name)
       
        % clock time of each event
        clocktimes = [];
        for ifile = 1 : size(MuseStruct{ipart},2)
            if isfield(MuseStruct{ipart}{ifile}.markers,cfg.muse.startmarker.(markername))
                if isfield(MuseStruct{ipart}{ifile}.markers.(cfg.muse.startmarker.(markername)),'clock')
                    clocktimes = [clocktimes, MuseStruct{ipart}{ifile}.markers.(cfg.muse.startmarker.(markername)).clock];
                end
            end
        end
        
        % create Fieldtrip trl based on concatinated files by adding nr of
        % samples of each file
        Startsample = [];
        Endsample   = [];
        Offset      = [];
        Trialnr     = [];
        Trialnr_dir    = [];
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
            
            %Check than data has events
            if ~isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startmarker.(markername))
                fprintf('No events starting with %s found in filenr %d\n',cfg.muse.startmarker.(markername),idir);
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)),'synctime')
                continue
            end
            if isempty(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime,2))
                continue
            end
            trialcount_dir = 1;
            
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime,2)
                %end of trial : take the following marker, no need to have the same index has begining marker
                ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) * hdr.Fs);
                idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime * hdr.Fs) >= ss,1,'first');
                es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(idx) * hdr.Fs);
                
                if isempty(es)
                    continue
                end
                
                Startsample  = [Startsample; ss + cfg.epoch.toi.(markername)(1) * hdr.Fs + dirOnset(idir)];
                Endsample    = [Endsample;   es + cfg.epoch.toi.(markername)(2) * hdr.Fs + dirOnset(idir)];
                Offset       = [Offset; cfg.epoch.toi.(markername)(1) * hdr.Fs];
                Trialnr_dir  = [Trialnr_dir; trialcount_dir];
                Trialnr      = [Trialnr; trialcount];
                Filenr       = [Filenr; idir];
                FileOffset   = [FileOffset; dirOnset(idir)];
                trialcount   = trialcount + 1;
                trialcount_dir = trialcount_dir+1;
            end
        end %idir
        full_trial = Startsample > 0 & Endsample < hdr.nSamples;% so not to read before BOF or after EOFs
        cfgtemp                         = [];
        cfgtemp.trl                     = [Startsample, Endsample, Offset];
%         cfgtemp.trl(:,4)                = ones(size(cfgtemp.trl,1),1) * idir;
%         cfgtemp.trl(:,5)                = Trialdir;                % trialnr. restart for each dir
%         cfgtemp.trl(:,6)                = Startsample;                          % startsample
%         cfgtemp.trl(:,7)                = Endsample;                            % endsample
%         cfgtemp.trl(:,8)                = Offset;                               % offset
%         cfgtemp.trl(:,9)                = Endsample-Startsample+1;              % duration in samples
%         cfgtemp.trl(:,10)               = Trialnr;                              % Trial nr, find missing trials afterwards
%         cfgtemp.trl(:,11)               = Filenr;                               % nr of the nlx dir
%         cfgtemp.trl(:,12)               = FileOffset;                           % starting sample of the nlx dir
%         
        cfgtemp.trl                     = cfgtemp.trl(full_trial,:); % so not to read before BOF or after EOFs
        
        % create spiketrials timelocked to events
        cfgtemp.trlunit                          = 'samples';
        cfgtemp.hdr                              = hdr;
        SpikeTrials{ipart}.(markername)               = ft_spike_maketrials(cfgtemp,SpikeRaw{ipart});
        SpikeTrials{ipart}.(markername).clocktimes    = clocktimes;
        SpikeTrials{ipart}.(markername).hdr           = hdr;
        
        SpikeTrials{ipart}.(markername).trialinfo               = table;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr_dir   = Trialnr_dir(full_trial);% trialnr which restarts for each dir
        SpikeTrials{ipart}.(markername).trialinfo.begsample     = Startsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.endsample     = Endsample(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.offset        = Offset(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.trialduration = Endsample(full_trial)-Startsample(full_trial)+1;
        SpikeTrials{ipart}.(markername).trialinfo.trialnr       = Trialnr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.idir        = Filenr(full_trial);
        SpikeTrials{ipart}.(markername).trialinfo.fileoffset    = FileOffset(full_trial);
        
    end % markername
    
end % ipart

save(fname,'SpikeTrials');

end

