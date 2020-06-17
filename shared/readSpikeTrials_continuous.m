function [SpikeTrials] = readSpikeTrials_continuous(cfg,MuseStruct,SpikeRaw,force,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [SpikeTrials] = readSpikeTrials_continuous(cfg,MuseStruct,SpikeRaw,force,varargin)
% Cut a Fieldtrip raw spike structure into consecutive equal trials of 
% length defined in cfg. Make a trialinfo array similar to
% readSpikeTrials_MuseMarkers.m so the same analysis can be applied.
% Trials not complete at the end of the file is not created
%
% ### Necessary input:
%
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files
% cfg.circus.channel    = micro electrode names
% cfg.spike.triallength = length of the trials, in seconds. All the data
%                         will be cut in consecutive trials of this length.
%
% MuseStruct{ipart}     = info (e.g. events, files) of original data,
%                         used to  make a trialinfo array similar to the 
%                         one created in readSpikeTrials_MuseMarkers.m (eg 
%                         for later removing of artefacts based on Muse 
%                         Markers).
%
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% ### Output:
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% SpikeTrials           = spike data epoched in FieldTrip trial data structure
%
% Paul Baudin (paul.baudin@live.fr) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) || strcmp(varargin{1},'all')
    parts_to_read = 1:size(cfg.directorylist,2);
else
    parts_to_read = varargin{1};
end


fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeTrials_continuous.mat']);
if exist(fname,'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname,'SpikeTrials');
    return
else
    fprintf('(re-)computing SpikeTrials_continuous for %s\n', cfg.prefix(1:end-1));
end

for ipart = parts_to_read
    
    temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
    hdr_fname   = fullfile(temp(1).folder,temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
    
    % create trials
    clear Trials
    for ilabel = 1 : size(cfg.name,2)
        fprintf('for %s \n', cfg.name{ilabel});
        % create Fieldtrip trl by cutting data in trials of
        % cfg.spike.triallength length.
        %make a trialinfo similar as readSpikeTrials_MuseMarkers.m
        Startsample = [];
        Endsample   = [];
        Offset      = [];
        Trialnr     = [];
        Trialdir    = [];
        Filenr      = [];
        FileOffset  = [];
        clocktimes  = [];
        
        dirOnset(1) = 1;
        trialcount = 1;
        for idir = 1 : size(MuseStruct{ipart},2)
            clear ss es 
            
            temp                = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*', cfg.circus.channel{1}(1:end-2),'*.ncs']));
            hdrtemp             = ft_read_header(fullfile(temp(1).folder,temp(1).name));
            dirOnset(idir+1)    = dirOnset(idir) + hdrtemp.nSamples; % assuming all channels have same sampleinfo
            
            ss                      = dirOnset(idir) : hdrtemp.Fs * cfg.spike.triallength : dirOnset(idir) + hdrtemp.nSamples - hdrtemp.Fs * cfg.spike.triallength;%remove one cycle size so endsample is still in the file
            es                      = ss + hdrtemp.Fs * cfg.spike.triallength - 1;
            Startsample             = [Startsample, ss];
            Endsample               = [Endsample, es];
            Offset                  = [Offset, zeros(1,size(ss,2))];
            Trialdir                = [Trialdir; (1:size(ss,2))'];
            Trialnr                 = [Trialnr; (trialcount:trialcount+size(ss,2)-1)'];
            Filenr                  = [Filenr; ones(size(ss,2),1)*idir];
            FileOffset              = [FileOffset; ones(size(ss,2),1)*dirOnset(idir)];
            clocktimes              = [clocktimes, MuseStruct{ipart}{idir}.starttime : seconds(cfg.spike.triallength): MuseStruct{ipart}{idir}.starttime + seconds(cfg.spike.triallength)*(size(ss,2)-1)];
            
            trialcount      = trialcount + size(ss,2);
            
        end
        
        cfgtemp                         = [];
        cfgtemp.trl                     = [Startsample', Endsample', Offset'];
        cfgtemp.trl(:,4)                = ones(size(cfgtemp.trl,1),1) * idir;
        cfgtemp.trl(:,5)                = Trialdir;                % trialnr. to try to find trials that are missing afterwards
        cfgtemp.trl(:,6)                = Startsample';                          % startsample
        cfgtemp.trl(:,7)                = Endsample';                            % endsample
        cfgtemp.trl(:,8)                = Offset';                               % offset
        cfgtemp.trl(:,9)                = Endsample-Startsample+1;              % duration in samples
        cfgtemp.trl(:,10)               = Trialnr;                              % trialnr. to try to find trials that are missing afterwards
        cfgtemp.trl(:,11)               = Filenr;                               % trialnr. to try to find trials that are missing afterwards
        cfgtemp.trl(:,12)               = FileOffset;                           % trialnr. to try to find trials that are missing afterwards
               
        % create spiketrials 
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



