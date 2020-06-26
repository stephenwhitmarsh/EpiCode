function [SpikeTrials] = readSpikeTrials_continuous(cfg,SpikeRaw,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [SpikeTrials] = readSpikeTrials_continuous(cfg,SpikeRaw,force,varargin)
% Cut a Fieldtrip raw spike structure into consecutive equal trials of 
% length defined in cfg. Make a trialinfo array similar to
% readSpikeTrials_MuseMarkers.m so the same analysis can be applied.
% Note : last trial, if not complete at the end of the file, is not created
%
% ### Necessary input:
%
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.channel    = micro electrode names
% cfg.spike.triallength = length of the trials, in seconds. All the data
%                         will be cut in consecutive trials of this length.
%
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% force                 = whether to redo analyses or read previous save
%                         (true/false)
% 
% ### Optional cfg field
% cfg.circus.postfix     = string postfix appended to spike data results. 
%                         Default = []. 
% cfg.circus.part_list  = list of parts to analyse. Can be an array of
%                         integers, or 'all'. Default = 'all'.
%
% ### Output:
% SpikeTrials           = spike data epoched in FieldTrip trial data structure
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list    = ft_getopt(cfg.circus, 'part_list', 'all');

if strcmp(cfg.circus.part_list,'all')
    cfg.circus.part_list = 1:size(cfg.directorylist,2);
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeTrials_continuous.mat']);
if exist(fname,'file') && force == false
    fprintf('Loading %s\n', fname);
    load(fname,'SpikeTrials');
    return
else
    fprintf('(re-)computing SpikeTrials_continuous for %s\n', cfg.prefix(1:end-1));
end

for ipart = cfg.circus.part_list
    
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
        for idir = 1 : size(cfg.directorylist{ipart},2)
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



