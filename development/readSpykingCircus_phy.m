function [SpikeRaw, SpikeTrials] = readSpykingCircus_phy(cfg,MuseStruct,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [SpikeRaw, SpikeTrials] = readSpykingCircus(cfg,MuseStruct,force)
% Read results from Spyking-Circus analysis; spike-times and templates.
%
%
% Necessary input: 
%
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files 
% cfg.circus.channel    = micro electrode names
%
% MuseStruct            = info (e.g. events, files) of original data,
%                         used to segment the spikes into trials
%
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
%
% Output:
%
% SpikeRaw = raw spike data in FieldTrip raw spike data structure
% SpikeTrials = spike data epoched in FieldTrip trial data structure
%
% (c) Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikedata.mat']);
if exist(fname,'file') && force == false
    load(fname,'SpikeRaw','SpikeTrials');
else
    
    %%%%% PHY ******
    
    d = dir(fullfile(cfg.circus.outputdir,[cfg.prefix,'*.GUI']));
    datadir = fullfile(d.folder,d.name);
    info = loadKSdir(datadir);
    temps = readNPY(fullfile(d.folder,d.name,'templates.npy'));

    SpikeRaw.template = temps;

    filename = fullfile(d.folder,d.name,'cluster_group.tsv'); 
    [cids, cgs] = readClusterGroupsCSV(filename);
    % cids is length nClusters, the cluster ID numbers
    % cgs is length nClusters, the "cluster group":
    % - 0 = noise
    % - 1 = mua
    % - 2 = good
    % - 3 = unsorted
    
    [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(fullfile(d.folder,d.name,'spike_times.npy'));
    
    spiketimes_indx     = readNPY(fullfile(d.folder,d.name,'spike_times.npy'));
    clusternr           = readNPY(fullfile(d.folder,d.name,'spike_templates.npy'))+1;
    
    temp        = dir(fullfile(cfg.datasavedir,[cfg.prefix,'all_data_*.ncs']));
    hdr_fname   = fullfile(temp(1).folder,temp(1).name);
    hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        
        for itemp = unique(clusternr)'
            % read spike timings (in samples)
            SpikeRaw.sample{itemp} = spiketimes_indx(clusternr==itemp);
            % convert to timestamps (for Fieldtrip when epoching raw data)            
            SpikeRaw.timestamp{itemp} = SpikeRaw.sample{itemp} * hdr.TimeStampPerSample + double(hdr.FirstTimeStamp);         

            [~,i] = max(mean(abs(SpikeRaw.template(itemp,:,:)),2));
            SpikeRaw.template_maxchan(itemp) = i;
            
            SpikeRaw.label{itemp} = sprintf('Unit %d',i);

        end
        
        for ilabel = 1 : size(cfg.name,2)
            
            % clock time of each event
            clocktimes = [];
            for ifile = 1 : size(MuseStruct,2)
                if isfield(MuseStruct{ifile}.markers,cfg.muse.startend{ilabel})
                    if isfield(MuseStruct{ifile}.markers.(cfg.muse.startend{ilabel}),'clock')
                        clocktimes = [clocktimes, MuseStruct{ifile}.markers.(cfg.muse.startend{ilabel}).clock];
                    end
                end
            end
            
            % create Fieldtrip trl based on concatinated files by adding nr of
            % samples of each file
            
            Startsample = [];
            Endsample = [];
            Offset = [];
            Trialnr = [];
            
            dirOnset = 0;
            trialcount = 1;
            for idir = 1 : size(MuseStruct,2)
                try
                    for ievent = 1 : size(MuseStruct{idir}.markers.(cfg.muse.startend{ilabel}).events,2)
                        Startsample  = [Startsample; MuseStruct{idir}.markers.(cfg.muse.startend{ilabel,1}).offset(ievent) + cfg.epoch.toi{ilabel}(1) * hdr.Fs + dirOnset];
                        Endsample    = [Endsample;   MuseStruct{idir}.markers.(cfg.muse.startend{ilabel,2}).offset(ievent) + cfg.epoch.toi{ilabel}(2) * hdr.Fs + dirOnset];
                        Offset       = [Offset; cfg.epoch.toi{ilabel}(1) * hdr.Fs];
                        Trialnr      = [Trialnr; trialcount];
                        trialcount   = trialcount + 1;
                    end
                catch
                    fprintf('No events starting with %s found in filenr %d, or more serious issue \n',cfg.muse.startend{ilabel},idir);
                end
                dirOnset = dirOnset + cfg.sampleinfo(idir,2);
            end
            
            cfgtemp                         = [];
            cfgtemp.trl                     = [Startsample, Endsample, Offset];
            cfgtemp.trl(:,4)                = ones(size(cfgtemp.trl,1),1) * idir;
            cfgtemp.trl(:,5)                = 1:size(cfgtemp.trl,1);                % trialnr. to try to find trials that are missing afterwards
            cfgtemp.trl(:,6)                = Startsample;                          % startsample
            cfgtemp.trl(:,7)                = Endsample;                            % endsample
            cfgtemp.trl(:,8)                = Offset;                               % offset
            cfgtemp.trl(:,9)                = Endsample-Startsample+1;              % duration in samples
            cfgtemp.trl(:,10)               = Trialnr;                              % trialnr. to try to find trials that are missing afterwards
            maxsamples                      = cumsum(cfg.sampleinfo(:,2));
            cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < maxsamples(end),:); % so not to read before BOF or after EOFs
            
            % create spiketrials timelocked to events
            cfgtemp.trlunit                 = 'samples';
            cfgtemp.hdr                     = hdr;
            SpikeTrials{ilabel}             = ft_spike_maketrials(cfgtemp,SpikeRaw);
            SpikeTrials{ilabel}.clocktimes  = clocktimes;
            
        end % patterns
        SpikeRaw.trialtime = [0 hdr.nSamples / hdr.Fs];
        
    end % if exists spike file.
    
    save(fname,'SpikeRaw','SpikeTrials');

end % save / force


