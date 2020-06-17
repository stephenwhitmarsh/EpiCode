function [SpikeWaveforms] = readSpikeWaveforms(cfg,SpikeTrials,force,parts_to_read)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [SpikeWaveforms] = readSpikeWaveforms(cfg,SpikeTrials,force,parts_to_read)
% Cut continuous data into trials for each spike sample. 
%
% ### INPUT
% cfg.name                      = label of the analysis
% cfg.prefix                    = prefix to output files
% cfg.datasavedir               = data directory of results
% cfg.circus.channel            = channel names which were analyzed by Spyking-Circus
% cfg.spikewaveform.toi         = in seconds, time to load before and after the peak of the spike
% cfg.spikewaveform.cutoff      = high pass filter frequency to apply to raw data
% cfg.spikewaveform.nspikes     = maximum number of spike waveforms to load. Can be 'all'. If there are more 
%                                 spikes than nspikes, a random selection of spikes is done.
% force                         = whether to redo analyses or read previous save (true/false)
%
% ## OUTPUT 
% SpikeWaveforms                = Fieldtrip raw structure. One cell per
%                                 unit, one trial per spike.
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(parts_to_read,'all')
    parts_to_read = 1:size(cfg.directorylist,2);
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spike_waveform.mat']);

if exist(fname,'file') && force == false
    fprintf('*******************************************\n');
    fprintf('*** Loading precomputed spike waveforms ***\n');
    fprintf('*******************************************\n\n');
    
    load(fname,'SpikeWaveforms');
    return;
    
elseif exist(fname,'file') && force == true
    fprintf('*******************************************\n');
    fprintf('*** Forced  recomputing spike waveforms ***\n');
    fprintf('*******************************************\n\n');
    
else
    fprintf('*********************************\n');
    fprintf('*** Computing spike waveforms ***\n');
    fprintf('*********************************\n\n');
end


for ipart = parts_to_read
    
    for ilabel = 1:size(cfg.name,2)
        
        for ichan = 1:size(cfg.circus.channel,2)
            
            % find concatenate channel used by Spyking Circus
            temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{ichan},'.ncs']));
            if ~isempty(temp)
                
                datafile = fullfile(temp.folder,temp.name);
                fprintf('Found datafile: %s\n',datafile);
                
                hdr = ft_read_header(datafile);
                
                clusters_idx = [];
                clusters_idx = find(SpikeTrials{ipart}{ilabel}.template_maxchan == ichan - 1); %-1 because Phy starts counting at zero
                
                if ~isempty(clusters_idx)
                    for icluster = clusters_idx
                        
                        
                        %Select random spike if required
                        if strcmp(cfg.spikewaveform.nspikes, 'all')
                            spikes_idx_sel = 1:size(SpikeTrials{ipart}{ilabel}.trial{icluster},2);
                        else
                            if size(SpikeTrials{ipart}{ilabel}.trial{icluster},2) > cfg.spikewaveform.nspikes
                                spikes_idx_sel = randperm(size(SpikeTrials{ipart}{ilabel}.trial{icluster},2), cfg.spikewaveform.nspikes);
                            else
                                spikes_idx_sel = 1:size(SpikeTrials{ipart}{ilabel}.trial{icluster},2);
                            end
                        end
                        
                        if ~isempty(spikes_idx_sel)
                            
                            %define Fieldtrip trials
                            trialcount  = 0;
                            Startsample = [];
                            Endsample   = [];
                            Offset      = [];
                            Trialnr     = [];
                            
                            for ispike = spikes_idx_sel
                                trialcount   = trialcount + 1;
                                Startsample  = [Startsample; round(SpikeTrials{ipart}{ilabel}.timestamp{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(1) * hdr.Fs)];
                                Endsample    = [Endsample;   round(SpikeTrials{ipart}{ilabel}.timestamp{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(2) * hdr.Fs)];
                                Offset       = [Offset; round(cfg.spikewaveform.toi(1) * hdr.Fs)];
                                Trialnr      = [Trialnr; trialcount];
                            end
                            
                            cfgtemp                         = [];
                            
                            cfgtemp.trl                     = [Startsample, Endsample, Offset];
                            cfgtemp.trl(:,4)                = Startsample;                          % startsample
                            cfgtemp.trl(:,5)                = Endsample;                            % endsample
                            cfgtemp.trl(:,6)                = Offset;                               % offset
                            cfgtemp.trl(:,7)                = Endsample-Startsample+1;              % duration in samples
                            cfgtemp.trl(:,8)                = Trialnr;                              % trialnr. to try to find trials that are missing afterwards
                            
                            cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs
                            cfgtemp.trlunit                 = 'samples';
                            
                            %filter data
                            cfgtemp.hpfilter                = 'yes';
                            cfgtemp.hpfilttype              = 'but'; %same as Spiking Circus
                            cfgtemp.hpfiltord               = 3;     %same as Spyking Circus
                            cfgtemp.hpfreq                  = cfg.spikewaveform.cutoff;
                            cfgtemp.padding                 = 1 / cfg.spikewaveform.cutoff * 5; %5 cycles for filtering
                            
                            cfgtemp.dataset = datafile;
                            
                            SpikeWaveforms{ipart}{ilabel}{icluster}                     = ft_preprocessing(cfgtemp);
                            SpikeWaveforms{ipart}{ilabel}{icluster}.label               = [];
                            SpikeWaveforms{ipart}{ilabel}{icluster}.label{1}            = SpikeTrials{ipart}{ilabel}.label{icluster};
                            SpikeWaveforms{ipart}{ilabel}{icluster}.template_maxchan    = SpikeTrials{ipart}{ilabel}.template_maxchan(icluster);
                        else %if there are no spike in all trials of ilabel for this cluster
                            SpikeWaveforms{ipart}{ilabel}{icluster} = [];
                        end
                        
                    end %icluster
                end %~isempty(cluster_idx)
            end
        end %ichan
    end % ilabel
end %ipart

save(fname,'SpikeWaveforms','-v7.3');

end




