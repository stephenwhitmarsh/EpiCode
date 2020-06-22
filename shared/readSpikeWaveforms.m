function [SpikeWaveforms] = readSpikeWaveforms(cfg,SpikeTrials,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [SpikeWaveforms] = readSpikeWaveforms(cfg,SpikeTrials,force,parts_to_read)
% Cut continuous data into trials for each spike sample.
%
% ### Necessary input
% cfg.name                      = label(s) of the analysis
% cfg.prefix                    = prefix to output files
% cfg.datasavedir               = data directory to the ncs data used by
%                                 Spyking-Circus, and to save output data
% cfg.circus.channel            = channel names which were analyzed by
%                                 Spyking-Circus
% force                         = whether to redo analyses or read previous
%                                 save (true/false)
%
% ### Optional cfg fields
% cfg.spikewaveform.toi         = in seconds, time to load before and after
%                                 the peak of the spike. Default =
%                                 [-0.0015 0.0015]
% cfg.spikewaveform.nspikes     = maximum number of spike waveforms to load.
%                                 Can be 'all'. If there are more spikes
%                                 than nspikes, a random selection of
%                                 spikes is done. Default = 1000.
% cfg.spikewaveform.hpfreq      = high pass filter frequency to apply to
%                                 raw data. Default = 300.
% cfg.spikewaveform.hpfilttype  = Default = 'but', same as Spiking Circus
% cfg.spikewaveform.hpfiltord   = Default = 3, same as Spyking Circus
% cfg.spikewaveform.part_list   = list of parts to analyse. Can be an array
%                                 of integers, or 'all'. Default = 'all'.
%
% ### OUTPUT
% SpikeWaveforms{parts}{labels} = Fieldtrip raw structure. For each labels
%                                 (corresponding to cfg.name), there is one
%                                 cell per unit, one trial per spike.
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get defaults cfg parameters
cfg.spikewaveform               = ft_getopt(cfg, 'spikewaveform', []);
cfg.spikewaveform.toi           = ft_getopt(cfg.spikewaveform, 'toi'    	, [-0.0015 0.0015]);
cfg.spikewaveform.nspikes       = ft_getopt(cfg.spikewaveform, 'nspikes' 	, 1000);
cfg.spikewaveform.hpfilttype    = ft_getopt(cfg.spikewaveform, 'hpfilttype'	, 'but');
cfg.spikewaveform.hpfiltord     = ft_getopt(cfg.spikewaveform, 'hpfiltord'  , 3);
cfg.spikewaveform.hpfreq        = ft_getopt(cfg.spikewaveform, 'hpfreq'     , 300);
cfg.spikewaveform.part_list    = ft_getopt(cfg.spikewaveform , 'part_list' , 'all');


if strcmp(cfg.spikewaveform.part_list,'all')
    cfg.spikewaveform.part_list = 1:size(SpikeTrials,2);
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spike_waveform.mat']);

if exist(fname,'file') && force == false
    fprintf('*******************************************\n');
    fprintf('*** Loading precomputed spike waveforms ***\n');
    fprintf('*******************************************\n\n');
    
    load(fname,'SpikeWaveforms');
    return
    
elseif exist(fname,'file') && force == true
    fprintf('*******************************************\n');
    fprintf('*** Forced  recomputing spike waveforms ***\n');
    fprintf('*******************************************\n\n');
    
else
    fprintf('*********************************\n');
    fprintf('*** Computing spike waveforms ***\n');
    fprintf('*********************************\n\n');
end


for ipart = cfg.spikewaveform.part_list
    
    for ilabel = 1:size(cfg.name,2)
        
        for ichan = 1:size(cfg.circus.channel,2)
            
            % find concatenated channel used by Spyking Circus
            temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{ichan},'.ncs']));
            if isempty(temp)
                continue
            end
            
            datafile = fullfile(temp.folder,temp.name);
            fprintf('Found datafile: %s\n',datafile);
            
            hdr = ft_read_header(datafile);
            
            clusters_idx = [];
            clusters_idx = find(SpikeTrials{ipart}{ilabel}.template_maxchan == ichan - 1); %-1 because Phy starts counting at zero
            
            if isempty(clusters_idx)
                continue
            end
            
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
                
                if isempty(spikes_idx_sel)
                    SpikeWaveforms{ipart}{ilabel}{icluster} = [];
                    continue
                end
                
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
                cfgtemp.hpfilttype              = cfg.spikewaveform.hpfilttype;
                cfgtemp.hpfiltord               = cfg.spikewaveform.hpfiltord;
                cfgtemp.hpfreq                  = cfg.spikewaveform.hpfreq;
                cfgtemp.padding                 = 1 / cfgtemp.hpfreq * 5; %5 cycles for filtering
                
                cfgtemp.dataset = datafile;
                
                SpikeWaveforms{ipart}{ilabel}{icluster}                     = ft_preprocessing(cfgtemp);
                SpikeWaveforms{ipart}{ilabel}{icluster}.label               = [];
                SpikeWaveforms{ipart}{ilabel}{icluster}.label{1}            = SpikeTrials{ipart}{ilabel}.label{icluster};
                SpikeWaveforms{ipart}{ilabel}{icluster}.template_maxchan    = SpikeTrials{ipart}{ilabel}.template_maxchan(icluster);
                
            end %icluster
        end %ichan
    end % ilabel
end %ipart

save(fname,'SpikeWaveforms','-v7.3');

end




