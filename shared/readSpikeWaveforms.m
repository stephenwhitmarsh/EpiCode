function [SpikeWaveforms] = readSpikeWaveforms(cfg, SpikeRaw, force)

% [SpikeWaveforms] = readSpikeWaveforms(cfg, SpikeTrials, force, parts_to_read)
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
% cfg.spikewaveform.lpfreq      = low pass filter frequency to apply to
%                                 raw data. Default = 3000.
% cfg.spikewaveform.part_list   = list of parts to analyse. Can be an array
%                                 of integers, or 'all'. Default = 'all'.
%
% ### OUTPUT
% SpikeWaveforms{parts}{labels} = Fieldtrip raw structure. For each labels
%                                 (corresponding to cfg.name), there is one
%                                 cell per unit, one trial per spike.
%

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'spike_waveform.mat']);

% 
% if exist(fname, 'file') && force == false
%     fprintf('Loading precomputed spike waveforms\n');
%     load(fname, 'SpikeWaveforms');
%     return
%     
% elseif exist(fname, 'file') && force == true
%     fprintf('Forced recomputing of spike waveforms\n');
%     
% else
%     fprintf('Computing spike waveforms\n');
% end
% 

if nargin == 1
    if exist(fname, 'file')
        fprintf('Reading %s\n', fname);
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname, 'SpikeWaveforms');
            catch ME
                err_count = err_count + 1;
                disp('Something went wrong loading the file. Trying again...')    
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        stats = {};
        return
    end
end

if exist(fname, 'file') && force == false
    fprintf('Loading %s\n', fname);
    count = 0;
    err_count = 0;
    while count == err_count
        try
            load(fname, 'SpikeWaveforms');
        catch ME
            err_count = err_count + 1;
            disp('Something went wrong loading the file. Trying again...')            
        end
        count = count + 1;
    end
    return
end

% get references from files to samplenumbers (need to re-read in case of
% change of OS due to different paths
cfg.circus.timestamps           = false;
[filelist, sampleinfo, ~, ~]    = writeSpykingCircusFileList(cfg, true);

%get defaults cfg parameters
cfg.spikewaveform               = ft_getopt(cfg, 'spikewaveform', []);
cfg.spikewaveform.toi           = ft_getopt(cfg.spikewaveform, 'toi'    	, [-0.0015 0.0015]);
cfg.spikewaveform.nspikes       = ft_getopt(cfg.spikewaveform, 'nspikes' 	, 1000);
cfg.spikewaveform.hpfilttype    = ft_getopt(cfg.spikewaveform, 'hpfilttype'	, 'but');
cfg.spikewaveform.hpfiltord     = ft_getopt(cfg.spikewaveform, 'hpfiltord'  , 3);
cfg.spikewaveform.hpfreq        = ft_getopt(cfg.spikewaveform, 'hpfreq'     , 300);
cfg.spikewaveform.lpfreq        = ft_getopt(cfg.spikewaveform, 'lpfreq'     , 6000);
cfg.spikewaveform.lpfilttype    = ft_getopt(cfg.spikewaveform, 'lpfilttype'	, 'but');
cfg.spikewaveform.lpfiltord     = ft_getopt(cfg.spikewaveform, 'lpfiltord'  , 3);
cfg.spikewaveform.part_list     = ft_getopt(cfg.spikewaveform , 'part_list' , 'all');
cfg.circus.channelname          = ft_getopt(cfg.circus, 'channelname', []);
cfg.circus.correct_chunk        = ft_getopt(cfg.circus, 'correct_chunk', false);

if strcmp(cfg.spikewaveform.part_list, 'all')
    cfg.spikewaveform.part_list = 1:size(SpikeRaw, 2);
end

for ipart = cfg.spikewaveform.part_list
    
    if isempty(SpikeRaw{ipart})
        continue
    end
    
    try
        if ~isempty(fieldnames(filelist{ipart}))
            temp = [];
            for fn = string(fieldnames(filelist{ipart}))'
                temp = [temp, filelist{ipart}.(fn)];
            end
            filelist{ipart} = temp;
            clear temp
        end
    catch
    end
    
    for iunit = 1 : size(SpikeRaw{ipart}.label, 2)
        
        chanindx        = SpikeRaw{ipart}.template_maxchan(iunit) + 1;
        
        if cfg.circus.correct_chunk{ipart} %% FIXME: WORKAROUND FOR BUG IN SPYKING CIRCUS
            cumsumsample = cumsum(sampleinfo{ipart}(:, 2)-512);
        else
            cumsumsample = cumsum(sampleinfo{ipart}(:, 2));
        end
        diroffset       = [0; cumsumsample(1:end-1)];
        cumsumfile      = filelist{ipart}(:, chanindx);
        spikes_idx_sel  = sort(randperm(min(size(SpikeRaw{ipart}.sample{iunit}, 2), cfg.spikewaveform.nspikes)), cfg.spikewaveform.nspikes);
        fileidx         = nan(1, length(spikes_idx_sel));
  
        for idx = 1 : length(spikes_idx_sel)
            fileidx(idx) = find(SpikeRaw{ipart}.sample{iunit}(spikes_idx_sel(idx)) < cumsumsample, 1, 'first');
        end
        
        filecnt = 1;
        for ifile = unique(fileidx)

            datafile = char(cumsumfile(ifile));
            fprintf('Loading spike waveforms from channel %s\n', datafile);
            hdr = ft_read_header(datafile);
            
            cfgtemp                         = [];
            cfgtemp.dataset                 = datafile;
            
            cfgtemp.hpfilter                = 'yes';
            cfgtemp.hpfilttype              = cfg.spikewaveform.hpfilttype;
            cfgtemp.hpfiltord               = cfg.spikewaveform.hpfiltord;
            cfgtemp.hpfreq                  = cfg.spikewaveform.hpfreq;
            
            cfgtemp.lpfilter                = 'yes';
            cfgtemp.lpfreq                  = cfg.spikewaveform.lpfreq;
            cfgtemp.lpfilttype              = cfg.spikewaveform.lpfilttype;
            cfgtemp.lpfiltord               = cfg.spikewaveform.lpfiltord;
            
            chandata                        = ft_preprocessing(cfgtemp);
            
            % define Fieldtrip trials
            Startsample = [];
            Endsample   = [];
            Offset      = [];
            Trialnr     = [];
            
            % note that the sample field has to been added by ft_spike_maketrials, so use the adapted version in /external/fieldtrip:
            % note that the indexing of the channels can be different.
            for ispike = spikes_idx_sel(fileidx == ifile)
                Startsample  = [Startsample; round(SpikeRaw{ipart}.sample{iunit}(ispike) + cfg.spikewaveform.toi(1) * hdr.Fs) - diroffset(ifile)];
                Endsample    = [Endsample;   round(SpikeRaw{ipart}.sample{iunit}(ispike) + cfg.spikewaveform.toi(2) * hdr.Fs) - diroffset(ifile)];
                Offset       = [Offset;      cfg.spikewaveform.toi(1) * hdr.Fs];
            end
            full_trial = Startsample > 0 & Endsample < length(chandata.trial{1}); % don't read before BOF or after EOF
            
            if sum(full_trial) == 0
                temp{ipart}.(char(markername)){iunit} = [];
                continue
            end
            
            chandata.hdr.TimeStampPerSample = 1;
            chandata.hdr.FirstTimeStamp     = 0;
            
            cfgtemp                                 = [];
            cfgtemp.trl                             = [int64(Startsample), int64(Endsample), int64(Offset)];
            cfgtemp.trl                             = cfgtemp.trl(full_trial, :); % so not to read before BOF or after EOFs
            temp{filecnt}                           = ft_redefinetrial(cfgtemp, chandata);
            
            temp{filecnt}.label                     = [];
            temp{filecnt}.label{1}                  = SpikeRaw{ipart}.label{iunit};
            temp{filecnt}.template_maxchan          = SpikeRaw{ipart}.template_maxchan(iunit);
            
            temp{filecnt}.trialinfo                 = table;
            temp{filecnt}.trialinfo.begsample       = Startsample(full_trial);
            temp{filecnt}.trialinfo.endsample       = Endsample(full_trial);
            temp{filecnt}.trialinfo.offset          = Offset(full_trial);
            temp{filecnt}.trialinfo.trialduration   = double(Endsample(full_trial) - Startsample(full_trial)+1) * 1 / hdr.Fs;
            temp{filecnt}.trialinfo.datafile        = repmat(datafile, height(temp{filecnt}.trialinfo), 1);            
            temp{filecnt}                           = rmfield(temp{filecnt}, 'cfg');
            
            filecnt = filecnt + 1;
        end % ifile

        SpikeWaveforms{ipart}{iunit} = ft_appenddata([], temp{:});
        SpikeWaveforms{ipart}{iunit}.cluster_group = SpikeRaw{ipart}.cluster_group{iunit};
        
%         if iunit == 2
%             figure; plot(SpikeWaveforms{ipart}{iunit}.time{1}, mean(vertcat(SpikeWaveforms{ipart}{iunit}.trial{:})));
%         end
        clear temp
    end %iunit
end %ipart

save(fname, 'SpikeWaveforms', '-v7.3');
