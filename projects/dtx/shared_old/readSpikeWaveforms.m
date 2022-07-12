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
% cfg.spikewaveform.part_list   = list of parts to analyse. Can be an array
%                                 of integers, or 'all'. Default = 'all'.
%
% ### OUTPUT
% SpikeWaveforms{parts}{labels} = Fieldtrip raw structure. For each labels
%                                 (corresponding to cfg.name), there is one
%                                 cell per unit, one trial per spike.
%

%get defaults cfg parameters
cfg.spikewaveform               = ft_getopt(cfg, 'spikewaveform', []);
cfg.spikewaveform.toi           = ft_getopt(cfg.spikewaveform, 'toi'    	, [-0.0015 0.0015]);
cfg.spikewaveform.nspikes       = ft_getopt(cfg.spikewaveform, 'nspikes' 	, 1000);
cfg.spikewaveform.hpfilttype    = ft_getopt(cfg.spikewaveform, 'hpfilttype'	, 'but');
cfg.spikewaveform.hpfiltord     = ft_getopt(cfg.spikewaveform, 'hpfiltord'  , 3);
cfg.spikewaveform.hpfreq        = ft_getopt(cfg.spikewaveform, 'hpfreq'     , 300);
cfg.spikewaveform.part_list     = ft_getopt(cfg.spikewaveform , 'part_list' , 'all');
cfg.circus.channelname          = ft_getopt(cfg.circus, 'channelname', []);

if strcmp(cfg.spikewaveform.part_list, 'all')
    cfg.spikewaveform.part_list = 1:size(SpikeRaw, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'spike_waveform.mat']);

if exist(fname, 'file') && force == false
    fprintf('Loading precomputed spike waveforms\n');
    load(fname, 'SpikeWaveforms');
    return

elseif exist(fname, 'file') && force == true
    fprintf('Forced recomputing of spike waveforms\n');

else
    fprintf('Computing spike waveforms\n');
end

for ipart = cfg.spikewaveform.part_list

    for markername = string(fieldnames(SpikeRaw{ipart}))'

        if isempty(SpikeRaw{ipart}.(char(markername)))
            continue
        end

        for icluster = 1 : size(SpikeRaw{ipart}.(char(markername)).label, 2)

            % only load new data when needed
            loadnew = false;
            if icluster > 1
                % if channel is different from loaded data
                if ~isempty(cfg.circus.channelname)
                    if ~strcmp(SpikeRaw{ipart}.(markername).channelname{icluster}, SpikeRaw{ipart}.(markername).channelname{icluster-1}); loadnew = true; end
                end
                % if max electrode plot is different from loaded data
                if SpikeRaw{ipart}.(char(markername)).template_maxchan(icluster) ~= SpikeRaw{ipart}.(char(markername)).template_maxchan(icluster-1); loadnew = true; end
            else
                loadnew = true;
            end

            if loadnew
                %FIXME : toremove:
                [~,datadir] = fileparts(cfg.datasavedir);
                if ~strcmp(datadir, 'data')
                    datasavedir = fileparts(cfg.datasavedir);
                else
                    datasavedir = cfg.datasavedir;
                end
                
                if isempty(cfg.circus.channelname)
                    chanindx_cfg    = SpikeRaw{ipart}.(char(markername)).template_maxchan(icluster) + 1;
                    temp            = dir(fullfile(datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], [cfg.prefix, 'p', num2str(ipart), '-multifile-', cfg.circus.channel{chanindx_cfg},'.ncs']));
                else
                    % make sure to pick the right file
                    channelname     = SpikeRaw{ipart}.(char(markername)).channelname{icluster};
                    maxchan         = SpikeRaw{ipart}.(char(markername)).template_maxchan(icluster);
                    chanoffset      = find(strcmp(cfg.circus.channelname, channelname), 1, 'first');
                    chanindx_cfg    = maxchan + chanoffset;
                    temp            = dir(fullfile(datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], channelname, [cfg.prefix, 'p', num2str(ipart), '-multifile-', cfg.circus.channel{chanindx_cfg} ,'.ncs']));          
                end
                datafile        = fullfile(temp.folder, temp.name);

                if isfield(SpikeRaw{ipart}.(char(markername)),'hdr')
                    hdr         = SpikeRaw{ipart}.(char(markername)).hdr;
                else
                    hdr         = ft_read_header(datafile);
                end

                clear chandata
                fprintf('Loading spike waveforms from channel %s\n', datafile);
                cfgtemp                         = [];
                cfgtemp.dataset                 = datafile;
                cfgtemp.hpfilter                = 'yes';
                cfgtemp.hpfilttype              = cfg.spikewaveform.hpfilttype;
                cfgtemp.hpfiltord               = cfg.spikewaveform.hpfiltord;
                cfgtemp.hpfreq                  = cfg.spikewaveform.hpfreq;
                chandata                        = ft_preprocessing(cfgtemp);
            end

            % Select random spike if required
            if strcmp(cfg.spikewaveform.nspikes, 'all')
                spikes_idx_sel = 1:length(SpikeRaw{ipart}.(markername).sample{icluster});
            else
                if size(SpikeRaw{ipart}.(markername).trial{icluster}, 2) > cfg.spikewaveform.nspikes
                    spikes_idx_sel = randperm(length(SpikeRaw{ipart}.(markername).sample{icluster}), cfg.spikewaveform.nspikes);
                else
                    spikes_idx_sel = 1:length(SpikeRaw{ipart}.(markername).sample{icluster});
                end
            end

            if isempty(spikes_idx_sel)
                SpikeWaveforms{ipart}.(markername){icluster} = [];
                continue
            end

            % equally spaced - make as option with nr. of waveforms?
            % spikes_idx_sel = 1 : round(size(SpikeRaw{ipart}.(markername).trial{icluster}, 2)/100) : size(SpikeRaw{ipart}.(markername).trial{icluster}, 2);

            % define Fieldtrip trials
            trialcount  = 0;
            Startsample = [];
            Endsample   = [];
            Offset      = [];
            Trialnr     = [];

            % note that the sample field has to been added by ft_spike_maketrials, so use the adapted version in /external/fieldtrip:
            % note that the indexing of the channels can be different.
            for ispike = spikes_idx_sel
                trialcount   = trialcount + 1;
                %                     Startsample  = [Startsample; round(SpikeRaw{ipart}.(char(markername)).sample{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(1) * hdr.Fs)];
                %                     Endsample    = [Endsample;   round(SpikeRaw{ipart}.(char(markername)).sample{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(2) * hdr.Fs)];
                Startsample  = [Startsample; round(SpikeRaw{ipart}.(char(markername)).sample{icluster}(ispike) + cfg.spikewaveform.toi(1) * hdr.Fs)];
                Endsample    = [Endsample;   round(SpikeRaw{ipart}.(char(markername)).sample{icluster}(ispike) + cfg.spikewaveform.toi(2) * hdr.Fs)];
                Offset       = [Offset; round(cfg.spikewaveform.toi(1) * hdr.Fs)];
                Trialnr      = [Trialnr; trialcount];
            end

            full_trial = Startsample > 0 & Endsample < length(chandata.trial{1});% don't read before BOF or after EOF
            if sum(full_trial) == 0
                SpikeWaveforms{ipart}.(char(markername)){icluster} = [];
                continue
            end

            cfgtemp                                                                     = [];
            cfgtemp.trl                                                                 = [Startsample, Endsample, Offset];
            cfgtemp.trl                                                                 = cfgtemp.trl(full_trial, :); % so not to read before BOF or after EOFs
            cfgtemp.trlunit                                                             = 'samples';
            SpikeWaveforms{ipart}.(char(markername)){icluster}                          = ft_redefinetrial(cfgtemp, chandata);
            SpikeWaveforms{ipart}.(char(markername)){icluster}.label                    = [];
            SpikeWaveforms{ipart}.(char(markername)){icluster}.label{1}                 = SpikeRaw{ipart}.(char(markername)).label{icluster};
            SpikeWaveforms{ipart}.(char(markername)){icluster}.template_maxchan         = SpikeRaw{ipart}.(char(markername)).template_maxchan(icluster);
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo                = table;
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo.begsample      = Startsample(full_trial);
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo.endsample      = Endsample(full_trial);
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo.offset         = Offset(full_trial);
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo.trialduration  = Endsample(full_trial)-Startsample(full_trial)+1;
            SpikeWaveforms{ipart}.(char(markername)){icluster}.trialinfo.trialnr        = Trialnr(full_trial);

            %remove cfg field which takes looooot of space on disk
            %(several Go per patient)
            SpikeWaveforms{ipart}.(char(markername)){icluster} = rmfield(SpikeWaveforms{ipart}.(char(markername)){icluster}, 'cfg');

        end %icluster
    end %markername
end %ipart

save(fname, 'SpikeWaveforms', '-v7.3');
