function [SpikeRaw, SpikeTrials] = readSpykingCircus(cfg,MuseStruct,force)

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
    
    % find spiking-circus output path, which is based on the name of the
    % first datafile

    temp = dir(fullfile(cfg.circus.outputdir,[cfg.prefix,'multifile-',cfg.circus.channel{1},'.result',cfg.circus.suffix,'.hdf5']));

    if isempty(temp)
        fprintf('Could not find Spyking-Circus results: %s\n',fullfile(cfg.datasavedir,cfg.circus.outputdir,[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.result.hdf5']));
        return
    else
        fname_spikes = fullfile(temp.folder,temp.name);
    end
%     temp = dir(fullfile(cfg.circus.outputdir,[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.templates',cfg.circus.suffix,'.hdf5']));

    temp = dir(fullfile(cfg.circus.outputdir,[cfg.prefix,'multifile-',cfg.circus.channel{1},'.templates.hdf5']));
    if isempty(temp)
        fprintf('Could not find Spyking-Circus templates: %s\n',fullfile(cfg.datasavedir,'SpykingCircus',[cfg.prefix,'all_data_',cfg.circus.channel{1}(1:end-2),'_*.templates.hdf5']));
        return
    end
    fname_templates = fullfile(temp.folder,temp.name);
    
    % load spiking data
    if exist(fname_spikes,'file')
        fprintf('Loading spike data from: %s\n',fname_spikes);
        datinfo     = h5info(fname_spikes);
        temp        = dir(fullfile(cfg.datasavedir,[cfg.prefix,'multifile-',cfg.circus.channel{1},'.ncs']));
        hdr_fname   = fullfile(temp(1).folder,temp(1).name);
        hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        timestamps  = ft_read_data(hdr_fname,'timestamp','true');  % take the first concatinated file to extract the timestamps
%         timestamps  =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster
        
        % find hdf5 index for spike times
        for i = 1 : size(datinfo.Groups,1)
            names(i) = string(datinfo.Groups(i).Name);
            if strfind(names(i),'spiketimes')
                spiketimes_indx = i;
            end
        end
        
        % find index of clusters 
        for i = 1 : size(datinfo.Groups(spiketimes_indx).Datasets,1) % number of templates
            SpikeRaw.label{i} = datinfo.Groups(spiketimes_indx).Datasets(i).Name;
            temp = strsplit(SpikeRaw.label{i},'_');
            clusternr(i) = str2double(temp{2});
        end
        
        for i = 1:numel(clusternr)
            % read spike timings (in seconds)
            datasetname = char(strcat('/spiketimes/',SpikeRaw.label{i}));
            SpikeRaw.sample{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % read amplitudes
            datasetname = char(strcat('/amplitudes/',SpikeRaw.label{i}));
            SpikeRaw.amplitude{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % map samplenrs onto timestamps
            SpikeRaw.timestamp{clusternr(i)+1} = timestamps(SpikeRaw.sample{clusternr(i)+1});
           % SpikeRaw.timestamp{clusternr(i)+1} = int64(SpikeRaw.sample{clusternr(i)+1}) * int64(hdr.TimeStampPerSample) + int64(hdr.FirstTimeStamp);         
        end
        
        % load templates
        templates_size  = double(h5read(fname_templates, '/temp_shape'));
        N_e             = templates_size(2);
        N_t             = templates_size(1);
        temp_x          = double(h5read(fname_templates, '/temp_x') + 1);
        temp_y          = double(h5read(fname_templates, '/temp_y') + 1);
        temp_z          = double(h5read(fname_templates, '/temp_data'));
        templates       = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
        templates_size  = [templates_size(1) templates_size(2) templates_size(3)/2];
        
        for itemp = 1:numel(clusternr)
            template = full(reshape(templates(:, itemp), templates_size(2), templates_size(1)))';
            [~,i] = max(mean(abs(template),2));
            SpikeRaw.template(itemp,:,:) = template;
            SpikeRaw.template_maxchan(itemp) = i;
        end
               
        % capitalize all markernames
        for idir = 1 : size(MuseStruct,2)
            markernames_old = fields(MuseStruct{idir}.markers);
            markernames_new = upper(markernames_old);
            for i = 1 : size(markernames_old,1)
                MuseStruct{idir}.markers.(markernames_new{i}) = MuseStruct{idir}.markers.(markernames_old{i});
                MuseStruct{idir}.markers = rmfield(MuseStruct{idir}.markers,markernames_old{i});
            end
        end
        for ilabel = 1 : size(cfg.name,2) 
            cfg.muse.startend{ilabel,1} = upper(cfg.muse.startend{ilabel,1});
            cfg.muse.startend{ilabel,2} = upper(cfg.muse.startend{ilabel,2});
        end
        
        clear Trials
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
            Filenr = [];
            FileOffset = [];
            
            dirOnset = 0;
            trialcount = 1;
            for idir = 1 : size(MuseStruct,2)
                try
                    for ievent = 1 : size(MuseStruct{idir}.markers.(cfg.muse.startend{ilabel}).events,2)
                        Startsample  = [Startsample; MuseStruct{idir}.markers.(cfg.muse.startend{ilabel,1}).offset(ievent) + cfg.epoch.toi{ilabel}(1) * hdr.Fs + dirOnset];
                        Endsample    = [Endsample;   MuseStruct{idir}.markers.(cfg.muse.startend{ilabel,2}).offset(ievent) + cfg.epoch.toi{ilabel}(2) * hdr.Fs + dirOnset];
                        Offset       = [Offset; cfg.epoch.toi{ilabel}(1) * hdr.Fs];
                        Trialnr      = [Trialnr; trialcount];
                        Filenr       = [Filenr; idir];
                        trialcount   = trialcount + 1;
                        FileOffset   = [FileOffset; dirOnset];
                    end
                catch
                    fprintf('No events starting with %s found in filenr %d\n',cfg.muse.startend{ilabel},idir);
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
            cfgtemp.trl(:,11)               = Filenr;                               % trialnr. to try to find trials that are missing afterwards
            cfgtemp.trl(:,12)               = FileOffset;                           % trialnr. to try to find trials that are missing afterwards

            maxsamples                      = sum(cfg.sampleinfo(:,2));
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


