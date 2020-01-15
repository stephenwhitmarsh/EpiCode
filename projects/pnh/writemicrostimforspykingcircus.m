function [trialinfo,dat] = writemicrostimforspykingcircus(cfg, markerdat, micromed_markers)

%% WILL ONLY RETURN TRIALINFO OF 50Hz

for stimrate = [1, 50]
    
    fname_output = fullfile(cfg.datasavedir,[stimrate,'_SC_trialinfo.mat']);
    
    % fname_output = fullfile(cfg.datasavedir,'SC_stim_trialinfo.mat');
    if exist(fname_output,'file') && cfg.force == false
        fprintf('Loading trialinfo: %s',fname_output);
        load(fname_output,'trialinfo');
    else
        
        % select trials
        triallist = find(micromed_markers.Frequency == stimrate)';
        
        % concatinate over repetitions, on preallocated memory
        l = 0;
        for itrial = triallist
            l = l + size(markerdat{itrial}.trial{1},2);
        end
        
        dat = rmfield(markerdat{itrial},{'trial','time','cfg','sampleinfo'});
        dat.trial{1} = zeros(size(markerdat{itrial}.trial{1},1),l);
        
        isample = 1;
        icounter = 1;
        for itrial = triallist
            fprintf('Concatinating trial %d of %d at %d Hz \n',icounter,length(triallist),stimrate)
            dat.trial{1}(:,isample:isample+size(markerdat{itrial}.trial{1},2)-1) = markerdat{itrial}.trial{1};
            trialinfo(icounter,1)   = isample;
            trialinfo(icounter,2)   = isample + size(markerdat{itrial}.trial{1},2)-1;
            isample                 = isample + size(markerdat{itrial}.trial{1},2);
            icounter                = icounter + 1;
        end
        
        % use headerinfo of first data file
        neuralynx_dirlist_temp      = dir2(fullfile(cfg.patientdir,'eeg'));
        neuralynx_datafiles         = dir2(fullfile(neuralynx_dirlist_temp(1).folder,neuralynx_dirlist_temp(1).name,'*_m*.ncs'));
        hdr_temp                    = ft_read_header(fullfile(neuralynx_dirlist_temp(1).folder,neuralynx_dirlist_temp(1).name,neuralynx_datafiles(1).name)); % read header from first ncs file
        
        % select data
        dat.trialdimord             = '{rpt}_chan_time';
        dat.time{1}                 = (0:size(dat.trial{1},2)-1)/hdr_temp.Fs;
        
        cfgtemp                     = [];
        cfgtemp.channel             = cfg.channel;
        dat                         = ft_selectdata(cfgtemp,dat);
        
        % write deadtime file
        fname           = fullfile(cfg.datasavedir,[num2str(stimrate),'Hz_concatinated_deadtime.txt']);
        deadtimes       = [];
        deadtimes(:,1)  = trialinfo(:,1) + cfg.start_deadtime * hdr_temp.Fs;
        deadtimes(:,2)  = trialinfo(:,1) + cfg.end_deadtime * hdr_temp.Fs;
        dlmwrite(fname,deadtimes,'delimiter','	','precision','%5.f');

        % write data in .ncs format
        hdr                         = [];
        hdr.Fs                      = hdr_temp.Fs;
        hdr.nSamples                = size(dat.trial{1},2);
        hdr.nSamplePre              = 0;
        hdr.nChans                  = 1;
        hdr.FirstTimeStamp          = 0;
        hdr.TimeStampPerSample      = hdr_temp.TimeStampPerSample;
        for ichan = 1 : size(dat.trial{1},1)
            hdr.label{1} = dat.label{ichan};
            fname = fullfile(cfg.datasavedir,[num2str(stimrate),'Hz_concatinated_',dat.label{ichan},'.ncs']);
            ft_write_data(fname,dat.trial{1}(ichan,:),'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
        end

    end
    save(fname_output,'trialinfo');
end

fprintf('Done.\n');