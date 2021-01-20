function [trialinfo,circuspath] = writemicroforspykingcircus(cfg, MuseStruct_micro)

fname_output = fullfile(cfg.datasavedir,[cfg.prefix,'single_pattern_',cfg.label,'_SC_trialinfo.mat']);
    
if exist(fname_output,'file') && cfg.force == false
    fprintf('Loading trialinfo: %s \n',fname_output);
    load(fname_output,'trialinfo','circuspath');
else
    
    % read aligned micro data at full samplerate
    % Note: no artefacts based on correlations are removed
    cfgtemp                     = cfg;
    cfgtemp.force               = cfg.forcereload;
    dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);
    
    % workaround, readMicroFs will now add fsample
    fnametemp                   = fullfile(MuseStruct_micro{1}.directory, MuseStruct_micro{1}.filenames(1).name);
    hdrtemp                     = ft_read_header(fnametemp);
    dat_microFs.fsample         = hdrtemp.Fs;
    
    % concatinate over repetitions, on preallocated memory
    l = 0;
    for itrial = 1 : size(dat_microFs.trial,2)
        l = l + size(dat_microFs.trial{itrial},2);
    end
    
    dat                         = rmfield(dat_microFs,{'trial','time','cfg'});
    dat.trial{1}                = zeros(size(dat_microFs.trial{1},1),l);
    isample = 1;
    trialinfo = zeros(size(dat_microFs.trial,2),2);
    for itrial = 1 : size(dat_microFs.trial,2)
        fprintf('Concatinating trial %d of %d\n',itrial,size(dat_microFs.trial,2))
        dat.trial{1}(:,isample:isample+size(dat_microFs.trial{itrial},2)-1) = dat_microFs.trial{itrial};
        trialinfo(itrial,1) = isample;
        trialinfo(itrial,2) = isample+size(dat_microFs.trial{itrial},2)-1;
        isample = isample + size(dat_microFs.trial{itrial},2);
    end
        
    % select data
    dat.trialdimord             = '{rpt}_chan_time';
    dat.time{1}                 = (0:size(dat.trial{1},2)-1); % *hdr_temp.Fs;
    
    cfgtemp                     = [];
    cfgtemp.channel             = cellstr(cfg.channel);
    dat                         = ft_selectdata(cfgtemp,dat);
    
    % use headerinfo of first data file
    hdr_temp                    = ft_read_header(fullfile(MuseStruct_micro{1}.directory,MuseStruct_micro{1}.filenames(1).name));
    
    % write data in .ncs format
    hdr                         = [];
    hdr.Fs                      = hdr_temp.Fs;
    hdr.nSamples                = size(dat.trial{1},2);
    hdr.nSamplePre              = 0;
    hdr.nChans                  = 1;
    hdr.FirstTimeStamp          = 0;
    hdr.TimeStampPerSample      = hdr_temp.TimeStampPerSample;
    for ichan = 1 : size(dat.trial{1},1)
        hdr.label{1} = dat_microFs.label{ichan};
        fname = fullfile(cfg.datasavedir,[cfg.prefix,'single_pattern_',cfg.label,'_',dat.label{ichan},'.ncs']);
        ft_write_data(fname,dat.trial{1}(ichan,:),'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
        circuspath{ichan} = fname;
    end
    save(fname_output,'trialinfo','circuspath');
end
