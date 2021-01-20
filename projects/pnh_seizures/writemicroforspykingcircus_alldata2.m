function [trialinfo] = writemicroforspykingcircus_allmarkers(cfg)

fname_output = fullfile(cfg.datasavedir,['SC_trialinfo_allmarkers.mat']);

if exist(fname_output,'file') && cfg.force == false
    fprintf('Loading trialinfo: %s \n',fname_output);
    load(fname_output,'trialinfo','circuspath');
else
    
    for ichan = 1 : size(cfg.channel,2)
        
        l = 0;
        for imarker = 1 : size(cfg.circuspath,2)
            fname = fullfile(cfg.circuspath{imarker}{ichan});
            fprintf('Loading %s \n',fname);
            temp{imarker} = ft_read_data(fname);
            l = l + size( temp{imarker},2);
        end
        
        hdr             = ft_read_header(fullfile(cfg.circuspath{1}{1}));
        dat             = zeros(1,l);
        
        isample = 1;
        for imarker = 1 : size(temp,2)
            fprintf('Concatinating marker %d of %d\n',imarker,size(temp,2))
            dat(:,isample:isample+size(temp{imarker},2)-1) = temp{imarker};
            trialinfo(imarker,1) = isample;
            trialinfo(imarker,2) = isample+size(temp{imarker},2)-1;
            isample = isample + size(temp{imarker},2);
        end
        
        fprintf('Writing concatinated data to: %s \n',fname);
        
        % select data
        
        % write data in .ncs format
        hdr.nSamples                = size(dat,2);
        hdr.nSamplePre              = 0;
        hdr.nChans                  = 1;
        hdr.FirstTimeStamp          = 0;
        hdr.label{1}                = cfg.channel{ichan};
        fname                       = fullfile(cfg.datasavedir,['all_concatinated_',cfg.channel{ichan},'.ncs']);
        ft_write_data(fname,dat,'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
    end
    
    save(fname_output,'trialinfo');
end
% 
% data = [];
% data.trial{1} = dat;
% data.label = hdr.label;
% data.time{1} = (1:length(dat))./hdr.Fs;
% 
% 
% ft_databrowser([],data);