function [data_avg_allchans, data_avg_chanalign, data_avg_EMG] = dtx_get_patient_data(cfg, data, ipart, imarker)
%comparison EEG EMG for now only with motor cortex electrodes (C3 or C4)

data = data{ipart}{imarker};

%avg all channels
%enveloppe EMG
%chan align avg

%% avg all channels
data_avg_allchans               = ft_timelockanalysis([],data);
data_avg_allchans.ID            = cfg.prefix(1:end-1);


%% avg chan align
%select channel
cfgtemp                     = [];
cfgtemp.channel             = cfg.align.channel{imarker};
data_chanalign              = ft_selectdata(cfgtemp,data);

%avg
cfgtemp                     = [];
data_avg_chanalign          = ft_timelockanalysis(cfgtemp,data_chanalign);

data_avg_chanalign.label = {'chan_align'};
data_avg_chanalign.ID = cfg.prefix(1:end-1);


%% enveloppe EMG

if any(strfind(cfg.LFP.name{imarker}, 'EMG'))
    
    %select channel
    cfgtemp                     = [];
    cfgtemp.channel             = cfg.LFP.emg{imarker};
    data_EMG                    = ft_selectdata(cfgtemp,data);
    
    %compute envelope for each trial
    for itrial = 1 : size(data_EMG.trial,2)
        rect_emg = abs(data_EMG.trial{itrial}(1,:));
        [env{itrial}, ~] = envelope(rect_emg,cfg.EMG.envparam,cfg.EMG.envmethod);
    end

    %compute avg of envelope
    for ioffset = 1:length(env{1}) %all trials must have the same size
        for itrial = 1:length(env)
            env_by_offset(itrial) = env{itrial}(ioffset);
        end
        env_avg(ioffset) = mean(env_by_offset);
    end   
    
    %return good data structure
    %remove EMG channel and replace with env
    cfgtemp                     = [];
    cfgtemp.channel             = {'eeg1020'};
    data_avg_EMG                = ft_selectdata(cfgtemp,data_avg_allchans);

    data_avg_EMG.label{end+1} = 'EMG';   
    data_avg_EMG.avg(end+1,:) = env_avg; %first channel is EEG channel  
    %to be consistent, but those infos are not used :
    data_avg_EMG.var(end+1,:) = env_avg;
    data_avg_EMG.dof(end+1,:) = env_avg;
    
    data_avg_EMG.ID = cfg.prefix(1:end-1);
    
    data_avg_allchans = [];
    data_avg_chanalign = [];
    
else
    data_avg_EMG = [];
    
end



end