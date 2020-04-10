function [data_avg_allchans, data_avg_chantoplot, data_TFR, data_avg_EMG] = dtx_get_patient_data(cfg, data, ipart, imarker)
%comparison EEG EMG for now only with motor cortex electrodes (C3 or C4)

data = data{ipart}{imarker};

%avg all channels
%enveloppe EMG
%chan align avg

%% avg all channels

%keep only channels common for all the rat/patients
if isfield(cfg, 'commonchans')
    cfgtemp                         = [];
    cfgtemp.channel                 = cfg.commonchans;
    data_commonchans                = ft_selectdata(cfgtemp,data);
else
    data_commonchans                = data;
end

%do avg
data_avg_allchans               = ft_timelockanalysis([],data_commonchans);
data_avg_allchans.ID            = cfg.prefix(1:end-1);



%% avg chan align
%select channel
cfgtemp                      = [];
cfgtemp.channel              = cfg.LFP.electrodetoplot{imarker};
data_chantoplot              = ft_selectdata(cfgtemp,data);

%avg
cfgtemp                      = [];
data_avg_chantoplot          = ft_timelockanalysis(cfgtemp,data_chantoplot);

data_avg_chantoplot.label = cfg.LFP.name(imarker);
data_avg_chantoplot.ID = cfg.prefix(1:end-1);

%% TFR
if isfield(cfg.LFP, 'TFR')
    if cfg.LFP.TFR.doTFR == true
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all';
        cfgtemp.method                  = 'mtmconvol'; %relative change. 'mtmfft' : frequency content
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'hanning';
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.keeptrials              = 'no';
        cfgtemp.foi                     = 0:0.2:50;%80:0.2:125;
        cfgtemp.t_ftimwin               = 40./cfgtemp.foi;
        %cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*0.5;
        cfgtemp.toi                     = 'all';
        data_TFR                        = ft_freqanalysis(cfgtemp,data_chantoplot);
        
        data_TFR.label = cfg.LFP.name(imarker);
    else
        data_TFR = [];
    end
else
    data_TFR = [];
end

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
    cfgtemp.channel             = {'all','-*EMG*'};% {'eeg1020'};
    data_avg_EMG                = ft_selectdata(cfgtemp,data_avg_allchans);
 
    data_avg_EMG.label{end+1} = 'EMG';   
    data_avg_EMG.avg(end+1,:) = env_avg; %first channel is EEG channel  
    %to be consistent, but those infos are not used :
    data_avg_EMG.var(end+1,:) = env_avg;
    data_avg_EMG.dof(end+1,:) = env_avg;
    
    data_avg_EMG.ID = cfg.prefix(1:end-1);
%     
%     data_avg_allchans = [];
%     data_avg_chantoplot = [];
    
else
    data_avg_EMG = [];
    
end



end