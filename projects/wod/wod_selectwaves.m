if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparams;

for irat= 4: size(config,2)
    
    %% loading LFP data and Muse markers
    
    fprintf('Load LFP data for rat %d/%d\n', irat,size(config,2));
    temp=load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'LFP_WoD','.mat']));
    LFP = temp.LFP{1}.WoD;
    clear temp
    
    
    fprintf('Load Muse Markers for rat %d/%d\n', irat, size(config,2));
    temp= load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'MuseStruct','.mat']));
    MuseStruct =temp.MuseStruct{1};
    clear temp
    
  %rename chans according to their real deepness.
        %the name is in cfg.LFP.channel, and it is renamed with the name at
        %the same index in cfg.LFP.rename
        %16 is surface, 1 is the deepest. 0 is the respi.
        n_chans = size(config{irat}.LFP.allchannel,2);
        for ichan = 1:n_chans
            if any(strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan}))
                %search channel into config
                chan_idx = strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan});
                new_name = config{irat}.LFP.rename{chan_idx};
                %search channel into LFP data to remane it
                chan_idx = strcmp(LFP.label, config{irat}.LFP.allchannel{ichan});
                LFP.label{chan_idx} = new_name;
            end
        end  
   %% Preprocessing LFP for Wave detection 
    %remove breathing and ekg channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG','-Puff'};
        LFP             = ft_selectdata(cfgtemp, LFP);
        LFP_cleaned     = LFP; %save for later removing of artefacts
    
    %Separate trials 
        
    for itrial= 1: size(LFP_cleaned.trial,2)
        cfgtemp         = [];
        cfgtemp.trials = itrial;
        LFP_trial             = ft_selectdata(cfgtemp, LFP_cleaned);
        
        
        % Lowpass Filter LFP for detection
        cfgtemp             = [];
        cfgtemp.lpfilter    = 'yes';
        cfgtemp.lpfilttype  = 'fir';
        cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
        LFP_trial_lpfilt      = ft_preprocessing(cfgtemp, LFP_trial);
        
        %recover trial real timings to use it with muse markers
        starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
        endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
        offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
        
        
        
        %get hand-annotated wod timing
        wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);
        wod_delay  = MuseStruct{1}.markers.WOD.synctime(itrial) - MuseStruct{1}.markers.Vent_Off.synctime(itrial);
        
        %select times where to search WOD peak
        t = LFP_trial_lpfilt.time{1};
        t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
        t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
        t_sel = t_1 & t_2;
        
        [v_peak_wod, t_peak_wod] = findpeaks(-LFP_trial_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        wod_time_volt{irat}{itrial}= -v_peak_wod
        wod_time_time{irat}{itrial}= t_peak_wod
    end %itrial







end % irat