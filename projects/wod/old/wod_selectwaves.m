if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\Antoine\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/Antoine/EpiCode/projects/dtx'))
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
        
        %plot to control filtering
%         plot(LFP_trial_lpfilt.time{1}(1:640),LFP_trial_lpfilt.trial{1}(1,1:640))
%         plot(LFP.time{1}(1:640),LFP.trial{1}(1,1:640))
%         close all
        
        
        
        %recover trial real timings to use it with muse markers
        starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
        endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
        offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
        
        
        %% WOD detection
        
        for ichan= 1:size(LFP_cleaned.label,1)
            
        chan_name= LFP_cleaned.label{ichan};
        
        
        %select only one channel
        cfgtemp=[];
        cfgtemp.channel= chan_name;
        chandata_filt= ft_selectdata(cfgtemp, LFP_trial_lpfilt);
       
        %get hand-annotated wod timing
        wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);
        wod_delay  = MuseStruct{1}.markers.WOD.synctime(itrial) - MuseStruct{1}.markers.Vent_Off.synctime(itrial);
        
        %select times where to search WOD peak
        t = LFP_trial_lpfilt.time{1};
        t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
        t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
        t_sel = t_1 & t_2;
        
        [v_peak_wod, t_peak_wod] = findpeaks(-chandata_filt.trial{1}(t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        
        wod_timings{itrial}(ichan)= t_peak_wod; % to save by rat
        
        clear t_1 t_2 t_sel
        
        %% WOR detection
        
        %get hand-annotated wod timing
        wor_marker = MuseStruct{1}.markers.WOR.synctime(itrial);
        wor_delay  = MuseStruct{1}.markers.WOR.synctime(itrial) - MuseStruct{1}.markers.Vent_Off.synctime(itrial);
        
        t_1 = t > (wor_marker + config{irat}.LFP.wor_toisearch(1) - starttrial + offsettrial);
        t_2 = t < (wor_marker + config{irat}.LFP.wor_toisearch(2) - starttrial + offsettrial);
        t_sel = t_1 & t_2;
        
        [v_peak_wor, t_peak_wor] = findpeaks(chandata_filt.trial{1}(t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        wor_timings{itrial}(ichan) = t_peak_wor;
        
        clear t_1 t_2 t_sel
        
        end %ichan
        
        %% Cutting the data
        % use to cut the LFP data
        t_minwod= min(wod_timings{itrial});
        t_maxwod= max(wod_timings{itrial});
        t_minwor= min(wor_timings{itrial});
        t_maxwor= max(wor_timings{itrial});
        
        
        %cut LFP for WOD
        cfgtemp= [];
        cfgtemp.latency= [t_minwod-30 t_maxwod+30];
        cfgtemp.channels= 'all';
        LFP_wodcut= ft_selectdata(cfgtemp, LFP_trial);
        plot(LFP_wodcut.time{1},LFP_wodcut.trial{1})
        %cut LFP for WOD
        cfgtemp= [];
        cfgtemp.latency= [t_minwor-30 t_maxwor+30];
        cfgtemp.channels= 'all';
        LFP_worcut= ft_selectdata(cfgtemp, LFP_trial);
        plot(LFP_worcut.time{1},LFP_worcut.trial{1})

        %save LFP and wave timings
        
        %Save LFP
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'WOD_cut',sprintf('%d',itrial),'.mat']),'LFP_wodcut','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'WOR_cut',sprintf('%d',itrial),'.mat']),'LFP_worcut','-v7.3');
        
        %Save wave timings
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'WOD_time',sprintf('%d',itrial),'.mat']),'wod_timings','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix,'WOR_time',sprintf('%d',itrial),'.mat']),'wor_timings','-v7.3');
                

        
    end %itrial
end % irat