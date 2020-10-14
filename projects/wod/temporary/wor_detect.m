%% WOR DATA : find WOR peak per channel, and normalize time per channel
                    %Antoine
                    %+ normaliser le temps que selon la 1ere wod/wor, pour toutes les électrodes du rat ? pour
                    %mieux voir le délai entre électrodes
                    for ichan= 1 : size(LFP_trial_filt.label,1)
                    for ichan_name= LFP_trial_filt.label{ichan}
                    timefreq_wor{ifreq}{itrial}.(ichan_name)            = timefreq_ichan_temp;
                    
                    %select lfp channel with the same name as tfr channel (in
                    %case channel numbers were schuffled by fieldtrip)
                    chan_idx    = strcmp(LFP_trial_filt.label, ichan_name);
                    
                    %get hand-annotated wod timing
                    wor_marker = MuseStruct{1}{1}.markers.WOR.synctime(itrial);
                    
                    %select times where to search WOD peak 
                    t = LFP_trial_filt.time{1};
                    t_1 = t > (wor_marker + config{irat}.LFP.wor_toisearch(1) - starttrial + offsettrial);
                    t_2 = t < (wor_marker + config{irat}.LFP.wor_toisearch(2) - starttrial + offsettrial);
                    t_sel = t_1 & t_2;
                    %Search LFP maximum peak in this selected window. '-'LFP because wod is negative 
                    [~, t_peak_wor] = findpeaks(LFP_trial_filt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                    end_wor= LFP.trialinfo.endsample(itrial) / LFP.fsample;
                    %keep only data between wod and end
                    cfgtemp                                   = [];
                    cfgtemp.latency                           = [t_peak_wor end_wor];
                    timefreq_wor{ifreq}{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wor{ifreq}{itrial}.(ichan_name));
                    
                    %normalize time
                    timefreq_wor_timenorm{ifreq}{itrial}.(ichan_name)        = timefreq_wor{ifreq}{itrial}.(ichan_name);
                    timefreq_wor_timenorm{ifreq}{itrial}.(ichan_name).time   = timefreq_wor{ifreq}{itrial}.(ichan_name).time ./ t_peak_wod;
                    %check the location of the peak detection
                                    figure;hold
                                    plot(LFP_trial_filt.time{1}/t_peak_wor, LFP_trial_filt.trial{1}(ichan,:));
                                    xlim([0 2]);
                                    xlim([0.8 1.2]);
                    end
                    end
                    
                    
                    %resample data to have the same number of data points, for
                    %time-normalized data
                    t_old                                                = timefreq_wor_timenorm{ifreq}{itrial}.(ichan_name).time;
                    t_new                                                = linspace(0,1,1000);
                    powspctrm_new(1,1,:)                                 = pchip(t_old,squeeze(timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
                    timefreq_wor_timenorm{ifreq}{itrial}.(ichan_name).time         = t_new;
                    timefreq_wor_timenorm{ifreq}{itrial}.(ichan_name).powspctrm    = powspctrm_new;
                    
                   % plot(timefreq_wor{ifreq}{itrial}.time{1},timefreq_wor{ifreq}{itrial}.trial{1}); xlim([0 0.95])
                    