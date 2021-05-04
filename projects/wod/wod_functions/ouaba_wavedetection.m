function stats = ouaba_wavedetection(cfg, force)

fname_out = fullfile(cfg{4}.datasavedir,'Detection', sprintf('ouaba_wavedetection_allrats.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'stats_all');
    return
end

for irat= 1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    
    iratname= sprintf('Rat_%i',irat);
    
    %Load LFP and Muse markers
    temp= load(fullfile(cfg{irat}.datasavedir,sprintf('%s%s_%s.mat',cfg{irat}.prefix,'LFP',cfg{irat}.name{1})));
    LFP=temp.LFP{1,1}.WoD;
    clear temp
    MuseStruct               = readMuseMarkers(cfg{irat}, false);
    
    %vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
    startmarker = cfg{irat}.muse.startmarker.(cfg{irat}.LFP.name{1});
    if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
        error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', cfg{irat}.prefix(1:end-1));
    end
 
    %rename channels according to depth
    for ichan = 1:size(cfg{irat}.LFP.channel, 2)
        idx = strcmp(cfg{irat}.LFP.channel{ichan}, LFP.label);
        label_renamed{idx} = cfg{irat}.LFP.rename{ichan};
    end
    LFP.label = label_renamed';
    
    
    %remove breathing and ekg channel
    cfgtemp         = [];
    cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG'};
    LFP             = ft_selectdata(cfgtemp, LFP);
    LFP_cleaned     = LFP; %save for later removing of artefacts
    
    %remove 50Hz and interpolate with 49 and 51 Hz
    %LFP_cleaned= ft_preproc_dftfilter(LFP_cleaned,LFP_cleaned.fsample,50,'Flreplace','neighbour');
    
    
    %filter lfp to better recognize WOD/WOR peak
    cfgtemp             = [];
    cfgtemp.lpfilter    = 'yes';
    cfgtemp.lpfilttype  = 'fir';
    
    cfgtemp.lpfreq      = cfg{irat}.LFP.lpfilter_wod_detection;
    LFP_lpfilt      = ft_preprocessing(cfgtemp, LFP_cleaned);
    
    for itrial = 1:size(LFP.trial,2)
        
        
        %recover trial real timings to use it with muse markers
        starttrial              = LFP_lpfilt.trialinfo.begsample / LFP_lpfilt.fsample;
        endtrial                = LFP_lpfilt.trialinfo.endsample / LFP_lpfilt.fsample;
        offsettrial             = LFP_lpfilt.trialinfo.offset / LFP_lpfilt.fsample;
        
        
        
        for ichan= 1:size(LFP.label,1)
            
            ichan_name              = LFP_lpfilt.label{ichan};
            
            %% WOD and WOR peak detection
            
            %WOD detection
            %select lfp channel (in
            %case channel numbers were schuffled by fieldtrip)
            chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
            
            wod_marker = MuseStruct{1}{1}.markers.AD__START__.synctime(itrial);
            %select times where to search WOD peak
            t = LFP_lpfilt.time{itrial};
            t_1 = t > (wod_marker + cfg{irat}.LFP.wod_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
            t_2 = t < (wod_marker + cfg{irat}.LFP.wod_toisearch(2) - starttrial(itrial) + offsettrial(itrial));
            t_sel = t_1 & t_2;
            
            [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            
            %store peak timings per channel in structure
            stats_all{irat}.WoD.peak_time(ichan,itrial)= t_peak_wod
            
            %plot detection for visual control
            fig_wodpeak= figure;
            plot(LFP_lpfilt.time{itrial},LFP_lpfilt.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_wod,-v_peak_wod,'x');
            xlim([t_peak_wod-10 t_peak_wod+10]);
            
           
            detectsavedir=fullfile(cfg{irat}.imagesavedir,'detection');
            detectpeak_wod=fullfile(detectsavedir,'WoD','peak',sprintf('%s',cfg{irat}.prefix));
            
            if ~isfolder(detectsavedir)
                mkdir(detectsavedir);
            end
            
            if ~isfolder(detectpeak_wod)
                mkdir(detectpeak_wod);
            end
            
            if ~isfolder(detectpeak_wor)
                mkdir(detectpeak_wor);
            end
            
            fname_wodpeak=fullfile(detectpeak_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodpeak,fname_wodpeak,'png','pdf','close');
            
            clear fname_wodpeak fname_worpeak detectpeak_wod detectpeak_wor
            %% Determine minimum and maximum slopes and extract timings and values
            
            %WOD window selection
            t1= t_peak_wod-30;
            t2= t_peak_wod+30;
            t_sel= [t1 t2];
            
            %cut data to keep only WOD
            cfgtemp=[];
            cfgtemp.latency= t_sel;
            WOD_cut= ft_selectdata(cfgtemp,LFP_lpfilt);
            clear t1 t2 t_sel
            
            %Transform into slope
            WOD_cut_slope=WOD_cut;
            WOD_cut_slope.trial{itrial}= ft_preproc_derivative(WOD_cut.trial{itrial});
            %smooth slope
            WOD_cut_slope.trial{itrial}= movmean(WOD_cut_slope.trial{itrial},100,2);
            
            %Search for peaks in slope data
            %Determine time window to search
            t = WOD_cut.time{itrial};
            t_1 = t > (t_peak_wod - 10);
            t_2 = t < (t_peak_wod + 10);
            t_sel = t_1 & t_2;
            
            [v_peak_wodslope, t_peak_wodslope] = findpeaks(-WOD_cut_slope.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %save values
            
            stats_all{irat}.WoD.min_slope_time(ichan,itrial)=  t_peak_wodslope;
            stats_all{irat}.WoD.min_slope_value(ichan,itrial)=   -v_peak_wodslope;

            %plot for visual control
            
            fig_wodslope=figure;
            plot(WOD_cut_slope.time{itrial},WOD_cut_slope.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_wodslope,-v_peak_wodslope,'x');
            xlim([t_peak_wod-10 t_peak_wod+10]);
        
            detectslope_wod=fullfile(detectsavedir,'WoD','slope',sprintf('%s',cfg{irat}.prefix));
            
            if ~isfolder(detectsavedir)
                mkdir(detectsavedir);
            end
            
            if ~isfolder(detectslope_wod)
                mkdir(detectslope_wod);
            end
            
          
            fname_wodslope=fullfile(detectslope_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
           
            dtx_savefigure(fig_wodslope,fname_wodslope,'png','pdf','close');
           
            clear fname_wodslope detectslope_wod 
            
            %% Determine threshold and crossing points
            
            %calculate threshold 20% of max slope
            wod_thr= -0.2*v_peak_wodslope;
          
            %Determine crossing point
            %define time window
            t = WOD_cut.time{itrial};
            t_1 = t > (t_peak_wodslope - 5);
            t_2 = t < (t_peak_wodslope + 5);
            t_sel = t_1 & t_2;
            
            %Create curve and horizontal line
            x1 = WOD_cut_slope.time{itrial}(1,t_sel);
            y1 = WOD_cut_slope.trial{itrial}(ichan,t_sel);
            x2 = x1;
            y2 = ones(size(y1)) * wod_thr;
            %Find values of intersection of 2 curves
            [x_wodintersect, y_wodintersect] = intersections(x1, y1, x2, y2);
            
            time_start_wod= x_wodintersect(1);
            value_start_wod= y_wodintersect(1);
            
            clear t t_1 t_2 t_sel
            
           
            %store values
            
            stats_all{irat}.WoD.start_time(ichan,itrial)=time_start_wod;
            stats_all{irat}.WoD.start_slope_value(ichan,itrial)=value_start_wod;
            
            %plot for visual control
            
            fig_wodthr= figure;
            plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
            xline(time_start_wod);
                                   
            detectstart_wod=fullfile(detectsavedir,'WoD','start',sprintf('%s',cfg{irat}.prefix));
            
            if ~isfolder(detectstart_wod)
                mkdir(detectstart_wod);
            end
            
         
            fname_wodstart=fullfile(detectstart_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodthr,fname_wodstart,'png','pdf','close');
            
            
            clear x_wodintersect  detectstart_wod  fname_wodstart 
            %% Determine Half-width of waves
            
            %Calculate half amplitude of waves
            wod_amp= -v_peak_wod -y_wodintersect(1);
            half_wod= wod_amp/2;
            
            clear y_wodintersect  x_wodintersect 
            %Determine time window to search
            %WOD
            
            t = WOD_cut.time{itrial};
            t_1 = t > (t_peak_wod - 10);
            t_2 = t < (t_peak_wod + 10);
            t_sel = t_1 & t_2;
            
            x1 = WOD_cut.time{itrial}(1,t_sel);
            y1 = WOD_cut.trial{itrial}(ichan,t_sel);
            x2 = x1;
            y2 = ones(size(y1)) * half_wod;
            [x_wodintersect, y_wodintersect] = intersections(x1, y1, x2, y2);
            
            WOD_halfwi= x_wodintersect(2)- x_wodintersect(1);
            clear x1 y1 x2 y2
            
            %Store data
            
            stats_all{irat}.WoD.half_width(ichan,itrial)=WOD_halfwi;

            plot for visual control
            
            fig_wodhalf=figure;
            plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
            hold on
            scatter(x_wodintersect,y_wodintersect,'rx')
            yline(half_wod);

            
            detecthalf_wod=fullfile(detectsavedir,'WoD','half-width',sprintf('%s',cfg{irat}.prefix));
            
            if ~isfolder(detecthalf_wod)
                mkdir(detecthalf_wod);
            end

            fname_wodhalf=fullfile(detecthalf_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodhalf,fname_wodhalf,'png','pdf','close');
            
            
            clear x_wodintersect  y_wodintersect 
            
            %% Create structure with electrode depths

            stats_all{irat}.Depth(ichan,itrial)=cfg{irat}.LFP.chan_depth{ichan};
        end %ichan
        
    end %itrial
end %irat


%% Save structures

save(fname_out, 'stats_all');

