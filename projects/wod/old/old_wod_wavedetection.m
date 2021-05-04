function wod_wavedetection(configscript)

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = eval(configscript);

for irat= 4:size(config,2)
    %Load LFP and Muse markers
    temp= load(fullfile(config{irat}.datasavedir,sprintf('%s%s_%s.mat',config{irat}.prefix,'LFP',config{irat}.name{1})));
    LFP=temp.LFP{1,1}.WoD;
    clear temp
    MuseStruct               = readMuseMarkers(config{irat}, false);
    
    %vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
    startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
    if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
        error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
    end
    
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
    
    cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
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
            
            wod_marker = MuseStruct{1}{1}.markers.WOD.synctime(itrial);
            %select times where to search WOD peak
            t = LFP_lpfilt.time{itrial};
            t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
            t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial(itrial) + offsettrial(itrial));
            t_sel = t_1 & t_2;
            
            [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %WOR detection
            
            wor_marker = MuseStruct{1}{1}.markers.WOR.synctime(itrial);
            %select times where to search WOR peak
            t = LFP_lpfilt.time{itrial};
            t_1 = t > (wor_marker + config{irat}.LFP.wor_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
            t_2 = t < (wor_marker + config{irat}.LFP.wor_toisearch(2) - starttrial(itrial) + offsettrial(itrial));
            t_sel = t_1 & t_2;
            
            [v_peak_wor, t_peak_wor] = findpeaks(LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %store peak timings per channel in structure
            
            WOD_data.timings.(sprintf('Rat_%i',irat)).peak(ichan,itrial)= t_peak_wod;
            %express wor data compared to Vent_On
            t_VentOn= MuseStruct{1}{1}.markers.Vent_On.synctime(itrial)-starttrial(itrial) +offsettrial(itrial);
            WOR_data.timings.(sprintf('Rat_%i',irat)).peak(ichan,itrial)= t_peak_wor - t_VentOn;
            
            %plot detection for visual control
            fig_wodpeak= figure;
            plot(LFP_lpfilt.time{itrial},LFP_lpfilt.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_wod,-v_peak_wod,'x');
            xlim([t_peak_wod-10 t_peak_wod+10]);
            
            fig_worpeak= figure;
            plot(LFP_lpfilt.time{itrial},LFP_lpfilt.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_wor,v_peak_wor,'x');
            xlim([t_peak_wor-10 t_peak_wor+10]);
            
            detectsavedir=fullfile(config{irat}.imagesavedir,'detection');
            detectpeak_wod=fullfile(detectsavedir,'WoD','peak',sprintf('%s',config{irat}.prefix));
            detectpeak_wor=fullfile(detectsavedir,'WoR','peak',sprintf('%s',config{irat}.prefix));
            
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
            fname_worpeak=fullfile(detectpeak_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodpeak,fname_wodpeak,'png','pdf','close');
            dtx_savefigure(fig_worpeak,fname_worpeak,'png','pdf','close');
            
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
            WOD_data.timings.(sprintf('Rat_%i',irat)).min_slope(ichan,itrial)= t_peak_wodslope;
            WOD_data.values.(sprintf('Rat_%i',irat)).min_slope(ichan,itrial)= -v_peak_wodslope;
            
            %WOR threshold
            t1= t_peak_wor-30;
            t2= t_peak_wor+30;
            t_sel= [t1 t2];
            
            %cut data to keep only WOR
            cfgtemp=[];
            cfgtemp.latency= t_sel;
            WOR_cut= ft_selectdata(cfgtemp,LFP_lpfilt);
            clear t1 t2 t_sel
            
            %Transform into slope
            WOR_cut_slope=WOR_cut;
            WOR_cut_slope.trial{itrial}= ft_preproc_derivative(WOR_cut.trial{itrial});
            %smooth slope
            WOR_cut_slope.trial{itrial}= movmean(WOR_cut_slope.trial{itrial},1000,2);
            
            %Search for peaks in slope data
            %Determine time window to search
            t = WOR_cut.time{itrial};
            t_1 = t > (t_peak_wor - 3);
            t_2 = t < (t_peak_wor + 3);
            t_sel = t_1 & t_2;
            
            [v_peak_worslope, t_peak_worslope] = findpeaks(WOR_cut_slope.trial{itrial}(chan_idx,:),t,'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %save values
            %express timings of WoR according to Vent On
            t_VentOn= MuseStruct{1}{1}.markers.Vent_On.synctime(itrial)-starttrial(itrial) +offsettrial(itrial);
            WOR_data.timings.(sprintf('Rat_%i',irat)).min_slope(ichan,itrial)= t_peak_worslope - t_VentOn;
            
            WOR_data.values.(sprintf('Rat_%i',irat)).min_slope(ichan,itrial)= v_peak_worslope;
            
            %plot for visual control
            
            fig_wodslope=figure;
            plot(WOD_cut_slope.time{itrial},WOD_cut_slope.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_wodslope,-v_peak_wodslope,'x');
            xlim([t_peak_wod-10 t_peak_wod+10]);
            
            fig_worslope=figure;
            plot(WOR_cut_slope.time{itrial},WOR_cut_slope.trial{itrial}(ichan,:));
            hold on
            scatter(t_peak_worslope,v_peak_worslope,'x');
            xlim([t_peak_wor-10 t_peak_wor+10]);
            
            detectslope_wod=fullfile(detectsavedir,'WoD','slope',sprintf('%s',config{irat}.prefix));
            detectslope_wor=fullfile(detectsavedir,'WoR','slope',sprintf('%s',config{irat}.prefix));
            
            if ~isfolder(detectsavedir)
                mkdir(detectsavedir);
            end
            
            if ~isfolder(detectslope_wod)
                mkdir(detectslope_wod);
            end
            
            if ~isfolder(detectslope_wor)
                mkdir(detectslope_wor);
            end
            
            fname_wodslope=fullfile(detectslope_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            fname_worslope=fullfile(detectslope_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodslope,fname_wodslope,'png','pdf','close');
            dtx_savefigure(fig_worslope,fname_worslope,'png','pdf','close');
            
            clear fname_wodslope fname_worslope detectslope_wod detectslope_wor
            
            %% Determine threshold and crossing points
            
            %calculate threshold 20% of max slope
            wod_thr= -0.2*v_peak_wodslope;
            wor_thr= 0.2*v_peak_worslope;
            
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
            
            t = WOR_cut.time{itrial};
            t_1 = t > (t_peak_worslope - 10);
            t_2 = t < (t_peak_worslope + 10);
            t_sel = t_1 & t_2;
            
            %Create curve and horizontal line
            x1 = WOR_cut_slope.time{itrial}(1,t_sel);
            y1 = WOR_cut_slope.trial{itrial}(ichan,t_sel);
            x2 = x1;
            y2 = ones(size(y1)) * wor_thr;
            %Find values of intersection of 2 curves
            [x_worintersect, y_worintersect] = intersections(x1, y1, x2, y2);
            time_start_wor= x_worintersect(1);
            value_start_wor= y_worintersect(1);
            
            clear t t_1 t_2 t_sel
            
            %store values
            WOD_data.timings.(sprintf('Rat_%i',irat)).start(ichan,itrial)= time_start_wod;
            WOD_data.values.(sprintf('Rat_%i',irat)).slope_start(ichan,itrial)= value_start_wod;
            
            WOR_data.timings.(sprintf('Rat_%i',irat)).start(ichan,itrial)= time_start_wor-t_VentOn;
            WOR_data.values.(sprintf('Rat_%i',irat)).slope_start(ichan,itrial)= value_start_wod;
            
            %plot for visual control
            
            fig_wodthr= figure;
            plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
            xline(time_start_wod);
            
            fig_worthr=figure;
            plot(WOR_cut.time{itrial},WOR_cut.trial{itrial}(ichan,:));
            xline(time_start_wor);
            
            
            detectstart_wod=fullfile(detectsavedir,'WoD','start',sprintf('%s',config{irat}.prefix));
            detectstart_wor=fullfile(detectsavedir,'WoR','start',sprintf('%s',config{irat}.prefix));
            
            if ~isfolder(detectstart_wod)
                mkdir(detectstart_wod);
            end
            
            if ~isfolder(detectstart_wor)
                mkdir(detectstart_wor);
            end
            
            fname_wodstart=fullfile(detectstart_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            fname_worstart=fullfile(detectstart_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodthr,fname_wodstart,'png','pdf','close');
            dtx_savefigure(fig_worthr,fname_worstart,'png','pdf','close');
            
            
            clear x_wodintersect x_worintersect detectstart_wod detectstart_wor fname_wodstart fname_worstart
            %% Determine Half-width of waves
            
            %Calculate half amplitude of waves
            wod_amp= -v_peak_wod -y_wodintersect(1);
            wor_amp= v_peak_wor - y_worintersect(1);
            half_wod= wod_amp/2;
            half_wor= wor_amp/2;
            
            clear y_wodintersect y_worintersect x_wodintersect x_worintersect
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
            
            %WOR
            t = WOR_cut.time{itrial};
            t_1 = t > (t_peak_wor - 20);
            t_2 = t < (t_peak_wor + 20);
            t_sel = t_1 & t_2;
            
            x1 = WOR_cut.time{itrial}(1,t_sel);
            y1 = WOR_cut.trial{itrial}(ichan,t_sel);
            x2 = x1;
            y2 = ones(size(y1)) * half_wor;
            [x_worintersect, y_worintersect] = intersections(x1, y1, x2, y2);
            
            WOR_halfwi= x_worintersect(2)- x_worintersect(1);
            clear x1 y1 x2 y2
            
            %Store data
            WOD_data.values.(sprintf('Rat_%i',irat)).halfwidth(ichan,itrial)=WOD_halfwi;
            WOR_data.values.(sprintf('Rat_%i',irat)).halfwidth(ichan,itrial)=WOR_halfwi;
            
            %plot for visual control
            
            fig_wodhalf=figure;
            plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
            hold on
            scatter(x_wodintersect,y_wodintersect,'rx')
            yline(half_wod);
            
            
            fig_worhalf=figure;
            plot(WOR_cut.time{itrial},WOR_cut.trial{itrial}(ichan,:));
            hold on
            scatter(x_worintersect,y_worintersect,'rx')
            yline(half_wor);
            
            detecthalf_wod=fullfile(detectsavedir,'WoD','half-width',sprintf('%s',config{irat}.prefix));
            detecthalf_wor=fullfile(detectsavedir,'WoR','half-width',sprintf('%s',config{irat}.prefix));
            
            if ~isfolder(detecthalf_wod)
                mkdir(detecthalf_wod);
            end
            
            if ~isfolder(detecthalf_wor)
                mkdir(detecthalf_wor);
            end
            
            fname_wodhalf=fullfile(detecthalf_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            fname_worhalf=fullfile(detecthalf_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
            
            dtx_savefigure(fig_wodhalf,fname_wodhalf,'png','pdf','close');
            dtx_savefigure(fig_worhalf,fname_worhalf,'png','pdf','close');
            
            
            clear x_wodintersect x_worintersect y_wodintersect y_worintersect
            
            %% Create structure with electrode depths
            
            Electrode_depth.(sprintf('Rat_%i',irat))(ichan,itrial)=config{irat}.LFP.chan_depth{ichan};
            
            
            
        end %ichan
        
    end %itrial
end %irat

%% Arrange structures to have a version with all trials



%arrange structure to have all trials

depth_allrats = [];
icol = 0;
for idepth = 0:250:2500 %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for ratname = string(fieldnames(Electrode_depth))'
        for itrial = 1:size(Electrode_depth.(ratname), 2)
            irow = irow+1;
            sel = abs(Electrode_depth.(ratname)(:, itrial) - idepth) < 125;
            if sum(sel) == 1
                depth_allrats(irow,icol) = Electrode_depth.(ratname)(sel, itrial);
            elseif sum(sel) == 0
                depth_allrats(irow,icol) = nan;
            elseif sum(sel) > 1
                error('it should have only one electrode for one deepness');
            end
        end
    end
end

depth_allrats=depth_allrats';





%% Save structures

%save 'stats'

Detectionpath=fullfile(config{4}.datasavedir,'Detection');

if ~isfolder(Detectionpath)
    mkdir(Detectionpath);
end

save(fullfile(Detectionpath,'WoD_data.mat'),'WOD_data');
save(fullfile(Detectionpath,'WoR_data.mat'),'WOR_data');
save(fullfile(Detectionpath,'Depth_electrode'),'Electrode_depth');



save(fullfile(Detectionpath,'allrats_depth.mat'),'depth_allrats');
