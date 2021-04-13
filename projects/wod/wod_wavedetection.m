function wod_grandaverage(configscript)

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

config = eval('wod_setparams');

for irat= rat_list
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
        
        %select one trial
        cfgtemp         = [];
        cfgtemp.trials  = itrial;
        LFP_trial       = ft_selectdata(cfgtemp, LFP_cleaned);
        
        %recover trial real timings to use it with muse markers
        starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
        endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
        offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
        
        
        
        for ichan= 1:size(LFP.label,1)
            
            ichan_name              = LFP_trial.label{ichan};
            
            %% WOD and WOR peak detection
            
            %WOD detection
            %select lfp channel (in
            %case channel numbers were schuffled by fieldtrip)
            chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
            
            wod_marker = MuseStruct{1}{1}.markers.WOD.synctime(itrial);
            %select times where to search WOD peak
            t = LFP_lpfilt.time{1};
            t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
            t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
            t_sel = t_1 & t_2;
            
            [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %WOR detection
            
            wor_marker = MuseStruct{1}{1}.markers.WOR.synctime(itrial);
            %select times where to search WOR peak
            t = LFP_lpfilt.time{1};
            t_1 = t > (wor_marker + config{irat}.LFP.wor_toisearch(1) - starttrial + offsettrial);
            t_2 = t < (wor_marker + config{irat}.LFP.wor_toisearch(2) - starttrial + offsettrial);
            t_sel = t_1 & t_2;
            
            [v_peak_wor, t_peak_wor] = findpeaks(LFP_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
            clear t t_1 t_2 t_sel
            
            %store peak timings per channel in structure
            
            WOD_data.timings.(irat).peak.(ichan,itrial)= t_peak_wod;
            WOD_data.timings.(irat).peak.(ichan,itrial+2)=ichan_name;
            WOR_data.timings.(irat).peak.(ichan,trial)= t_peak_wor;
            WOR_data.timings.(irat).peak.(ichan,itrial+2)=ichan_name;
            
            
        end %ichan
        
        %% Determine threshold
            
            %WOD threshold
            t1= MuseStruct{1}{1}.markers.WOR.synctime(itrial)-30;
            t2= MuseStruct{1}{1}.markers.WOR.synctime(itrial)+30;
            t_sel= [t1 t2];
            
            %cut data to keep only WOD
            cfgtemp=[];
            cfgtemp.latency= t_sel;
            WOD_data= ft_selectdata(cfgtemp,LFP_lpfilt);
            
            

           
            
        
    end %itrial
end %irat