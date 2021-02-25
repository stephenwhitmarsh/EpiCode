function wod_short(slurm_task_id)

% time frequency analysis of WOD Neuralynx data for base line and WoD data
% parameters are set in wod_setparams.m

% slurm_task_id : input integer value used by the slurm script, to
% parallelise computation between rats.
% slurm_task_id is the number of the rat, as set in wod_setparams.m
% to compute/plot the average between rats, set slurm_task_id = 0

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


ipart = 1; %ipart is always 1 for this project

if slurm_task_id(1) > 0
    for irat = slurm_task_id
        [~,dir_name]                       = fileparts(config{irat}.rawdir);
        config{irat}.rawdir                = fullfile(config{irat}.datasavedir,'concatenated_LFP');
        config{irat}.directorylist{ipart}  = {dir_name};
        
         %read Muse markers
        MuseStruct               = readMuseMarkers(config{irat}, true);
        
        LFP = readLFP(config{irat}, MuseStruct, true);
        LFP = LFP{1}.(config{irat}.LFP.name{1});
        
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
        
        %analyse each trial independently
        for itrial = 1:size(LFP.trial,2)
            
            %select one trial
            cfgtemp         = [];
            cfgtemp.trials  = itrial;
            LFP_trial       = ft_selectdata(cfgtemp, LFP);
            
            %filter lfp to better recognize WOD/WOR peak
            cfgtemp             = [];           
            cfgtemp.lpfilter    = 'yes';
            cfgtemp.lpfilttype  = 'fir';
            cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
            LFP_trial_lpfilt      = ft_preprocessing(cfgtemp, LFP_trial);

            
            %hp filter lfp to exclude WoD and WoR
            cfgtemp             = [];           
            cfgtemp.hpfilter    = 'yes';
            cfgtemp.hpfilttype  = 'fir';
            cfgtemp.hpfreq      = config{irat}.LFP.hpfilter_wod_exclusion;
            LFP_trial           = ft_preprocessing(cfgtemp, LFP_trial);
            
            %recover trial real timings to use it with muse markers
            starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
            endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
            offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
            
            %%find index of t(wod) in LFP data
            Wod_indx= find(LFP_trial.time{1}== MuseStruct{1}{1}.markers.WOD.synctime(itrial));
            
            for iparam =['long' 'short']
            %do time frequency analysis
            cfgtemp                         = [];
            cfgtemp.channel                 = 'all';
            cfgtemp.method                  = 'mtmconvol';
            cfgtemp.output                  = 'pow';
            cfgtemp.taper                   = 'dpss'; %default = dpss
            cfgtemp.pad                     = 'nextpow2';
            cfgtemp.keeptrials              = 'no';
            cfgtemp.tapsmofrq               = config{irat}.timefreq.tapsmofrq_short;
            cfgtemp.foi                     = config{irat}.timefreq.foi;
            cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*config{irat}.timefreq.t_ftimwin_short;
            cfgtemp.toi                     = LFP_trial.time{1}(1) : config{irat}.timefreq.timestep : LFP_trial.time{1}(end);
            timefreq_alldata{itrial}{iparam}        = ft_freqanalysis(cfgtemp,LFP_trial);
                            %% DETECT WOD AND SELECT DATA AS BASELINE AND WOD
            for ichan = 1:size(LFP_trial_lpfilt.label,1)
                %select channel
                ichan_name              = LFP_trial_lpfilt.label{ichan};
                cfgtemp                 = [];
                cfgtemp.channel         = ichan_name;
                timefreq_ichan_temp   	= ft_selectdata(cfgtemp,LFP_trial_lpfilt.trial{ichan});
                
                %% WOD DATA : find WOD peak per channel, and normalize time per channel
                  %use filtered data to find wod
                  %select lfp channel with the same name as tfr channel (in
                  %case channel numbers were schuffled by fieldtrip)
                  chan_idx    = strcmp(LFP_trial_lpfilt.label, ichan_name);
                  
                   %get hand-annotated wod timing
                wod_marker = MuseStruct{1}{1}.markers.WOD.synctime(itrial);
                
                %select times where to search WOD peak
                t = LFP_trial_lpfilt.time{1};
                t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
                t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
                t_sel = t_1 & t_2;
                %Search LFP maximum peak in this selected window. '-'LFP because wod is negative
                [v_peak_wod, t_peak_wod] = findpeaks(-LFP_trial_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                clear t t_1 t_2 t_sel
                
                %selectime time between 0 and t(WoD)
                cfgtemp                                   = [];
                cfgtemp.latency                           = [0 t_peak_wod];
                LFP_trial_wod.(ichan_name)  = ft_selectdata(cfgtemp,LFP_trial.trial{ichan});
                
               
            end% ichan
        end% itrial
    end %irat
end % slurm task id