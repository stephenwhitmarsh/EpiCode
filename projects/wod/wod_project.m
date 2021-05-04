function wod_project(rat_list, configscript)

% time frequency analysis of WOD Neuralynx data
% parameters are set in wod_setparams.m

% slurm_task_id : input integer value used by the slurm script, to
% parallelise computation between rats.
% slurm_task_id is the number of the rat, as set in wod_setparams.m
% to compute/plot the average between rats, set slurm_task_id = 0

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = eval(configscript);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute data for each rat %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ipart = 1; %ipart is always 1 for this project

    for irat = rat_list
        %% load data and compute time freq
        %find concatenated LFP (see wod_concatenateLFP.m)
        [~,dir_name]                       = fileparts(config{irat}.rawdir);
        config{irat}.rawdir                = fullfile(config{irat}.concatdata_path,config{irat}.prefix);
        config{irat}.directorylist{ipart}  = config{irat}.prefix;%{dir_name};
        
        %read Muse markers
        MuseStruct               = readMuseMarkers(config{irat}, true);
        
        %read LFP. T0 = Vent_Off. Each trial is one protocol
        %if exist(name_ica, 'file')
        %   load(name_ica)
        %else
        LFP = readLFP(config{irat}, MuseStruct, true);
        LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
        %end
        
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
        %             plot(LFP_trial_lpfilt.time{1}(1:640),LFP_trial_lpfilt.trial{1}(1,1:640))
        %             plot(LFP.time{1}(1:640),LFP.trial{1}(1,1:640))
        %             close all
        
        
        %hp filter lfp to exclude WoD and WoR and Notch 50Hz
        cfgtemp             = [];
        cfgtemp.hpfilter    = 'yes';
        cfgtemp.hpfilttype  = 'fir';
        cfgtemp.hpfreq      = config{irat}.LFP.hpfilter_wod_exclusion;
        cfgtemp.bsfilter     = 'yes';
        cfgtemp.bsfilttype     = 'fir';
        cfgtemp.bsfreq          = [49 51];
        LFP_cleaned          = ft_preprocessing(cfgtemp, LFP_cleaned);
        
        
        %analyse each trial independently
        for itrial = 1:size(LFP.trial,2)
            
            %select one trial
            cfgtemp         = [];
            cfgtemp.trials  = itrial;
            LFP_trial       = ft_selectdata(cfgtemp, LFP_cleaned);
            
            %recover trial real timings to use it with muse markers
            starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;
            endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
            offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
            
            
            for iparam= ["long" "short"]
                
                
                %do time frequency analysis
                cfgtemp                         = [];
                cfgtemp.channel                 = 'all';
                cfgtemp.method                  = 'mtmconvol';
                cfgtemp.output                  = 'pow';
                cfgtemp.taper                   = 'dpss'; %default = dpss
                cfgtemp.pad                     = 'nextpow2';
                cfgtemp.keeptrials              = 'no';
                cfgtemp.tapsmofrq               = config{irat}.timefreq.tapsmofrq.(iparam);
                cfgtemp.foi                     = config{irat}.timefreq.foi;
                cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*config{irat}.timefreq.t_ftimwin.(iparam);
                cfgtemp.toi                     = config{irat}.timefreq.toi.(iparam)(1) : config{irat}.timefreq.timestep.(iparam) : config{irat}.timefreq.toi.(iparam)(end);
                timefreq_alldata{itrial}.(iparam) = ft_freqanalysis(cfgtemp,LFP_trial);
                
                %replace artifacts by nans
                %need to remove artefacts after time freq analysis, because
                %any nan in the lfp data creates a time freq with only nan values
                if isfield(MuseStruct{1}{1}.markers, 'BAD__START__')
                    if isfield(MuseStruct{1}{1}.markers.BAD__START__, 'synctime')
                        %get bad timings
                        bad_start                   = MuseStruct{1}{1}.markers.BAD__START__.synctime;
                        bad_end                     = MuseStruct{1}{1}.markers.BAD__END__.synctime;
                        if length(bad_start) ~= length(bad_end)
                            error('not the same amount of bad start and end markers');
                        end
                        %                         t_tfr                   = timefreq_alldata{itrial}.(iparam).time;
                        t_lfp                   = LFP_trial.time{1};
                        t_tfr                   = timefreq_alldata{itrial}.(iparam).time;
                        bad_sel                 = find(bad_start >= starttrial & bad_start <= endtrial);
                        %go through each bad timing
                        for ibad = bad_sel
                            %remove lfp artefacts
                            bad_period_lfp = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) & t_lfp <= (bad_end(ibad)- starttrial + offsettrial);
                            LFP_cleaned.trial{itrial}(:,bad_period_lfp) = nan(size(LFP_cleaned.trial{itrial},1),sum(bad_period_lfp));
                            %remove tfr artefacts: all window with at least one artefacted sample
                            bad_period_tfr = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - config{irat}.timefreq.t_ftimwin.(iparam)/2 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + config{irat}.timefreq.t_ftimwin.(iparam)/2;
                            timefreq_alldata{itrial}.(iparam).powspctrm(:,:,bad_period_tfr) = nan(size(timefreq_alldata{itrial}.(iparam).powspctrm,1),size(timefreq_alldata{itrial}.(iparam).powspctrm,2),sum(bad_period_tfr));
                        end
                        clear bad_sel t_tfr t_lfp
                    end
                end
                
            end %iparam
            
            for ichan = 1:size(timefreq_alldata{itrial}.long.label,1)
                
                %select channel for long param
                ichan_name              = timefreq_alldata{itrial}.long.label{ichan};
                cfgtemp                 = [];
                cfgtemp.channel         = ichan_name;
                timefreq_ichan_temp   	= ft_selectdata(cfgtemp,timefreq_alldata{itrial}.long);
                
                %% RECOVERY DATA : Change T0 from Vent_Off to Vent_On
                %MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_On.synctime(itrial);
                timeshift                                           = config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2) + config{irat}.epoch.pad.(config{irat}.LFP.name{1}) - timefreq_ichan_temp.time(end);%MuseStruct{1}{1}.markers.Vent_On.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial);
                timefreq_recovery{itrial}.(ichan_name)       = timefreq_ichan_temp;
                timefreq_recovery{itrial}.(ichan_name).time 	= timefreq_ichan_temp.time + timeshift;
                
                cfgtemp = [];
                cfgtemp.latency = [0 config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2)];
                timefreq_recovery{itrial}.(ichan_name) = ft_selectdata(cfgtemp,timefreq_recovery{itrial}.(ichan_name));
                
                
                %% WOR DATA : find WOR peak per channel, and normalize time per channel
                %Antoine
                %+ normaliser le temps que selon la 1ere wod/wor, pour toutes les électrodes du rat ? pour
                %mieux voir le délai entre électrodes
                
                timefreq_wor{itrial}.(ichan_name)            = timefreq_ichan_temp;
                
                %select lfp channel with the same name as tfr channel (in
                %case channel numbers were schuffled by fieldtrip)
                chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
                
                %get hand-annotated wod timing
                wor_marker = MuseStruct{1}{1}.markers.WOR.synctime(itrial);
                
                %select times where to search WOR peak
                t = LFP_lpfilt.time{1};
                t_1 = t > (wor_marker + config{irat}.LFP.wor_toisearch(1) - starttrial + offsettrial);
                t_2 = t < (wor_marker + config{irat}.LFP.wor_toisearch(2) - starttrial + offsettrial);
                t_sel = t_1 & t_2;
                %Search LFP maximum peak in this selected window. wor is positive
                [v_peak_wor, t_peak_wor] = findpeaks(LFP_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                end_wor= timefreq_wor{itrial}.(ichan_name).time(end);
                clear t t_1 t_2 t_sel
                
                %keep only data between wor and end
                cfgtemp                                   = [];
                cfgtemp.latency                           = [t_peak_wor end_wor];
                timefreq_wor{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wor{itrial}.(ichan_name));
                
                
                
                %normalize time
                timefreq_wor_timenorm{itrial}.(ichan_name)        = timefreq_wor{itrial}.(ichan_name);
                timefreq_wor_timenorm{itrial}.(ichan_name).time   = timefreq_wor{itrial}.(ichan_name).time-t_peak_wor;%./end_wor;
                timefreq_wor_timenorm{itrial}.(ichan_name).time   = timefreq_wor_timenorm{itrial}.(ichan_name).time./timefreq_wor_timenorm{itrial}.(ichan_name).time(end);%./end_wor;
                %check the location of the peak detection
                fig=figure;hold
                plot(LFP_lpfilt.time{1}, LFP_lpfilt.trial{1}(ichan,:));
                scatter(t_peak_wor, v_peak_wor,50,'red');
                
                xlim([t_peak_wor-10 t_peak_wor+10]);
                print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'Detection',sprintf('Rat%g_WOR%g_%s.png',irat,itrial,ichan_name)),'-r600');
                close all
                %
                for ifreq = 1:size(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm,2)
                    %resample data to have the same number of data points, for
                    %time-normalized data
                    t_old                                                = timefreq_wor_timenorm{itrial}.(ichan_name).time;
                    t_new                                                = linspace(0,1,10000);
                    %                     powspctrm_new(1,1,:)                                 = pchip(t_old,squeeze(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
                    powspctrm_new(1,ifreq,:)                                 = interp1(t_old,squeeze(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:)),t_new);%pchip(t_old,squeeze(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
                end
                timefreq_wor_timenorm{itrial}.(ichan_name).time         = t_new;
                timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm    = powspctrm_new;
                clear powspctrm_new
                
                
               
                %% baseline long
                %VOIR COMMENT
                
                cfgtemp = [];
                cfgtemp.latency = [-300 0];
                timefreq_baseline_long{itrial}.(ichan_name) = ft_selectdata(cfgtemp,timefreq_ichan_temp);
                
                
                %% MAKE BASELINE CORRECTION FOR LONG PERIODS
                
                %duplicate timefreq data
                timefreq_recovery_blcorrected{itrial}.(ichan_name)= timefreq_recovery{itrial}.(ichan_name);
                timefreq_wor_blcorrected{itrial}.(ichan_name)= timefreq_wor{itrial}.(ichan_name);
                timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name)= timefreq_wor_timenorm{itrial}.(ichan_name);
                
                
                %average power by frequency step and apply baseline
                %correction
                for ifreq = 1:size(timefreq_recovery{itrial}.(ichan_name).freq,2) %baseline
                    baseline_ifreq = nanmean(squeeze(timefreq_baseline_long{itrial}.(ichan_name).powspctrm(1,ifreq,:))); %baseline
                    timefreq_recovery_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_recovery{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq; %baseline
                    timefreq_wor_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wor{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
                    timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
                end %baseline
                
                clear timefreq_ichan_temp timefreq_baseline_long
                
                size(log10(squeeze(timefreq_recovery_blcorrected{itrial}.(ichan_name).powspctrm(1,1,:))))
                size(log10(timefreq_recovery_blcorrected{itrial}.(ichan_name).powspctrm(1,1,:)))
                %% Make Logarythm for long periods
                
                %duplicate fieldtrip structure
                log_timefreq_recovery{itrial}.(ichan_name)=timefreq_recovery{itrial}.(ichan_name);
                log_timefreq_recovery_blcorrected{itrial}.(ichan_name) = timefreq_recovery{itrial}.(ichan_name);
                log_timefreq_wor{itrial}.(ichan_name)=timefreq_wor{itrial}.(ichan_name);
                log_timefreq_wor_blcorrected{itrial}.(ichan_name)=timefreq_wor_blcorrected{itrial}.(ichan_name);
                log_timefreq_wor_timenorm{itrial}.(ichan_name)=timefreq_wor_timenorm{itrial}.(ichan_name);
                log_timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name)=timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name);
                
                
                %apply logarythm to powspctrm
                log_timefreq_recovery{itrial}.(ichan_name).powspctrm= log10(timefreq_recovery{itrial}.(ichan_name).powspctrm);
                log_timefreq_recovery_blcorrected{itrial}.(ichan_name).powspctrm =log10(timefreq_recovery_blcorrected{itrial}.(ichan_name).powspctrm) ;
                              
                log_timefreq_wor{itrial}.(ichan_name).powspctrm= log10(timefreq_wor{itrial}.(ichan_name).powspctrm);
                log_timefreq_wor_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wor_blcorrected{itrial}.(ichan_name).powspctrm);
                log_timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm= log10(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm);
                log_timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name).powspctrm);
                
                
                %% select channel for short param
                ichan_name              = timefreq_alldata{itrial}.short.label{ichan};
                cfgtemp                 = [];
                cfgtemp.channel         = ichan_name;
                timefreq_ichan_temp   	= ft_selectdata(cfgtemp,timefreq_alldata{itrial}.short);
                
                
                %% WOD DATA : find WOD peak per channel, and normalize time per channel
                %use filtered data to find wod
                
                timefreq_wod{itrial}.(ichan_name)            = timefreq_ichan_temp;
                
                %select lfp channel with the same name as tfr channel (in
                %case channel numbers were schuffled by fieldtrip)
                chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
                
                %get hand-annotated wod timing
                wod_marker = MuseStruct{1}{1}.markers.WOD.synctime(itrial);
                
                %select times where to search WOD peak
                t = LFP_lpfilt.time{1};
                t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
                t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
                t_sel = t_1 & t_2;
                %Search LFP maximum peak in this selected window. '-'LFP because wod is negative
                [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                clear t t_1 t_2 t_sel
                
                %keep only data between 0 and wod
                cfgtemp                                   = [];
                cfgtemp.latency                           = [0 t_peak_wod];
                timefreq_wod{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wod{itrial}.(ichan_name));
                
                %normalize time
                timefreq_wod_timenorm{itrial}.(ichan_name)        = timefreq_wod{itrial}.(ichan_name);
                timefreq_wod_timenorm{itrial}.(ichan_name).time   = timefreq_wod{itrial}.(ichan_name).time ./ t_peak_wod;
                %check the location of the peak detection
                %                 fig= figure;hold
                %                 plot(LFP_trial_lpfilt.time{1}, LFP_trial_lpfilt.trial{1}(ichan,:));
                %                 scatter(t_peak_wod, -v_peak_wod,50,'xr');
                %                 xlim([t_peak_wod-10 t_peak_wod+10]);
                %                 print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'Detection',sprintf('Rat%g_WOD%g_%s.png',irat,itrial,ichan_name)),'-r600');
                %                 close all
                
                for ifreq = 1:size(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm,2)
                    %resample data to have the same number of data points, for
                    %time-normalized data
                    t_old                                                = timefreq_wod_timenorm{itrial}.(ichan_name).time;
                    t_new                                                = linspace(0,1,1000);
                    %powspctrm_new(1,1,:)                                = pchip(t_old,squeeze(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
                    powspctrm_new(1,ifreq,:)                             = interp1(t_old,squeeze(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:)),t_new);
                end
                timefreq_wod_timenorm{itrial}.(ichan_name).time         = t_new;
                timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm    = powspctrm_new;
                clear powspctrm_new
                %plot(timefreq_wod{itrial}.time{1},timefreq_wod{itrial}.trial{1}); xlim([0 0.95])
                
                
                %% baseline short
                
                cfgtemp = [];
                cfgtemp.latency = [-300 0];
                timefreq_baseline{itrial}.(ichan_name) = ft_selectdata(cfgtemp,timefreq_ichan_temp);
                
                
                %% MAKE BASELINE CORRECTION FOR SHORT PERIODS
                
                timefreq_baseline_blcorrected{itrial}.(ichan_name)= timefreq_baseline{itrial}.(ichan_name);
                timefreq_wod_blcorrected{itrial}.(ichan_name)= timefreq_wod{itrial}.(ichan_name);
                timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name) = timefreq_wod_timenorm{itrial}.(ichan_name);
                
                for ifreq = 1:size(timefreq_wod{itrial}.(ichan_name).freq,2) %baseline
                    baseline_ifreq = nanmean(squeeze(timefreq_baseline{itrial}.(ichan_name).powspctrm(1,ifreq,:))); %baseline
                    timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_baseline{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq; %baseline
                    timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wod{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
                    timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
                end %baseline
                
                clear timefreq_ichan_temp
                
                %% Make Logarythm for short periods
                
                %duplicate fieldtrip structure
                log_timefreq_baseline_blcorrected{itrial}.(ichan_name)= timefreq_baseline_blcorrected{itrial}.(ichan_name);
                log_timefreq_baseline{itrial}.(ichan_name)=timefreq_baseline{itrial}.(ichan_name);
                log_timefreq_wod_blcorrected{itrial}.(ichan_name)=timefreq_wod_blcorrected{itrial}.(ichan_name);
                log_timefreq_wod{itrial}.(ichan_name)=timefreq_wod{itrial}.(ichan_name);
                log_timefreq_wod_timenorm{itrial}.(ichan_name)=timefreq_wod_timenorm{itrial}.(ichan_name);
                log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name)=timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name);
                
                %apply logarythm to powspctrm
                log_timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm);
                log_timefreq_baseline{itrial}.(ichan_name).powspctrm= log10(timefreq_baseline{itrial}.(ichan_name).powspctrm);
                log_timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm);
                log_timefreq_wod{itrial}.(ichan_name).powspctrm= log10(timefreq_wod{itrial}.(ichan_name).powspctrm);
                log_timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm);
                log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm);
                
                
                
                %% remove cfg fields as it is what takes the most of place on disk, whereas we do not use it later
                timefreq_wod{itrial}.(ichan_name)            = rmfield(timefreq_wod{itrial}.(ichan_name),{'cfg'});
                timefreq_wod_blcorrected{itrial}.(ichan_name)            = rmfield(timefreq_wod_blcorrected{itrial}.(ichan_name),{'cfg'});
                timefreq_wod_timenorm{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm{itrial}.(ichan_name),{'cfg'});
                timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
                timefreq_wor{itrial}.(ichan_name)            = rmfield(timefreq_wor{itrial}.(ichan_name),{'cfg'});
                timefreq_wor_blcorrected{itrial}.(ichan_name)            = rmfield(timefreq_wor_blcorrected{itrial}.(ichan_name),{'cfg'});
                timefreq_wor_timenorm{itrial}.(ichan_name) 	 = rmfield(timefreq_wor_timenorm{itrial}.(ichan_name),{'cfg'});
                timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name) 	 = rmfield(timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
                timefreq_recovery{itrial}.(ichan_name)                   = rmfield(timefreq_recovery{itrial}.(ichan_name),{'cfg'});
                timefreq_recovery_blcorrected{itrial}.(ichan_name)    	 = rmfield(timefreq_recovery_blcorrected{itrial}.(ichan_name),{'cfg'});
                timefreq_baseline{itrial}.(ichan_name)                  = rmfield(timefreq_baseline{itrial}.(ichan_name),{'cfg'});
                timefreq_baseline_blcorrected{itrial}.(ichan_name)       = rmfield(timefreq_baseline_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_baseline_blcorrected{itrial}.(ichan_name)    = rmfield(log_timefreq_baseline_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_baseline{itrial}.(ichan_name)                = rmfield(log_timefreq_baseline{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wod_blcorrected{itrial}.(ichan_name)         = rmfield(log_timefreq_wod_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wod{itrial}.(ichan_name)                     = rmfield(log_timefreq_wod{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wod_timenorm{itrial}.(ichan_name)            = rmfield(log_timefreq_wod_timenorm{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name)= rmfield(log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_recovery_blcorrected{itrial}.(ichan_name)    = rmfield(log_timefreq_recovery_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_recovery{itrial}.(ichan_name)               = rmfield(log_timefreq_recovery{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wor{itrial}.(ichan_name)                     = rmfield(log_timefreq_wor{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wor_blcorrected{itrial}.(ichan_name)         = rmfield(log_timefreq_wor_blcorrected{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wor_timenorm{itrial}.(ichan_name)            = rmfield(log_timefreq_wor_timenorm{itrial}.(ichan_name),{'cfg'});
                log_timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name)= rmfield(log_timefreq_wor_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
                
                
            end %ichan
        end %itrial
        
        
        % add empty missing channels channels to have the same channels between rats
        for itrial = 1:size(timefreq_recovery,2)
            for chan_name = string(config{irat}.LFP.allchannel(1:end-1))
                chan_renamed = sprintf('E%d',str2num(regexp(chan_name,'\d*','Match')));
                if ~isfield(timefreq_recovery{itrial}, chan_renamed)
                    timefreq_recovery{itrial}.(chan_renamed)            = [];
                    timefreq_recovery_blcorrected{itrial}.(chan_renamed)            = [];
                    
                    timefreq_wod{itrial}.(chan_renamed)                 = [];
                    timefreq_wod_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    timefreq_wod_timenorm{itrial}.(chan_renamed)        = [];
                    timefreq_wod_timenorm_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    timefreq_wor{itrial}.(chan_renamed)                 = [];
                    timefreq_wor_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    timefreq_wor_timenorm{itrial}.(chan_renamed)        = [];
                    timefreq_wor_timenorm_blcorrected{itrial}.(chan_renamed)        = [];
                    
                    timefreq_baseline{itrial}.(chan_renamed)            = [];
                    timefreq_baseline_blcorrected{itrial}.(chan_renamed)            = [];
                    
                    log_timefreq_recovery{itrial}.(chan_renamed)            = [];
                    log_timefreq_recovery_blcorrected{itrial}.(chan_renamed)            = [];
                    
                    log_timefreq_wod{itrial}.(chan_renamed)                 = [];
                    log_timefreq_wod_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    log_timefreq_wod_timenorm{itrial}.(chan_renamed)        = [];
                    log_timefreq_wod_timenorm_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    log_timefreq_wor{itrial}.(chan_renamed)                 = [];
                    log_timefreq_wor_blcorrected{itrial}.(chan_renamed)                 = [];
                    
                    log_timefreq_wor_timenorm{itrial}.(chan_renamed)        = [];
                    log_timefreq_wor_timenorm_blcorrected{itrial}.(chan_renamed)        = [];
                    
                    log_timefreq_baseline{itrial}.(chan_renamed)            = [];
                    log_timefreq_baseline_blcorrected{itrial}.(chan_renamed)            = [];
                    
                end
            end
        end
        
        %save time freq data to disk :
        
        %WOD data (entre vent_off et pic de wod), power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod.mat']),'timefreq_wod','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_blcorrected.mat']),'timefreq_wod_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wod.mat']),'log_timefreq_wod','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wod_blcorrected.mat']),'log_timefreq_wod_blcorrected','-v7.3');
        
        %WOD data : time normalized between vent_off and wod peak, per channel. power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_timenorm.mat']),'timefreq_wod_timenorm','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_timenorm_blcorrected.mat']),'timefreq_wod_timenorm_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wod_timenorm.mat']),'log_timefreq_wod_timenorm','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wod_timenorm_blcorrected.mat']),'log_timefreq_wod_timenorm_blcorrected','-v7.3');
        
        %Recovery data : t0 at Vent_On, t_end at Vent_On+3600s:
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_recovery.mat']),'timefreq_recovery','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_recovery_blcorrected.mat']),'timefreq_recovery_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_recovery.mat']),'log_timefreq_recovery','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_recovery_blcorrected.mat']),'log_timefreq_recovery_blcorrected','-v7.3');
        
        %WOR data : power normalized with baseline period (time not normalized) :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor.mat']),'timefreq_wor','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor_blcorrected.mat']),'timefreq_wor_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wor.mat']),'log_timefreq_wor','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wor_blcorrected.mat']),'log_timefreq_wor_blcorrected','-v7.3');
        
        %WOR data : time normalized between vent_off and wor peak, per channel. power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor_timenorm.mat']),'timefreq_wor_timenorm','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor_timenorm_blcorrected.mat']),'timefreq_wor_timenorm_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wor_timenorm.mat']),'log_timefreq_wor_timenorm','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_wor_timenorm_blcorrected.mat']),'log_timefreq_wor_timenorm_blcorrected','-v7.3');
        
        %Baseline data: t0 at Vent_Off analysis made before Vent_Off
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_baseline.mat']),'timefreq_baseline','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_baseline_blcorrected.mat']),'timefreq_baseline_blcorrected','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_baseline.mat']),'log_timefreq_baseline','-v7.3');
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'log_timefreq_baseline_blcorrected.mat']),'log_timefreq_baseline_blcorrected','-v7.3');
        
        
        
    end %irat
    
end % wod_project
