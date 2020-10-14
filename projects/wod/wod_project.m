function wod_project(slurm_task_id)

% time frequency analysis of WOD Neuralynx data
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute data for each rat %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ipart = 1; %ipart is always 1 for this project

if slurm_task_id(1) > 0
    for irat = slurm_task_id
        
        %find concatenated LFP (see wod_concatenateLFP.m)
        [~,dir_name]                       = fileparts(config{irat}.rawdir);
        config{irat}.rawdir                = fullfile(config{irat}.datasavedir,'concatenated_LFP');
        config{irat}.directorylist{ipart}  = {dir_name};
        
        %read Muse markers
        MuseStruct               = readMuseMarkers(config{irat}, true);
        
        %read LFP. T0 = Vent_Off. Each trial is one protocol
        LFP = readLFP(config{irat}, MuseStruct, true);
        LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
        
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
        
        %remove breathing channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0', '-Respi'};
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
            cfgtemp.trials      = itrial;
            cfgtemp.lpfilter    = 'yes';
            cfgtemp.lpfilttype  = 'fir';
            cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
            LFP_trial_filt      = ft_preprocessing(cfgtemp, LFP);
            plot(LFP_trial_filt.time{1}(1:640),LFP_trial_filt.trial{1}(1,1:640))
            plot(LFP.time{1}(1:640),LFP.trial{1}(1,1:640))
            
            %recover trial real timings to use it with muse markers
            starttrial              = LFP.trialinfo.begsample(itrial) / LFP.fsample;
            endtrial                = LFP.trialinfo.endsample(itrial) / LFP.fsample;
            offsettrial             = LFP.trialinfo.offset(itrial) / LFP.fsample;
            
            
            %do time frequency analysis
            cfgtemp                         = [];
            cfgtemp.channel                 = 'all';
            cfgtemp.method                  = 'mtmconvol';
            cfgtemp.output                  = 'pow';
            cfgtemp.taper                   = 'hanning';
            cfgtemp.pad                     = 'nextpow2';
            cfgtemp.keeptrials              = 'no';
            cfgtemp.foi                     = config{irat}.timefreq.foi;
            cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*config{irat}.timefreq.t_ftimwin;
            cfgtemp.toi                     = LFP_trial.time{1}(1) : config{irat}.timefreq.timestep : LFP_trial.time{1}(end);
            timefreq_alldata{itrial} = ft_freqanalysis(cfgtemp,LFP_trial);
            
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
                    t_tfr                   = timefreq_alldata{itrial}.time;
                    t_lfp                   = LFP.time{itrial};
                    bad_sel                 = find(bad_start >= starttrial & bad_start <= endtrial);
                    %go through each bad timing
                    for ibad = bad_sel
                        %remove lfp artefacts
                        bad_period_lfp = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) & t_lfp <= (bad_end(ibad)- starttrial + offsettrial);
                        LFP_cleaned.trial{itrial}(:,bad_period_lfp) = nan(size(LFP_cleaned.trial{itrial},1),sum(bad_period_lfp));
                        %remove tfr artefacts: all window with at least one artefacted sample
                        bad_period_tfr = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - config{irat}.timefreq.t_ftimwin/2 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + config{irat}.timefreq.t_ftimwin/2;
                        timefreq_alldata{itrial}.powspctrm(:,:,bad_period_tfr) = nan(size(timefreq_alldata{itrial}.powspctrm,1),size(timefreq_alldata{itrial}.powspctrm,2),sum(bad_period_tfr));
                        %check artefact removal
                        %                             figure;hold;
                        %                             sel = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) - 10 & t_lfp <= (bad_end(ibad)- starttrial + offsettrial) + 10;
                        %                             plot(LFP.time{itrial}(sel),LFP.trial{itrial}(1,sel),'r');
                        %                             plot(LFP_cleaned.time{itrial}(sel),LFP_cleaned.trial{itrial}(1,sel),'k');
                        %                             sel = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - 10 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + 10;
                        %                             plot(timefreq_alldata{itrial}.time(sel),squeeze(timefreq_alldata{itrial}.powspctrm(1,1,sel)),'k');
                    end
                end
            end
            
            %check TFR results, first channel on top
            %                              figure;hold;
            %                             h = 50000;
            %                             for ichan = 1:size(timefreq_alldata{itrial}.label,1)
            %                                  plot(timefreq_alldata{itrial}.time,squeeze(timefreq_alldata{itrial}.powspctrm(ichan,:,:))+h*(size(timefreq_alldata{itrial}.label,1)-ichan));
            %                              end
            %                              ylim([0 (ichan+1)*h]);
            %
            %normalize time freq with the baseline
            cfgtemp                         = [];
            cfgtemp.baseline                = [config{irat}.epoch.toi.(config{irat}.LFP.name{1})(1) 0];
            cfgtemp.baselinetype            = 'relchange'; % /mean. relchange : -mean / mean
            timefreq_alldata{itrial} = ft_freqbaseline(cfgtemp, timefreq_alldata{itrial});
            
            %check normalization
            %                             figure;hold;
            %                             h = 10;
            %                             for ichan = 1:size(timefreq_alldata{itrial}.label,1)
            %                                 plot(timefreq_alldata{itrial}.time,squeeze(timefreq_alldata{itrial}.powspctrm(ichan,:,:))+h*(size(timefreq_alldata{itrial}.label,1)-ichan));
            %                             end
            %                             ylim([0 (ichan+1)*h]);
            %
            %             average freq values for the whole frequency band
            %             cfgtemp                     = [];
            %             cfgtemp.frequency           = 'all';
            %             cfgtemp.avgoverfreq         = 'yes';
            %             cfgtemp.nanmean             = 'yes';
            %             timefreq_alldata{itrial} = ft_selectdata(cfgtemp, timefreq_alldata{itrial});
            %                             plot(timefreq_alldata{itrial}.time, squeeze(timefreq_alldata{itrial}.powspctrm(:,1,:))); ylim([0 10]); xlim([-900 0])
            %
            %go trough each channel to output results, as there will not
            %be the same amout of data between each channel (WOD/WOR do not
            %occur at the same time) :
            for ichan = 1:size(timefreq_alldata{itrial}.label,1)
                
                %select channel
                ichan_name              = timefreq_alldata{itrial}.label{ichan};
                cfgtemp                 = [];
                cfgtemp.channel         = ichan_name;
                timefreq_ichan_temp   	= ft_selectdata(cfgtemp,timefreq_alldata{itrial});
                
                %% RECOVERY DATA : Change T0 from Vent_Off to Vent_On
                %MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_On.synctime(itrial);
                timeshift                                           = config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2) + config{irat}.epoch.pad.(config{irat}.LFP.name{1}) - timefreq_ichan_temp.time(end);%MuseStruct{1}{1}.markers.Vent_On.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial);
                timefreq_recovery{itrial}.(ichan_name)       = timefreq_ichan_temp;
                timefreq_recovery{itrial}.(ichan_name).time 	= timefreq_ichan_temp.time + timeshift;
                
                cfgtemp = [];
                cfgtemp.latency = [0 config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2)];
                timefreq_recovery{itrial}.(ichan_name) = ft_selectdata(cfgtemp,timefreq_recovery{itrial}.(ichan_name));
                
                %% WOD DATA : find WOD peak per channel, and normalize time per channel
                %use filtered data to find wod
                
                timefreq_wod{itrial}.(ichan_name)            = timefreq_ichan_temp;
                
                %select lfp channel with the same name as tfr channel (in
                %case channel numbers were schuffled by fieldtrip)
                chan_idx    = strcmp(LFP_trial_filt.label, ichan_name);
                
                %get hand-annotated wod timing
                wod_marker = MuseStruct{1}{1}.markers.WOD.synctime(itrial);
                
                %select times where to search WOD peak
                t = LFP_trial_filt.time{1};
                t_1 = t > (wod_marker + config{irat}.LFP.wod_toisearch(1) - starttrial + offsettrial);
                t_2 = t < (wod_marker + config{irat}.LFP.wod_toisearch(2) - starttrial + offsettrial);
                t_sel = t_1 & t_2;
                %Search LFP maximum peak in this selected window. '-'LFP because wod is negative
                [v_peak_wod, t_peak_wod] = findpeaks(-LFP_trial_filt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                
                %keep only data between 0 and wod
                cfgtemp                                   = [];
                cfgtemp.latency                           = [0 t_peak_wod];
                timefreq_wod{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wod{itrial}.(ichan_name));
                
                %normalize time
                timefreq_wod_timenorm{itrial}.(ichan_name)        = timefreq_wod{itrial}.(ichan_name);
                timefreq_wod_timenorm{itrial}.(ichan_name).time   = timefreq_wod{itrial}.(ichan_name).time ./ t_peak_wod;
                %check the location of the peak detection
                fig= figure;hold
                plot(LFP_trial_filt.time{1}, LFP_trial_filt.trial{1}(ichan,:));
                scatter(t_peak_wod, -v_peak_wod,50,'xr');
                xlim([t_peak_wod-10 t_peak_wod+10]);
                print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'Detection',sprintf('Rat%g_WOD%g_%s.png',irat,itrial,ichan_name)),'-r600');
                
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
                
                
                %% WOR DATA : find WOR peak per channel, and normalize time per channel
                %Antoine
                %+ normaliser le temps que selon la 1ere wod/wor, pour toutes les électrodes du rat ? pour
                %mieux voir le délai entre électrodes
                
                timefreq_wor{itrial}.(ichan_name)            = timefreq_ichan_temp;
                
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
                %Search LFP maximum peak in this selected window. wor is positive
                [v_peak_wor, t_peak_wor] = findpeaks(LFP_trial_filt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                end_wor= timefreq_wor{itrial}.(ichan_name).time(end);
                %keep only data between wod and end
                cfgtemp                                   = [];
                cfgtemp.latency                           = [t_peak_wor end_wor];
                timefreq_wor{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wor{itrial}.(ichan_name));
                
                %normalize time
                timefreq_wor_timenorm{itrial}.(ichan_name)        = timefreq_wor{itrial}.(ichan_name);
                timefreq_wor_timenorm{itrial}.(ichan_name).time   = timefreq_wor{itrial}.(ichan_name).time-t_peak_wor;%./end_wor;
                timefreq_wor_timenorm{itrial}.(ichan_name).time   = timefreq_wor_timenorm{itrial}.(ichan_name).time./timefreq_wor_timenorm{itrial}.(ichan_name).time(end);%./end_wor;
                %check the location of the peak detection
                fig=figure;hold
                plot(LFP_trial_filt.time{1}, LFP_trial_filt.trial{1}(ichan,:));
                scatter(t_peak_wor, v_peak_wor,50,'red');
                xlim([t_peak_wor-10 t_peak_wor+10]);
                print(fig, '-dpng', fullfile(config{irat}.imagesavedir,'Detection',sprintf('Rat%g_WOR%g_%s.png',irat,itrial,ichan_name)),'-r600');
                
                for ifreq = 1:size(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm,2)
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
                %                     plot(timefreq_wor_timenorm{itrial}.(ichan_name).time,squeeze(timefreq_wor_timenorm{itrial}.(ichan_name).powspctrm),'r'); %xlim([0 0.95])
                
                
                %%%%%%%%%%%%%%%%do general time frequency map
                %%%%%%%%%%%%%%%%from Vent_Off to peak WoD
                %                 cfgteam                         = [];
                %                 cfgtemp.channel                 = 'all';
                %                 cfgtemp.method                  = 'mtmconvol';
                %                 cfgtemp.output                  = 'pow';
                %                 cfgtemp.taper                   = 'hanning';
                %                 cfgtemp.pad                     = 'nextpow2';
                %                 cfgtemp.keeptrials              = 'no';
                %                 cfgtemp.foi                     = 0.1 : 1 : 50.1;
                %                 cfgtemp.toi                     =
                
                
                
                %% remove cfg fields as it is what takes the most of place on disk, whereas we do not use it later
                timefreq_wod{itrial}.(ichan_name)             = rmfield(timefreq_wod{itrial}.(ichan_name),{'cfg'});
                timefreq_wod_timenorm{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm{itrial}.(ichan_name),{'cfg'});
                timefreq_wor{itrial}.(ichan_name)             = rmfield(timefreq_wor{itrial}.(ichan_name),{'cfg'});
                timefreq_wor_timenorm{itrial}.(ichan_name) 	 = rmfield(timefreq_wor_timenorm{itrial}.(ichan_name),{'cfg'});
                timefreq_recovery{itrial}.(ichan_name)      	 = rmfield(timefreq_recovery{itrial}.(ichan_name),{'cfg'});
                
            end %ichan
        end
        
        
        % add empty missing channels channels to have the same channels between rats
        for itrial = 1:size(timefreq_recovery,2)
            for chan_name = string(config{irat}.LFP.allchannel(1:end-1))
                chan_renamed = sprintf('E%d',str2num(regexp(chan_name,'\d*','Match')));
                if ~isfield(timefreq_recovery{itrial}, chan_renamed)
                    timefreq_recovery{itrial}.(chan_renamed)            = [];
                    timefreq_wod{itrial}.(chan_renamed)                 = [];
                    timefreq_wod_timenorm{itrial}.(chan_renamed)        = [];
                    timefreq_wor{itrial}.(chan_renamed)                 = [];
                    timefreq_wor_timenorm{itrial}.(chan_renamed)        = [];
                end
            end
        end
        
        %save time freq data to disk :
        
        %WOD data (entre vent_off et pic de wod), power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod.mat']),'timefreq_wod','-v7.3');
        %WOD data : time normalized between vent_off and wod peak, per channel. power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_timenorm.mat']),'timefreq_wod_timenorm','-v7.3');
        %Recovery data : t0 at Vent_On, t_end at Vent_On+3600s:
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_recovery.mat']),'timefreq_recovery','-v7.3');
        %WOR data : power normalized with baseline period (time not normalized) :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor.mat']),'timefreq_wor','-v7.3');
        %WOR data : time normalized between vent_off and wor peak, per channel. power normalized with baseline period :
        save(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor_timenorm.mat']),'timefreq_wor_timenorm','-v7.3');
        
    end %irat
    
end %slurm_task_id >0



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% average over rats, plot, and compute stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if slurm_task_id(1) == 0
    
    cfgcommon = config{4}; %same common config for all rats
    
    %initialize data structures
    timefreq_recovery_allrats       = [];
    timefreq_wod_allrats            = [];
    timefreq_wod_timenorm_allrats 	= [];%cell(1,size(cfgcommon.timefreq.foi,2));
    timefreq_wor_allrats            = [];%cell(1,size(cfgcommon.timefreq.foi,2));
    timefreq_wor_timenorm_allrats 	= [];%cell(1,size(cfgcommon.timefreq.foi,2));
    count_trials                    = 0;
    
    % load timefreq data from all rats
    for irat = 4:6%[4 5 6 7 8 9 10 11] %4:size(config,2)
        
        %load data
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_recovery.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_timenorm.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wor_timenorm.mat']));
        
        %pool all rat data in the same structure
        for itrial = 1:size(timefreq_recovery{1},2)
            count_trials = count_trials + 1;
            wod_rat(count_trials) = irat; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(timefreq_wod{itrial});
            for ichan = 1:numel(chan_list)
                chan_name                                                       = chan_list{ichan};
                timefreq_recovery_allrats.(chan_name){count_trials}     	= timefreq_recovery{itrial}.(chan_name);
                timefreq_wod_allrats.(chan_name){count_trials}           = timefreq_wod{itrial}.(chan_name);
                timefreq_wod_timenorm_allrats.(chan_name){count_trials} 	= timefreq_wod_timenorm{itrial}.(chan_name);
                timefreq_wor_allrats.(chan_name){count_trials}           = timefreq_wor{itrial}.(chan_name);
                timefreq_wor_timenorm_allrats.(chan_name){count_trials} 	= timefreq_wor_timenorm{itrial}.(chan_name);
            end
        end
        
        clear timefreq_recovery timefreq_wod timefreq_wod_timenorm timefreq_wor timefreq_wor_timenorm
    end
    
    % pool data for each frequency band and channel
    % if not the same number of points, keep only the number of points of
    % the shorter trial
    chan_list                                   = fieldnames(timefreq_wod_timenorm_allrats);
    for ichan = 1:numel(chan_list)
        chan_name                               = chan_list{ichan};
        
        %replace empty channels by nan channels
        hasdata = find(~cellfun(@isempty,timefreq_wod_timenorm_allrats.(chan_name)));
        for iwod = 1:size(timefreq_wod_allrats.(chan_name),2)
            if isempty(timefreq_wod_allrats.(chan_name){iwod})
                timefreq_wod_allrats.(chan_name){iwod}                       = timefreq_wod_allrats.(chan_name){hasdata(1)};
                timefreq_wod_allrats.(chan_name){iwod}.powspctrm             = nan(size(timefreq_wod_allrats.(chan_name){iwod}.powspctrm));
                timefreq_wor_allrats.(chan_name){iwod}                       = timefreq_wor_allrats.(chan_name){hasdata(1)};
                timefreq_wor_allrats.(chan_name){iwod}.powspctrm             = nan(size(timefreq_wor_allrats.(chan_name){iwod}.powspctrm));
                timefreq_recovery_allrats.(chan_name){iwod}                  = timefreq_recovery_allrats.(chan_name){hasdata(1)};
                timefreq_recovery_allrats.(chan_name){iwod}.powspctrm        = nan(size(timefreq_recovery_allrats.(chan_name){iwod}.powspctrm));
                timefreq_wod_timenorm_allrats.(chan_name){iwod}              = timefreq_wod_timenorm_allrats.(chan_name){hasdata(1)};
                timefreq_wod_timenorm_allrats.(chan_name){iwod}.powspctrm    = nan(size(timefreq_wod_timenorm_allrats.(chan_name){iwod}.powspctrm));
                timefreq_wor_timenorm_allrats.(chan_name){iwod}            	 = timefreq_wor_timenorm_allrats.(chan_name){hasdata(1)};
                timefreq_wor_timenorm_allrats.(chan_name){iwod}.powspctrm    = nan(size(timefreq_wor_timenorm_allrats.(chan_name){iwod}.powspctrm));
            end
        end
        
        %data separated for each channel
        cfgtemp                             = [];
        cfgtemp.keepindividual              = 'yes';
        timefreq_wod_grandavg.(chan_name)               = ft_freqgrandaverage(cfgtemp, timefreq_wod_allrats.(chan_name){:});
        timefreq_wod_timenorm_grandavg.(chan_name)     	= ft_freqgrandaverage(cfgtemp, timefreq_wod_timenorm_allrats.(chan_name){:});
        timefreq_wor_grandavg.(chan_name)               = ft_freqgrandaverage(cfgtemp, timefreq_wor_allrats.(chan_name){:});
        timefreq_wor_timenorm_grandavg.(chan_name)     	= ft_freqgrandaverage(cfgtemp, timefreq_wor_timenorm_allrats.(chan_name){:});
        timefreq_recovery_grandavg.(chan_name)          = ft_freqgrandaverage(cfgtemp, timefreq_recovery_allrats.(chan_name){:});
  
    end
    
    clear timefreq_recovery_allrats timefreq_wod_allrats timefreq_wod_timenorm_allrats
    
    
    if ~isfolder(cfgcommon.imagesavedir)
        fprintf('Creating directory %s\n', cfgcommon.imagesavedir);
        mkdir(cfgcommon.imagesavedir)
    end
    
    analysis_names = {'Timefreq_WOD', 'Timefreq_WOD_timenorm', 'Timefreq_recovery','Timefreq_WOR','Timefreq_WOR_timenorm'};
    data           = {timefreq_wod_grandavg, timefreq_wod_timenorm_grandavg, timefreq_recovery_grandavg, timefreq_wor_grandavg, timefreq_wor_timenorm_grandavg};
    
    %go through each data structure
    for idata = 1:size(data,2)
        
        %% TFR average over chans and over rats
        
        %gather powspctrm data for each channel
        powspctrm_allchans = nan(size(data{idata}.(chan_name).powspctrm));
        ichan=0;
        for chan_name = string(fieldnames(data{idata}))'
            ichan=ichan+1;
            
            if size(data{idata}.(chan_name).powspctrm,4) > size(powspctrm_allchans,4)
                data{idata}.(chan_name).powspctrm = data{idata}.(chan_name).powspctrm(:,:,:,1:size(powspctrm_allchans,4));
            elseif size(data{idata}.(chan_name).powspctrm,4) < size(powspctrm_allchans,4)
                powspctrm_allchans = powspctrm_allchans(:,:,:,1:size(powspctrm_allchans,4));
            end
            powspctrm_allchans(:,ichan,:,:) = data{idata}.(chan_name).powspctrm;
        end

        %re-create fieldtrip structure
        datafreq_allchans.label = fieldnames(data{idata})';
        datafreq_allchans.powspctrm = powspctrm_allchans;
        datafreq_allchans.freq = data{idata}.(chan_name).freq;
        datafreq_allchans.time = data{idata}.(chan_name).time(1:size(powspctrm_allchans,4));
        datafreq_allchans.dimord = 'subj_chan_freq_time';
        
        %average over protocols
        data_plot = datafreq_allchans;
        data_plot.label = {'average channels'};
        data_plot.powspctrm = squeeze(nanmean(datafreq_allchans.powspctrm,1));
        %data_plot.powspctrm = squeeze((datafreq_allchans.powspctrm(iwod,:,:,:)));
        data_plot.powspctrm = nanmean(data_plot.powspctrm,1);
        data_plot.dimord = 'chan_freq_time';
        
        fig=figure;
        
        sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        subplot(2,1,1)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfg.fontsize = 12;
        cfgtemp.ylim= [0 49]
        cfgtemp.masknans    = 'yes';
        ft_singleplotTFR(cfgtemp, data_plot);
        ylim([0 49]);
        caxis([-2 2]);
        

           
        
        subplot(2,1,2)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfgtemp.masknans    = 'yes';
        cfgtemp.ylim= [51 100]
        cfg.fontsize = 12;
        ft_singleplotTFR(cfgtemp, data_plot);
        ylim([51 100]);
        caxis([-2 2]);
        
   %Title and label plot then save plot with name

            
            
        


%figure settings (fontsize,fontweight, ticks dir, renderer etc.)
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            ylabel(sprintf('Frequency (Hz)'));
             axis tight;
            if contains(analysis_names{idata}, 'WOD')
                ylim([0 max(data_max)]);
            end 
% print image to file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,['AllRats_TFR_',analysis_names{idata},'.pdf']),'-r600');
print(fig, '-dpng', fullfile(cfgcommon.imagesavedir,['AllRats_TFR_',analysis_names{idata},'.png']),'-r600');
savefig(fig,fullfile(cfgcommon.imagesavedir,['AllRats_TFR_',analysis_names{idata}]));
              
        
        print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,sprintf('%s_AllRats_TFR.pdf',analysis_names{idata})),'-r600');
             
        
        %% TFR average over chans by rat
%get wod_rat and wod_rat_number for titles
%gather powspctrm data for each channel
        powspctrm_allchans = nan(size(data{idata}.(chan_name).powspctrm));
        ichan=0;
        for chan_name = string(fieldnames(data{idata}))'
            ichan=ichan+1;
            
            if size(data{idata}.(chan_name).powspctrm,4) > size(powspctrm_allchans,4)
                data{idata}.(chan_name).powspctrm = data{idata}.(chan_name).powspctrm(:,:,:,1:size(powspctrm_allchans,4));
            elseif size(data{idata}.(chan_name).powspctrm,4) < size(powspctrm_allchans,4)
                powspctrm_allchans = powspctrm_allchans(:,:,:,1:size(powspctrm_allchans,4));
            end
            powspctrm_allchans(:,ichan,:,:) = data{idata}.(chan_name).powspctrm;
        end

        %re-create fieldtrip structure
        datafreq_allchans.label = fieldnames(data{idata})';
        datafreq_allchans.powspctrm = powspctrm_allchans;
        datafreq_allchans.freq = data{idata}.(chan_name).freq;
        datafreq_allchans.time = data{idata}.(chan_name).time(1:size(powspctrm_allchans,4));
        datafreq_allchans.dimord = 'subj_chan_freq_time';
        
        %average over protocols
        data_plot = datafreq_allchans;
        data_plot.label = {'average channels'};
        data_plot.powspctrm = squeeze(nanmean(datafreq_allchans.powspctrm,1));
        %data_plot.powspctrm = squeeze((datafreq_allchans.powspctrm(iwod,:,:,:)));
        data_plot.powspctrm = nanmean(data_plot.powspctrm,1);
        data_plot.dimord = 'chan_freq_time';
        
        fig=figure;
        
        sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        subplot(2,1,1)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfg.fontsize = 12;
        cfgtemp.ylim= [0 49]
        cfgtemp.masknans    = 'yes';
        ft_singleplotTFR(cfgtemp, data_plot);
        caxis([-2 2]);
        

           
        
        subplot(2,1,2)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfg.fontsize = 12;
        cfgtemp.ylim= [51 100]
        cfgtemp.masknans    = 'yes';
        ft_singleplotTFR(cfgtemp, data_plot);
        caxis([-2 2]);
        
   %Title and label plot then save plot with name

            
            
        


%figure settings (fontsize,fontweight, ticks dir, renderer etc.)
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            ylabel(sprintf('Frequency (Hz)'));
             axis tight;
        

%         
        %% plot frequency bands over time for all rats
        fig = figure;
        sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
            
            subplot(2,2,ifreq);hold;
            C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
            
            for chan_name = string(fieldnames(data{idata}))'
                %go trhough each electrode in the ascending order
                chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
                legend_name{chan_nr} = chan_name;
                %                 if isempty(data{idata}{ifreq}.(chan_name))
                %                     continue
                %                 end
                color   = C(chan_nr,:);
                
                x = data{idata}.(chan_name).time;
                
                %find frequencies of the band (comparer des valeurs : > <
                %>= <= ==)
                freq_band_idx = data{idata}.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data{idata}.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
                pow_band= data{idata}.(chan_name).powspctrm(:,:,freq_band_idx,:);
                %average over frequencies and over protocols
                pow_band = nanmean(pow_band,1); %average over protocols
                pow_band = nanmean(pow_band,3); %average over frequencies
                y = squeeze(pow_band);
                
                if ismember(idata, [3 4 5]) %only for recovery, wor and wor timenorm
                    y = movmean(y,100,'omitnan');
                end
                
                %plot std
                %                 std = nanstd(squeeze(data{idata}{ifreq}.(chan_name).powspctrm(:,1,1,:)));
                %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
                %                 filled_SD = area(x,y_area);
                %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
                %                 filled_SD(1).ShowBaseLine = 'off';
                
                %find maximum of the first half of data to adapt y scale
                data_max(chan_nr) = max(y(x>0 & x<x(end)/2)+std(y(x>0 & x<x(end)/2)));
                
                %plot avg
                leg{chan_nr} = plot(x, y, 'Color', color);
            end
            
            %set figure display :
            %set(gca, 'YScale', 'log');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(end)));
            leg = flip(leg); legend_name = flip(legend_name); %set E16 on top
            legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
            legend('boxoff');
            axis tight;
            if contains(analysis_names{idata}, 'WOD')
                ylim([0 max(data_max)]);
            end
            
        end
        
        %save figure :
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,['AllRats_',analysis_names{idata},'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfgcommon.imagesavedir,['AllRats_',analysis_names{idata},'.png']),'-r600');
        savefig(fig,fullfile(cfgcommon.imagesavedir,['AllRats_',analysis_names{idata}]));
        
        
        
        %% plot data by ferquency band by protocol
        for iwod = 1 :size(data{idata}.(chan_name).powspctrm,1)
            
            fig = figure;
            rat_wod_sum = sum(wod_rat == wod_rat(iwod));
            rat_wod_nr = 1;
            if iwod >1
                if wod_rat(iwod-1) == wod_rat(iwod)
                    rat_wod_nr = 2;
                end
            end
            sgtitle(sprintf('%s Rat %d (%d/%d)',analysis_names{idata},wod_rat(iwod), rat_wod_nr, rat_wod_sum), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
            
            for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
                
                subplot(2,2,ifreq);hold;
                C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
                
                for chan_name = string(fieldnames(data{idata}))'
                    %go trhough each electrode in the ascending order
                    chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
                    legend_name{chan_nr} = chan_name;
                    color   = C(chan_nr,:);
                    
                    x = data{idata}.(chan_name).time;
                    %select wod and frequency band
                    freq_band_idx = data{idata}.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data{idata}.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
                    pow_band= data{idata}.(chan_name).powspctrm(iwod,:,freq_band_idx,:);
                    %average over frequencies and over protocols
                    pow_band = nanmean(pow_band,3); %average over frequencies
                    y = squeeze(pow_band);
                    
                    if ismember(idata, [3 4 5]) %only for recovery
                        y = movmean(y,100,'omitnan');
                    end
                    
                    %plot std
                    %                 std = nanstd(squeeze(data{idata}{ifreq}.(chan_name).powspctrm(:,1,1,:)));
                    %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
                    %                 filled_SD = area(x,y_area);
                    %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                    %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                    %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
                    %                 filled_SD(1).ShowBaseLine = 'off';
                    
                    %find maximum of the first half of data to adapt y scale
                    data_max(chan_nr) = max(y(x>0 & x<x(end)/2)+std(y(x>0 & x<x(end)/2)));
                    
                    %plot avg
                    leg{chan_nr} = plot(x, y, 'Color', color);
                end
                
                %set figure display :
                %set(gca, 'YScale', 'log');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time');
                ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(2)));
                leg = flip(leg); legend_name = flip(legend_name); %set E16 on top
                legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
                legend('boxoff');
                axis tight;
                if contains(analysis_names{idata}, 'WOD')
                    ylim([0 max(data_max)]);
                end
                
            end
            
            %save figure :
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,sprintf('%s_Rat%d_%d.pdf',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
            print(fig, '-dpng', fullfile(cfgcommon.imagesavedir,sprintf('%s_Rat%d_%d.pdf',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
            savefig(fig,fullfile(cfgcommon.imagesavedir,sprintf('%s_WOD%g',analysis_names{idata},iwod)));
            
            close all
        end
        
        % compute statistiques : y a-t-il des profondeurs qui récupèrent plus vite ? Lesquelles ?
        
        % trouver période stable de recovery
        
    end %idata
    
end %slurm task id == 0


