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
        
        %rename chans according to their real deepness. Create nan channel if
        %not channel at this deepness. 16 is surface, 1 is the deepest. 0 is the respi.
        n_chans = size(config{irat}.LFP.channel,2);
        for ichan = 1:n_chans
            if isempty(config{irat}.LFP.channel{ichan})
                continue
            else
                %rename the channel to have a standard name between each rat, to allow averaging
                new_name = ['E', num2str(n_chans-ichan)];
                %find idx of this chan :
                chan_idx = strcmp(LFP.label, config{irat}.LFP.channel{ichan});
                LFP.label{chan_idx} = new_name;
            end
        end
        
        %remove breathing channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0'};
        LFP             = ft_selectdata(cfgtemp, LFP);
        LFP_cleaned     = LFP; %save for later removing of artefacts
        
        %analyse each trial independently
        for itrial = 1:size(LFP.trial,2)
            
            %select one trial
            cfgtemp         = [];
            cfgtemp.trials  = itrial;
            LFP_trial       = ft_selectdata(cfgtemp, LFP);
            
            %filter lfp to better recognise WOD/WOR peak
            cfgtemp             = [];
            cfgtemp.trials      = itrial;
            cfgtemp.lpfilter    = 'yes';
            cfgtemp.lpfilttype  = 'fir';
            cfgtemp.lpfreq      = config{irat}.LFP.lpfilter_wod_detection;
            LFP_trial_filt      = ft_preprocessing(cfgtemp, LFP);
            %plot(LFP_trial_filt.time{1}(1:640),LFP_trial_filt.trial{1}(1,1:640))
            %plot(LFP.time{1}(1:640),LFP.trial{1}(1,1:640))
            
            %recover trial real timings to use it with muse markers
            starttrial              = LFP.trialinfo(itrial,5) / LFP.fsample;
            endtrial                = LFP.trialinfo(itrial,6) / LFP.fsample;
            offsettrial             = LFP.trialinfo(itrial,7) / LFP.fsample;
            
            %go trhough each frequency band
            for ifreq = 1:size(config{irat}.timefreq.foi,2)
                
                %do time frequency analysis
                cfgtemp                         = [];
                cfgtemp.channel                 = 'all';
                cfgtemp.method                  = 'mtmconvol';
                cfgtemp.output                  = 'pow';
                cfgtemp.taper                   = 'hanning';
                cfgtemp.pad                     = 'nextpow2';
                cfgtemp.keeptrials              = 'no';
                cfgtemp.foi                     = config{irat}.timefreq.foi{ifreq};
                cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*config{irat}.timefreq.t_ftimwin;
                cfgtemp.toi                     = LFP_trial.time{1}(1) : config{irat}.timefreq.timestep : LFP_trial.time{1}(end);
                timefreq_alldata{ifreq}{itrial} = ft_freqanalysis(cfgtemp,LFP_trial);
                
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
                        t_tfr                   = timefreq_alldata{ifreq}{itrial}.time;
                        t_lfp                   = LFP.time{itrial};
                        bad_sel                 = find(bad_start >= starttrial & bad_start <= endtrial);
                        %go through each bad timing
                        for ibad = bad_sel
                            %remove lfp artefacts
                            bad_period_lfp = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) & t_lfp <= (bad_end(ibad)- starttrial + offsettrial);
                            LFP_cleaned.trial{itrial}(:,bad_period_lfp) = nan(size(LFP_cleaned.trial{itrial},1),sum(bad_period_lfp));
                            %remove tfr artefacts: all window with at least one artefacted sample
                            bad_period_tfr = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - config{irat}.timefreq.t_ftimwin/2 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + config{irat}.timefreq.t_ftimwin/2;
                            timefreq_alldata{ifreq}{itrial}.powspctrm(:,:,bad_period_tfr) = nan(size(timefreq_alldata{ifreq}{itrial}.powspctrm,1),size(timefreq_alldata{ifreq}{itrial}.powspctrm,2),sum(bad_period_tfr));
                            %check artefact removal
                            %figure;hold;
                            %sel = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) - 10 & t_lfp <= (bad_end(ibad)- starttrial + offsettrial) + 10;
                            %plot(LFP.time{itrial}(sel),LFP.trial{itrial}(1,sel),'r');
                            %plot(LFP_cleaned.time{itrial}(sel),LFP_cleaned.trial{itrial}(1,sel),'k');
                            %sel = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - 10 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + 10;
                            %plot(timefreq_alldata{ifreq}{itrial}.time(sel),squeeze(timefreq_alldata{ifreq}{itrial}.powspctrm(1,1,sel)),'k');
                        end
                    end
                end
                
                %check TFR results, first channel on top
%                             figure;hold;
%                             h = 50000;
%                             for ichan = 1:size(timefreq_alldata{ifreq}{itrial}.label,1)
%                                 plot(timefreq_alldata{ifreq}{itrial}.time,squeeze(timefreq_alldata{ifreq}{itrial}.powspctrm(ichan,:,:))+h*(size(timefreq_alldata{ifreq}{itrial}.label,1)-ichan));
%                             end
%                             ylim([0 (ichan+1)*h]);
                
                %normalize time freq with the baseline
                cfgtemp                         = [];
                cfgtemp.baseline                = [config{irat}.epoch.toi.(config{irat}.LFP.name{1})(1) 0];
                cfgtemp.baselinetype            = 'relative'; % /mean. relchange : -mean / mean
                timefreq_alldata{ifreq}{itrial} = ft_freqbaseline(cfgtemp, timefreq_alldata{ifreq}{itrial});
                
                %check normalization
%                             figure;hold;
%                             h = 10;
%                             for ichan = 1:size(timefreq_alldata{ifreq}{itrial}.label,1)
%                                 plot(timefreq_alldata{ifreq}{itrial}.time,squeeze(timefreq_alldata{ifreq}{itrial}.powspctrm(ichan,:,:))+h*(size(timefreq_alldata{ifreq}{itrial}.label,1)-ichan));
%                             end
%                             ylim([0 (ichan+1)*h]);
                
                %average freq values for the whole frequency band
                cfgtemp                     = [];
                cfgtemp.frequency           = 'all';
                cfgtemp.avgoverfreq         = 'yes';
                timefreq_alldata{ifreq}{itrial} = ft_selectdata(cfgtemp, timefreq_alldata{ifreq}{itrial});
                %plot(timefreq_alldata{ifreq}{itrial}.time, squeeze(timefreq_alldata{ifreq}{itrial}.powspctrm(:,1,:))); ylim([0 10]); xlim([-900 0])
                
                %go trough each channel to output results, as there will not
                %be the same amout of data between each channel (WOD/WOR do not
                %occur at the same time) :
                for ichan = 1:size(timefreq_alldata{ifreq}{itrial}.label,1)
                    
                    %select channel
                    ichan_name              = timefreq_alldata{ifreq}{itrial}.label{ichan};
                    cfgtemp                 = [];
                    cfgtemp.channel         = ichan_name;
                    timefreq_ichan_temp   	= ft_selectdata(cfgtemp,timefreq_alldata{ifreq}{itrial});
                    
                    %% RECOVERY DATA : Change T0 from Vent_Off to Vent_On
                    %MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_On.synctime(itrial);
                    timeshift                                           = config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2) + config{irat}.epoch.pad.(config{irat}.LFP.name{1}) - timefreq_ichan_temp.time(end);%MuseStruct{1}{1}.markers.Vent_On.synctime(itrial) - MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial);
                    timefreq_recovery{ifreq}{itrial}.(ichan_name)       = timefreq_ichan_temp;
                    timefreq_recovery{ifreq}{itrial}.(ichan_name).time 	= timefreq_ichan_temp.time + timeshift;
                    
                    %% WOD DATA : find WOD peak per channel, and normalize time per channel
                    %use filtered data to find wod

                    timefreq_wod{ifreq}{itrial}.(ichan_name)            = timefreq_ichan_temp;
                    
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
                    [~, t_peak_wod] = findpeaks(-LFP_trial_filt.trial{1}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
                    
                    %keep only data between 0 and wod
                    cfgtemp                                   = [];
                    cfgtemp.latency                           = [0 t_peak_wod];
                    timefreq_wod{ifreq}{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wod{ifreq}{itrial}.(ichan_name));
                    
                    %normalize time
                    timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name)        = timefreq_wod{ifreq}{itrial}.(ichan_name);
                    timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).time   = timefreq_wod{ifreq}{itrial}.(ichan_name).time ./ t_peak_wod;
                    %check the location of the peak detection
%                                     figure;hold
%                                     plot(LFP_trial_filt.time{1}./t_peak_wod, LFP_trial_filt.trial{1}(ichan,:));
%                                     xlim([0 2]);
%                                     xlim([0.8 1.2]);
                    
                    %resample data to have the same number of data points, for
                    %time-normalized data
                    t_old                                                = timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).time;
                    t_new                                                = linspace(0,1,1000);
                    powspctrm_new(1,1,:)                                 = pchip(t_old,squeeze(timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
                    timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).time         = t_new;
                    timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name).powspctrm    = powspctrm_new;
                    
                    %plot(timefreq_wod{ifreq}{itrial}.time{1},timefreq_wod{ifreq}{itrial}.trial{1}); xlim([0 0.95])
                    
                    %% WOR DATA : find WOR peak per channel, and normalize time per channel
                    %Antoine
                    %+ normaliser le temps que selon la 1ere wod/wor, pour toutes les électrodes du rat ? pour
                    %mieux voir le délai entre électrodes
                    
                    %% remove cfg fields as it is what takes the most of place on disk, whereas we do not use it later
                    timefreq_wod{ifreq}{itrial}.(ichan_name)             = rmfield(timefreq_wod{ifreq}{itrial}.(ichan_name),{'cfg'});
                    timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm{ifreq}{itrial}.(ichan_name),{'cfg'});
                    timefreq_recovery{ifreq}{itrial}.(ichan_name)      	 = rmfield(timefreq_recovery{ifreq}{itrial}.(ichan_name),{'cfg'});
                    
                end %ichan
            end
        end
        
        % add empty missing channels as nan channels to have the same channels between rats
        for ifreq = 1 : size(timefreq_recovery,2)
            for itrial = 1:size(timefreq_recovery{ifreq},2)
                chan_list = fieldnames(timefreq_recovery{ifreq}{itrial});
                for ichan = 1:n_chans-1
                    chan_name = ['E',num2str(n_chans-ichan)];
                    if ~any(strcmp(chan_list, chan_name))                       
                        timefreq_recovery{ifreq}{itrial}.(chan_name)            = [];
                        timefreq_wod{ifreq}{itrial}.(chan_name)                 = [];
                        timefreq_wod_timenorm{ifreq}{itrial}.(chan_name)        = [];
                    end
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
        %WOR data : time normalized between vent_off and wor peak, per channel. power normalized with baseline period :
        
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
    timefreq_recovery_allrats       = cell(1,size(cfgcommon.timefreq.foi,2));
    timefreq_wod_allrats            = cell(1,size(cfgcommon.timefreq.foi,2));
    timefreq_wod_timenorm_allrats 	= cell(1,size(cfgcommon.timefreq.foi,2));
    count_trials                    = 0;
    
    % load timefreq data from all rats
    for irat = [4 5 6 7 9 10 11] %4:size(config,2)
        
        %load data
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_recovery.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod.mat']));
        load(fullfile(config{irat}.datasavedir,[config{irat}.prefix, 'timefreq_wod_timenorm.mat']));
        
        %pool all rat data in the same structure
        for itrial = 1:size(timefreq_recovery{1},2)
            count_trials = count_trials + 1;
            for ifreq = 1:size(config{irat}.timefreq.foi,2)
                chan_list                                   = fieldnames(timefreq_wod{ifreq}{itrial});
                for ichan = 1:numel(chan_list)
                    chan_name                                                       = chan_list{ichan};
                    timefreq_recovery_allrats{ifreq}.(chan_name){count_trials}     	= timefreq_recovery{ifreq}{itrial}.(chan_name);
                    timefreq_wod_allrats{ifreq}.(chan_name){count_trials}           = timefreq_wod{ifreq}{itrial}.(chan_name);
                    timefreq_wod_timenorm_allrats{ifreq}.(chan_name){count_trials} 	= timefreq_wod_timenorm{ifreq}{itrial}.(chan_name);
                end
            end
        end
        
        clear timefreq_recovery timefreq_wod timefreq_wod_timenorm
    end
    
    % pool data for each frequency band and channel
    % if not the same number of points, keep only the number of points of
    % the shorter trial
    for ifreq = 1:size(config{irat}.timefreq.foi,2)
        chan_list                                   = fieldnames(timefreq_wod_timenorm_allrats{ifreq});
        for ichan = 1:numel(chan_list)
            chan_name                               = chan_list{ichan};
            
            %replace empty channels by nan channels
            hasdata = find(~cellfun(@isempty,timefreq_wod_timenorm_allrats{ifreq}.(chan_name)));
            for iwod = 1:size(timefreq_wod_allrats{ifreq}.(chan_name),2)
                if isempty(timefreq_wod_allrats{ifreq}.(chan_name){iwod})
                    timefreq_wod_allrats{ifreq}.(chan_name){iwod}                       = timefreq_wod_allrats{ifreq}.(chan_name){hasdata(1)};
                    timefreq_wod_allrats{ifreq}.(chan_name){iwod}.powspctrm             = nan(size(timefreq_wod_allrats{ifreq}.(chan_name){iwod}.powspctrm));
                    timefreq_recovery_allrats{ifreq}.(chan_name){iwod}                  = timefreq_recovery_allrats{ifreq}.(chan_name){hasdata(1)};
                    timefreq_recovery_allrats{ifreq}.(chan_name){iwod}.powspctrm        = nan(size(timefreq_recovery_allrats{ifreq}.(chan_name){iwod}.powspctrm));
                    timefreq_wod_timenorm_allrats{ifreq}.(chan_name){iwod}            	= timefreq_wod_timenorm_allrats{ifreq}.(chan_name){hasdata(1)};
                    timefreq_wod_timenorm_allrats{ifreq}.(chan_name){iwod}.powspctrm    = nan(size(timefreq_wod_timenorm_allrats{ifreq}.(chan_name){iwod}.powspctrm));
                end
            end
            
            cfgtemp                             = [];
            cfgtemp.keepindividual              = 'yes';
            timefreq_wod_grandavg{ifreq}.(chan_name)                = ft_freqgrandaverage(cfgtemp, timefreq_wod_allrats{ifreq}.(chan_name){:});
            timefreq_wod_timenorm_grandavg{ifreq}.(chan_name)     	= ft_freqgrandaverage(cfgtemp, timefreq_wod_timenorm_allrats{ifreq}.(chan_name){:});
            
            cfgtemp.toilim                      = [0 config{irat}.epoch.toi.(config{irat}.LFP.name{1})(2)];
            timefreq_recovery_grandavg{ifreq}.(chan_name)           = ft_freqgrandaverage(cfgtemp, timefreq_recovery_allrats{ifreq}.(chan_name){:});

        end
    end
    
    clear timefreq_recovery_allrats timefreq_wod_allrats timefreq_wod_timenorm_allrats
    
    %% compute cluster-based permutation statistics
    
    
    %% plot all data
    if ~isfolder(cfgcommon.imagesavedir)
        fprintf('Creating directory %s\n', cfgcommon.imagesavedir);
        mkdir(cfgcommon.imagesavedir)
    end
    
    analysis_names = {'Timefreq_WOD', 'Timefreq_WOD_timenorm', 'Timefreq_recovery'};
    data           = {timefreq_wod_grandavg, timefreq_wod_timenorm_grandavg, timefreq_recovery_grandavg};
    
    for idata = 1:size(data,2)
        fig = figure;
        sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(cfgcommon.timefreq.foi,2)
            
            subplot(2,2,ifreq);hold;
            chan_list 	= fieldnames(timefreq_wod_grandavg{ifreq});% donne le nom des champs de ma variable, ici le nom des electrodes
            C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
            
            for ichan = 1:numel(chan_list)
                
                chan_nr          = numel(chan_list)-ichan+1;% recréé le nom de l'électrode pour éviter un mélange des électrodes
                chan_label{ichan} = ['E', num2str(chan_nr)];
                %                 if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                %                     continue
                %                 end
                color   = C(chan_nr,:);
                
                x = data{idata}{ifreq}.(chan_label{ichan}).time;
                y = nanmean(squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(:,1,1,:)));
                
                if idata == 3 %only for recovery
                    y = movmean(y,10,'omitnan');
                end
                
                %plot std
                %                 std = nanstd(squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(:,1,1,:)));
                %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
                %                 filled_SD = area(x,y_area);
                %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
                %                 filled_SD(1).ShowBaseLine = 'off';
                
                %find maximum of the first half of data to adapt y scale
                data_max(ichan) = max(y(x>0 & x<x(end)/2)+std(x>0 & x<x(end)/2));
                
                %plot avg
                leg{ichan} = plot(x, y, 'Color', color);
            end
            
            %set figure display :
            %set(gca, 'YScale', 'log');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi{ifreq}(1),config{irat}.timefreq.foi{ifreq}(end)));
            legend([leg{:}], chan_label{:}, 'Fontsize',8,'location','eastoutside');
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
        
    end
    
    close all
    
    % plot data by protocol
    for idata = 1:size(data,2)
        for iwod = 1 :size(data{idata}{ifreq}.(chan_label{ichan}).powspctrm,1)
            
            fig = figure;
            sgtitle(sprintf('%s WOD%g',analysis_names{idata},iwod), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
            
            for ifreq = 1:size(cfgcommon.timefreq.foi,2)
               
                subplot(2,2,ifreq);hold;
                chan_list 	= fieldnames(timefreq_wod_grandavg{ifreq});% donne le nom des champs de ma variable, ici le nom des electrodes
                C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
                
                for ichan = 1:numel(chan_list)
                    chan_nr          = numel(chan_list)-ichan+1;% recréé le nom de l'électrode pour éviter un mélange des électrodes
                    chan_label{ichan} = ['E', num2str(chan_nr)];
                    %                 if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                    %                     continue
                    %                 end
                    color   = C(chan_nr,:);
                    
                    x = data{idata}{ifreq}.(chan_label{ichan}).time;
                    y = squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(iwod,1,1,:));
                    
                    if idata == 3 %only for recovery
                        y = movmean(y,10,'omitnan');
                    end
                    
                    %plot std
                    %                 std = nanstd(squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(:,1,1,:)));
                    %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
                    %                 filled_SD = area(x,y_area);
                    %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                    %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                    %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
                    %                 filled_SD(1).ShowBaseLine = 'off';
                    
                    %find maximum of the first half of data to adapt y scale
                    data_max(ichan) = max(y(x>0 & x<x(end)/2)+std(x>0 & x<x(end)/2));
                    
                    %plot avg
                    leg{ichan} = plot(x, y, 'Color', color);
                end
                
                %set figure display :
                %set(gca, 'YScale', 'log');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time');
                ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi{ifreq}(1),config{irat}.timefreq.foi{ifreq}(end)));
                legend([leg{:}], chan_label{:}, 'Fontsize',8,'location','eastoutside');
                legend('boxoff');
                axis tight;
                if contains(analysis_names{idata}, 'WOD')
                    ylim([0 max(data_max)]);
                end
                
            end
            handle_plot(1)= fig
            %save figure :
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,sprintf('%s_WOD%g.pdf',analysis_names{idata},iwod)),'-r600');
            print(fig, '-dpng', fullfile(cfgcommon.imagesavedir,sprintf('%s_WOD%g.png',analysis_names{idata},iwod)),'-r600');
            savefig(fig,fullfile(cfgcommon.imagesavedir,sprintf('%s_WOD%g',analysis_names{idata},iwod)));
            
            close all
        end
    end
    
    % compute statistiques : y a-t-il des profondeurs qui récupèrent plus vite ? Lesquelles ?
        
    % trouver période stable de recovery
    
end %slurm_task_id == 0

end


