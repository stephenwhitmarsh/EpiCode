function dtx_Spiky(slurm_task_id)

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\SPIKY_2_Dez_2019\SPIKY'))
    cd \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\SPIKY_2_Dez_2019\SPIKY\SPIKY_MEX
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/SPIKY_2_Dez_2019/SPIKY'))
    cd /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/SPIKY_2_Dez_2019/SPIKY/SPIKY_MEX
    
end

ft_defaults

% %This is the first file that should be run once the zip-package has been extracted.
% %Once the MEX-files have been compiled you can run the main program SPIKY.
% SPIKY_compile_MEX
% %GUI
% test_Spiky = spikes{itrial};
% SPIKY;

%Setting parameters
config      = dtx_setparams_probe_spikes([]);
ipart       =1;
ilabel      = 4; %cfg.name{4} = 'SlowWave_Larger' : need to have a long period of padding

for irat = slurm_task_id
    
    MuseStruct      = readMuseMarkers(config{irat}, false);
    MuseStruct      = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true); 
    SpikeTrials     = readSpikeTrials_MuseMarkers(config{irat}, [],[], false);
    
    %remove noisy units
    units_infos     = dtx_read_unit_table(config{irat},SpikeTrials);
    for i_unit = 1:size(units_infos.label, 2)
        if contains(units_infos.group{i_unit}, 'noise')
            toremove(i_unit) = 1;
        else
            toremove(i_unit) = 0;
        end
    end
    cfgtemp = [];
    cfgtemp.spikechannel = SpikeTrials{ipart}{ilabel}.label(~toremove);
    SpikeTrials{ipart}{ilabel} = ft_spike_select(cfgtemp, SpikeTrials{ipart}{ilabel});
    
    %remove trials with too short interictal periods
    cfgtemp                         = config{irat};
    try cfgtemp = rmfield(cfgtemp, 'remove'); end %remove 'remove' field if it exists
    cfgtemp.remove.markerstart      = 'Crise_End';
    cfgtemp.remove.markerend        = 'Crise_End';
    cfgtemp.remove.timefrombegin    = -80; %-80 aucune chance de commencer avant le trial d'avant (-120). Et si complètement inclu, alors il est gardé, mais ça voudra dire que la distance est suffisante
    cfgtemp.remove.timefromend      = -50; %-50 : time before SW min 70s (SW-120)
    cfgtemp.remove.label_list       = ilabel;
    cfgtemp.remove.searchdirection  = 'before';
    cfgtemp.remove.write            = 'no';
    [SpikeTrials, ~]                = removetrials_MuseMarkers(cfgtemp, SpikeTrials, MuseStruct);
    
    %create spike structure to use with spiky
    for i_unit = 1:size(SpikeTrials{1}{ilabel}.label, 2)
        for itrial = 1:size(SpikeTrials{1}{ilabel}.trialinfo, 1)
            idx                         = SpikeTrials{1}{ilabel}.trial{i_unit} == itrial;
            spikedata{itrial}{i_unit}   = SpikeTrials{1}{ilabel}.time{i_unit}(idx);
        end
    end
    
    %control : merge units of different trials
    %each unit come from a different trial. There are not several units of
    %the same trial. Except if there are not enough trials
    for itrial = 1:length(spikedata)
        trial_list = 1:length(spikedata);
        trial_list(trial_list==itrial) = []; %remove the current trial of the random selection of trials
        for i_unit = 1:length(spikedata{itrial})
            if isempty(trial_list) %if more trials than units, start again the random selection (usefull for 1 rat from 5)
                trial_list = 1:length(spikedata);
                trial_list(trial_list==itrial) = [];
            end
            rand_trial = trial_list(randi(length(trial_list)));
            trial_list(trial_list==rand_trial) = []; %remove the selected trial for the next unit
            spikedata_ctrl{itrial}{i_unit} = spikedata{rand_trial}{i_unit};
        end
    end
            
    % remove empty units in each trial, as synchrony would be artificially maximal with those
    for itrial = 1:length(spikedata)
        spikedata{itrial}           = spikedata{itrial}(~cellfun('isempty', spikedata{itrial}));
        spikedata_ctrl{itrial}      = spikedata_ctrl{itrial}(~cellfun('isempty', spikedata_ctrl{itrial}));
    end
    
    %remove if interictal trials <70s (ilabel = 3) 
%     triallength = (SpikeTrials{1}{3}.trialtime(:, 2) - SpikeTrials{1}{3}.trialtime(:, 1));
%     allspikes   = allspikes((triallength>70)'); %NO NEED ANYMORE? USE SCRIPT RMTRIALS
    trialtime   = SpikeTrials{1}{ilabel}.trialtime;
    
    % for i_unit = 1:size(SpikeTrials{1}{ilabel}.label, 2)
    %     itrial = 100;
    %     spike_indx_onetrial = SpikeTrials{1}{ilabel}.trial{i_unit} == itrial;
    %     test_spiky{i_unit} = SpikeTrials{1}{ilabel}.time{i_unit}(spike_indx_onetrial);
    % end
    

    
    %% spiky loop
    %define parameters
    ft_progress('init','text', config{irat}.prefix(1:end-1));
    for itrial = 1:length(spikedata)
        
        ft_progress(itrial/length(spikedata), 'processing trial %d from %d', itrial, length(spikedata));
        
        % For real data
        spikes = spikedata{itrial};
        para.tmin = trialtime(itrial,1);
        para.tmax = trialtime(itrial,2);
        para.dts  = 1/SpikeTrials{1}{ilabel}.hdr.Fs;
        para.select_measures      =[0 1 0 0 0 0 1];  % Select measures (0-do not calculate,1-calculate)
        m_para.all_measures_string={'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_forward';'SPIKE_synchro';'SPIKE_order';'PSTH'};  % order of select_measures
        para.num_trains=length(spikes);
        
        %automatic setting of params
        d_para=para;
        SPIKY_check_spikes
%         if ret==1
%             continue
%         end
        para=d_para;
        
        trialdata.results{itrial}      = SPIKY_loop_f_distances(spikes,para); %FIXME une erreur corrigée
        
        % for control data with trials merged
        %same 'para'
        spikes = spikedata_ctrl{itrial};
        d_para=para;
        SPIKY_check_spikes
%         if ret==1
%             continue
%         end
        para=d_para;
        trialdata.results_ctrl_trials{itrial}      = SPIKY_loop_f_distances(spikes,para); %FIXME une erreur corrigée
        
        % for control data with same psth
        para.choice=4;%1;%spike number;%4; %select psth conservation
        ori_spikes=spikes;
        spikes=SPIKY_f_spike_train_surrogates(ori_spikes,para);
        trialdata.results_ctrl_psth{itrial} = SPIKY_loop_f_distances(spikes,para);
        
        clear spikes ori_spikes 
    end
    ft_progress('close');
    
    %gather results 
    for idata = ["results","results_ctrl_trials","results_ctrl_psth"]
        for itrial = 1:length(spikedata)
            spikyresults.(idata).profile_time{itrial}        = trialdata.(idata){itrial}.SPIKE.time;% - results{itrial}.SPIKE.time(end); %timelock to the end
            spikyresults.(idata).profile_values{itrial}      = trialdata.(idata){itrial}.SPIKE.profile;
            spikyresults.(idata).psth_time{itrial}           = trialdata.(idata){itrial}.PSTH.time;
            spikyresults.(idata).psth_values{itrial}         = trialdata.(idata){itrial}.PSTH.profile;
        end
    end
    clear trialdata
    
    % SPIKY_f_average_pi : moyenne sur plusieurs profils de dissimilarité.
    for idata = ["results","results_ctrl_trials","results_ctrl_psth"]
        fprintf('average dissimilarity profiles for %s\n', idata);
        [spikyresults.(idata).profileavg_x, spikyresults.(idata).profileavg_y] = SPIKY_f_average_pi(spikyresults.(idata).profile_time,spikyresults.(idata).profile_values,para.dts);
        fprintf('average psth for %s\n', idata);
        [spikyresults.(idata).psthavg_x, spikyresults.(idata).psthavg_y] = SPIKY_f_average_pi(spikyresults.(idata).psth_time,spikyresults.(idata).psth_values,para.dts);
    end
%        
%     figure;hold;
%     plot(x_avg_ctrl,y_avg_ctrl);
%     plot(x_avg,y_avg);
%     plot(x_psth,y_psth);
%     plot(x_psth_ctrl,y_psth_ctrl);
    
    fprintf('save data\n');
  
    fname_out = fullfile(config{irat}.datasavedir,[config{irat}.prefix,'spiky_results.mat']);
    save(fname_out,'spikyresults','-v7.3');
    
    clear SpikeTrials spikes profile* x* y* results para allspikes triallength spikyresults
end
return 

%% load data
rat_list = 1:5;

for irat = rat_list
    fprintf('Load spiky results of %s\n', config{irat}.prefix(1:end-1));
    fname = fullfile(config{irat}.datasavedir,[config{irat}.prefix,'spiky_results.mat']);
    load(fname);
    spiky_results{irat} = spikyresults;
    clear spikyresults;
end

%plot spike distance,for the 3 kinds of data
figure;
subplot(3,1,1);hold;
%plot ctrl psth data
for irat = rat_list
    x = spiky_results{irat}.results_ctrl_psth.profileavg_x;
    y = spiky_results{irat}.results_ctrl_psth.profileavg_y;
    ctrlpsth = plot(x,y, 'r');
end
%plot ctrl trials data
for irat = rat_list
    x = spiky_results{irat}.results_ctrl_trials.profileavg_x;
    y = spiky_results{irat}.results_ctrl_trials.profileavg_y;
    ctrltrials = plot(x,y, 'g');
end
%plot real data
for irat = rat_list
    x = spiky_results{irat}.results.profileavg_x;
    y = spiky_results{irat}.results.profileavg_y;
    dtx = plot(x,y, 'k');
end
legend([dtx ctrltrials ctrlpsth], 'Data', 'Trials schuffled', 'PSTH shuffled');
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('Time (s)');
xlim([-10 10]);
% xlim([-60 10]);

%plot psth, for the 3 kinds of data
subplot(3,1,2);hold;
%plot ctrl psth data
for irat = rat_list
    x = spiky_results{irat}.results_ctrl_psth.psthavg_x;
    y = spiky_results{irat}.results_ctrl_psth.psthavg_y;
    ctrlpsth = plot(x,y, 'r');
end
%plot ctrl trials data
for irat = rat_list
    x = spiky_results{irat}.results_ctrl_trials.psthavg_x;
    y = spiky_results{irat}.results_ctrl_trials.psthavg_y;
    ctrltrials = plot(x,y, 'g');
end
%plot real data
for irat = rat_list
    x = spiky_results{irat}.results.psthavg_x;
    y = spiky_results{irat}.results.psthavg_y;
    dtx = plot(x,y, 'k');
end
legend([dtx ctrltrials ctrlpsth], 'Data', 'Trials schuffled', 'PSTH shuffled');
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('Time (s)');
xlim([-10 10]);
% xlim([-60 10]);



% [dat_LFP ~] = removetrials_MuseMarkers([],dat_LFP, MuseStruct);
% %FIXME : load LFP before
% subplot(3,1,3);hold;
% cfgtemp = [];
% cfgtemp.channame = 'E12LFP';
% plot_morpho(cfgtemp, dat_LFP{1}{1});
% xlim([-60 10]);   
% title([]);
% xlim([-2 2])



%% plot avg over rats
% % find max x nr for interpolation
% spikyval_int = [];
% spikyctrl_int = [];
% for irat = rat_list
%     %remove dupplicated abscisses
%     [spiky_time, idx, ~]        = unique(spiky_results{irat}.spiky.spike_distance_results.time_avg);
%     spiky_val                   = spiky_results{irat}.spiky.spike_distance_results.profil_avg(idx);
%     spikytime_int               = linspace(spiky_time(1),spiky_time(end),10000);
%     spikyval_int(irat,:)        = pchip(spiky_time,spiky_val,spikytime_int);
%     
%     %same for ctrl
%     [spikyctrl_time, idx, ~]    = unique(spiky_results{irat}.spiky.spike_distance_control.time_avg);
%     spikyctrl_val               = spiky_results{irat}.spiky.spike_distance_control.profil_avg(idx);
%     spikytime_int               = linspace(spikyctrl_time(1),spikyctrl_time(end),10000);
%     spikyctrl_int(irat,:)       = pchip(spikyctrl_time,spikyctrl_val,spikytime_int);
% end
% 
% plot(spikytime_int, nanmean(spikyval_int,1), spikytime_int, nanmean(spikyctrl_int,1));

%% examples with the different kinds of results
% SPIKY_loop
% 
% %Spike distance
% imagesc(SPIKY_results.SPIKE.matrix);
% c = colorbar;
% c.Limits = [0 1];
% plot(SPIKY_results.SPIKE.time, SPIKY_results.SPIKE.profile);
% plot(results{itrial}.SPIKE.time, results{itrial}.SPIKE.profile);
% 
% %Spiek synchro
% fig = imagesc((SPIKY_results.SPIKE_synchro.matrix));
% c = colorbar;
% c.Limits = [0 1];
% colormap(flipud(parula));
% colorbar('Direction', 'reverse');
% 
% hold;
% x = SPIKY_results.SPIKE_synchro.time;
% y = SPIKY_results.SPIKE_synchro.profile;
% plot(x, y);
% 
% y = movmean(y,1000);
% plot(x, y);
% 
