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

% This is the first file that should be run once the zip-package has been extracted.
% Once the MEX-files have been compiled you can run the main program SPIKY.
% SPIKY_compile_MEX
config = dtx_setparams_probe_spikes([]);

ipart=1;
for irat = 1:5
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    MuseStruct = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true); 
    SpikeTrials = readSpikeTrials_MuseMarkers(config{irat}, [],[], false);
    cfgtemp = [];
    cfgtemp.remove.markerstart = 'Crise_End';
    cfgtemp.remove.markerend = 'Crise_End';
    cfgtemp.remove.timefrombegin = -80; %-80 aucune chance de commancer avant le trial d'avant (-120). Et si complètement inclu, alors il est gardé, mais ça vodra dire que la distance est suffisante
    cfgtemp.remove.timefromend = -50; %-50 : time before SW min 70s (SW-120)
    cfgtemp.remove.label_list = 4;
    cfgtemp.remove.searchdirection = 'before';
    [SpikeTrials, ~]                = removetrials_MuseMarkers(cfgtemp, SpikeTrials, MuseStruct);
    %create spike structure to use with spiky
    %control : add a random number between -2.5 and 2.5 to each spike idx
    for i_unit = 1:size(SpikeTrials{1}{4}.label, 2)
        for itrial = 1:size(SpikeTrials{1}{4}.trialinfo, 1)
            idx = SpikeTrials{1}{4}.trial{i_unit} == itrial;
            allspikes{itrial}{i_unit} = SpikeTrials{1}{4}.time{i_unit}(idx);
        end
    end
    % remove empty trials, see if results are better
    for itrial = 1:length(allspikes)
        allspikes{itrial}   = allspikes{itrial}(~cellfun('isempty', allspikes{itrial}));
    end
    
%     %remove if interictal trials <70s (ilabel = 3) 
    triallength = (SpikeTrials{1}{3}.trialtime(:, 2) - SpikeTrials{1}{3}.trialtime(:, 1));
%     allspikes   = allspikes((triallength>70)'); %NO NEED ANYMORE? USE SCRIPT RMTRIALS
    trialtime   = SpikeTrials{1}{4}.trialtime;
    
    % for i_unit = 1:size(SpikeTrials{1}{4}.label, 2)
    %     itrial = 100;
    %     spike_indx_onetrial = SpikeTrials{1}{4}.trial{i_unit} == itrial;
    %     test_spiky{i_unit} = SpikeTrials{1}{4}.time{i_unit}(spike_indx_onetrial);
    % end
    
    %GUI
    % test_Spiky = spikes{itrial};
    % SPIKY;
    
    %% loop
    %define parameters
    clear waitbar
    for itrial = 1:length(allspikes)
        try eval(['fprintf(repmat(''\b'', 1, length(waitbar)));']); end %error the first time. \b to remove previous text
        waitbar = sprintf('Trial %d/%d\n', itrial, length(allspikes));
        fprintf('%s', waitbar);
%         fprintf('Trial %d/%d\n', itrial, length(allspikes));
        spikes = allspikes{itrial};
        para.tmin = trialtime(itrial,1);
        para.tmax = trialtime(itrial,2);
        para.dts  = 1/SpikeTrials{1}{4}.hdr.Fs;
        para.select_measures      =[0 1 0 0 0 0 1];  % Select measures (0-do not calculate,1-calculate)
        m_para.all_measures_string={'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_forward';'SPIKE_synchro';'SPIKE_order';'PSTH'};  % order of select_measures
        para.num_trains=length(spikes);
        
        %automatic setting of params
        d_para=para;
        SPIKY_check_spikes
        if ret==1
            return
        end
        para=d_para;
        
        results{itrial}      = SPIKY_loop_f_distances(spikes,para); %FIXME une erreur corrigée
        
        %compute control values
        para.choice=4;%1;%spike number;%4; %select psth conservation
        ori_spikes=spikes;
        spikes=SPIKY_f_spike_train_surrogates(ori_spikes,para);
        results_ctrl{itrial} = SPIKY_loop_f_distances(spikes,para);
        
        clear spikes ori_spikes para
    end
    
    %gather results 
    for itrial = 1:length(allspikes)
        profile_time{itrial}        = results{itrial}.SPIKE.time;% - results{itrial}.SPIKE.time(end); %timelock to the end
        profile_values{itrial}      = results{itrial}.SPIKE.profile;        
        profile_time_ctrl{itrial}   = results_ctrl{itrial}.SPIKE.time;% - results_ctrl{itrial}.SPIKE.time(end); %timelock to the end
        profile_values_ctrl{itrial} = results_ctrl{itrial}.SPIKE.profile;
        psth_time{itrial}           = results{itrial}.PSTH.time;
        psth_values{itrial}         = results{itrial}.PSTH.profile;
        psth_time_ctrl{itrial}      = results_ctrl{itrial}.PSTH.time;
        psth_values_ctrl{itrial}    = results_ctrl{itrial}.PSTH.profile;
    end
    
    % SPIKY_f_average_pi : moyenne sur plusieurs profils de dissimilarité.
    fprintf('average dissimilarity profiles for results\n');
    [x_avg, y_avg] = SPIKY_f_average_pi(profile_time,profile_values,1/32000);
    fprintf('average dissimilarity profiles for results_ctrl\n');
    [x_avg_ctrl, y_avg_ctrl] = SPIKY_f_average_pi(profile_time_ctrl,profile_values_ctrl,1/32000);
    fprintf('average dissimilarity profiles for psth\n');
    [x_psth, y_psth] = SPIKY_f_average_pi(psth_time,psth_values,1/32000);
    fprintf('average dissimilarity profiles for psth\n');
    [x_psth_ctrl, y_psth_ctrl] = SPIKY_f_average_pi(psth_time_ctrl,psth_values_ctrl,1/32000);
    
    figure;hold;
    plot(x_avg_ctrl,y_avg_ctrl);
    plot(x_avg,y_avg);
    plot(x_psth,y_psth);
    plot(x_psth_ctrl,y_psth_ctrl);
    
    fprintf('save data\n');
    
    spiky.spike_distance_results.time         = profile_time;
    spiky.spike_distance_results.profile      = profile_values;
    spiky.spike_distance_results.time_avg     = x_avg;
    spiky.spike_distance_results.profil_avg   = y_avg;    
    spiky.spike_distance_control.time         = profile_time_ctrl;
    spiky.spike_distance_control.profile      = profile_values_ctrl;
    spiky.spike_distance_control.time_avg     = x_avg_ctrl;
    spiky.spike_distance_control.profil_avg   = y_avg_ctrl;
    
    spiky.psth_results.time         = psth_time;
    spiky.psth_results.profile      = psth_values;
    spiky.psth_results.time_avg     = x_psth;
    spiky.psth_results.profil_avg   = y_psth;
    spiky.psth_control.time         = psth_time_ctrl;
    spiky.psth_control.profile      = psth_values_ctrl;
    spiky.psth_control.time_avg     = x_psth_ctrl;
    spiky.psth_control.profil_avg   = y_psth_ctrl;
    
    fname_out = fullfile(config{irat}.datasavedir,[config{irat}.prefix,'spiky_results.mat']);
    save(fname_out,'spiky','-v7.3');
    
    clear SpikeTrials spikes profile* x* y* results para allspikes triallength spiky
end
error('stop here with the cluster'); %for cluster

%% load data
for irat = 1:5
    fprintf('Load spike distance for %s\n', config{irat}.prefix(1:end-1));
    fname = fullfile(config{irat}.datasavedir,[config{irat}.prefix,'spiky_results.mat']);
    spiky_results{irat} = load(fname);
end

figure;hold;
C = linspecer(5);
for irat = 1:5
    for itrial = 1:length(spiky_results{irat}.spiky.spike_distance_results.time)
        x = spiky_results{irat}.spiky.spike_distance_results.time{itrial};
        y = spiky_results{irat}.spiky.spike_distance_results.profile{itrial};
        ctrl = plot(x,y, 'Color', C(irat,:));
    end
end


figure;hold;
C = linspecer(5);
for irat = 1:5
    for itrial = 1:length(spiky_results{irat}.spiky.spike_distance_results.time)
        x = spiky_results{irat}.spiky.spike_distance_control.time{itrial};
        y = spiky_results{irat}.spiky.spike_distance_control.profile{itrial};
        ctrl = plot(x,y, 'Color', C(irat,:));
    end
end

% test fieldtrip psth
% cfgtemp = [];
% cfgtemp.binsize      = 1/50;
% cfg.keeptrials       = 'yes';
% test_psth = ft_spike_psth(cfgtemp,SpikeTrials{1}{1});
% 
% cfgtemp = [];
% cfgtemp.errorbars        = 'std';
% cfgtemp.spikechannel     = 1;
% ft_spike_plot_psth(cfgtemp, test_psth);


figure;
subplot(3,1,1);hold;
%plot ctrl data
for irat = 1:5
    x = spiky_results{irat}.spiky.spike_distance_control.time_avg;
    y = spiky_results{irat}.spiky.spike_distance_control.profil_avg;
    ctrl = plot(x,y, 'r');
end
%plot data
for irat = 1:5
    x = spiky_results{irat}.spiky.spike_distance_results.time_avg;
    y = spiky_results{irat}.spiky.spike_distance_results.profil_avg;
    dtx = plot(x,y, 'k');
end
xlim([-60 10]);

subplot(3,1,2);hold;
%plot psth data
for irat = 1:5
    x = spiky_results{irat}.spiky.psth_results.time_avg;
    y = spiky_results{irat}.spiky.psth_results.profil_avg;
    dtx = plot(x,y, 'k');
end
%plot ctrl psth
for irat = 1:5
    x = spiky_results{irat}.spiky.psth_control.time_avg;
    y = spiky_results{irat}.spiky.psth_control.profil_avg;
    dtx = plot(x,y, 'r');
end
xlim([-60 10]);
    
[dat_LFP ~] = removetrials_MuseMarkers([],dat_LFP, MuseStruct);
%FIXME : load LFP before
subplot(3,1,3);hold;
cfgtemp = [];
cfgtemp.channame = 'E12LFP';
plot_morpho(cfgtemp, dat_LFP{1}{1});
xlim([-60 10]);   
title([]);
xlim([-2 2])

legend([dtx ctrl], 'Data', 'Schuffled');
% xlim([-60 -20]); %0 = SlowWave for ilabel = 4
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('Time (s)');
% xticklabels(xticks+20);
xlim([-10 10]);

%plot avg over rats
% find max x nr for interpolation
spikyval_int = [];
spikyctrl_int = [];
for irat = 1:5
    %remove dupplicated abscisses
    [spiky_time, idx, ~]        = unique(spiky_results{irat}.spiky.spike_distance_results.time_avg);
    spiky_val                   = spiky_results{irat}.spiky.spike_distance_results.profil_avg(idx);
    spikytime_int               = linspace(spiky_time(1),spiky_time(end),10000);
    spikyval_int(irat,:)        = pchip(spiky_time,spiky_val,spikytime_int);
    
    %same for ctrl
    [spikyctrl_time, idx, ~]    = unique(spiky_results{irat}.spiky.spike_distance_control.time_avg);
    spikyctrl_val               = spiky_results{irat}.spiky.spike_distance_control.profil_avg(idx);
    spikytime_int               = linspace(spikyctrl_time(1),spikyctrl_time(end),10000);
    spikyctrl_int(irat,:)       = pchip(spikyctrl_time,spikyctrl_val,spikytime_int);
end

plot(spikytime_int, nanmean(spikyval_int,1), spikytime_int, nanmean(spikyctrl_int,1));

%% examples with the different kinds of results
SPIKY_loop

%Spike distance
imagesc(SPIKY_results.SPIKE.matrix);
c = colorbar;
c.Limits = [0 1];
plot(SPIKY_results.SPIKE.time, SPIKY_results.SPIKE.profile);
plot(results{itrial}.SPIKE.time, results{itrial}.SPIKE.profile);

%Spiek synchro
fig = imagesc((SPIKY_results.SPIKE_synchro.matrix));
c = colorbar;
c.Limits = [0 1];
colormap(flipud(parula));
colorbar('Direction', 'reverse');

hold;
x = SPIKY_results.SPIKE_synchro.time;
y = SPIKY_results.SPIKE_synchro.profile;
plot(x, y);

y = movmean(y,1000);
plot(x, y);

