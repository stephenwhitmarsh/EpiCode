function dtx_spike_beforeafter_inj
%run it on the cluster because stats over time is too big
% statsovertime{ipart}{ilabel}.(param){i_unit}{i_trial}
% stats_concat.(param){ipart}{irat}{i_unit}(itrial,:)
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    
end

ft_defaults

config = dtx_setparams_probe_spikes([]);

rat_list.dtx     = [1 3 4];
rat_list.control = [6 7]; %only rats for which a baseline period was recorded
ipart       = 1;
after_inj_index      = 3;
before_inj_index      = 9;

for irat = [rat_list.dtx rat_list.control]
    fprintf('loading data for irat = %d\n',irat);

      %align markers and remove wrong seizures
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);

    if strcmp(config{irat}.type, 'dtx') 
        MuseStruct = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true); 
    end

    %read spike data
    SpikeRaw                        = readSpikeRaw_Phy(config{irat},false);
    if strcmp(config{irat}.type, 'dtx')
        %make trials based on Muse Markers
        SpikeTrials{irat}                     = readSpikeTrials_MuseMarkers(config{irat}, MuseStruct,SpikeRaw, false);
    elseif strcmp(config{irat}.type, 'ctrl')
        %read spike data
        spikedata_temp                        = readSpikeTrials_continuous(config{irat},SpikeRaw, false);
        %get muse timings
        [~,baseline_start]      = concatenateMuseMarkers(MuseStruct,1,'Baseline_Start');
        [~,injection]      = concatenateMuseMarkers(MuseStruct,1,'Injection');
        [~,analysis_end]      = concatenateMuseMarkers(MuseStruct,1,'Analysis_End');
        injection       = injection{1,1}*spikedata_temp{ipart}{1}.hdr.Fs;
        baseline_start  = baseline_start{1,1}*spikedata_temp{ipart}{1}.hdr.Fs;
        analysis_end    = analysis_end{1,1}*spikedata_temp{ipart}{1}.hdr.Fs;
        %select data before injection 
        cfgtemp = [];
        cfgtemp.trials = spikedata_temp{ipart}{1}.trialinfo(:,3)>baseline_start & spikedata_temp{ipart}{1}.trialinfo(:,4)<injection;
        SpikeTrials{irat}{ipart}{before_inj_index} = ft_spike_select(cfgtemp,spikedata_temp{ipart}{1});
        %select data after injection
        cfgtemp = [];
        cfgtemp.trials = spikedata_temp{ipart}{1}.trialinfo(:,3)>injection & spikedata_temp{ipart}{1}.trialinfo(:,4)<analysis_end;
        SpikeTrials{irat}{ipart}{after_inj_index} = ft_spike_select(cfgtemp,spikedata_temp{ipart}{1});
    end
    Spikeinfos{irat}                = dtx_read_unit_table(config{irat},SpikeTrials{irat}{ipart}{before_inj_index});
end

for experiment = ["dtx", "control"]
    freq_afterinj.(experiment) = [];
    freq_beforeinj.(experiment)  = [];
    cv2_afterinj.(experiment)  = [];
    cv2_beforeinj.(experiment)   = [];
    for irat = rat_list.(experiment)
        isi{irat}{after_inj_index} = ft_spike_isi([], SpikeTrials{irat}{ipart}{after_inj_index});
        isi{irat}{before_inj_index} = ft_spike_isi([], SpikeTrials{irat}{ipart}{before_inj_index});
        for i_unit = 1:size(Spikeinfos{irat}.label,2)
            %select only sua data
            if contains(Spikeinfos{irat}.group{i_unit}, 'noise')
                continue
            end
            if contains(Spikeinfos{irat}.group{i_unit}, 'mua')
                continue
            end
            
            %compute freq for afterinj
            trials_duration = SpikeTrials{irat}{ipart}{after_inj_index}.trialtime(:,2) - SpikeTrials{irat}{ipart}{after_inj_index}.trialtime(:,1);
            total_duration = sum(trials_duration);
            freq_afterinj.(experiment)(end+1) = size(SpikeTrials{irat}{ipart}{after_inj_index}.trial{i_unit},2) / total_duration;
            
            %compute cv2 for afterinj
            isitemp = isi{irat}{after_inj_index}.isi{i_unit};
            cv2_data = [];
            for i = 1:length(isitemp)-1
                cv2_data(i) = 2*abs(isitemp(i)-isitemp(i+1))/(isitemp(i)+isitemp(i+1));
            end
            cv2_afterinj.(experiment)(end+1) = nanmean(cv2_data);
            
            %compute freq for beforeinj
            trials_duration = SpikeTrials{irat}{ipart}{before_inj_index}.trialtime(:,2) - SpikeTrials{irat}{ipart}{before_inj_index}.trialtime(:,1);
            total_duration = sum(trials_duration);
            freq_beforeinj.(experiment)(end+1) = size(SpikeTrials{irat}{ipart}{before_inj_index}.trial{i_unit},2) / total_duration;
            
            %compute cv2 for beforeinj
            isitemp = isi{irat}{before_inj_index}.isi{i_unit};
            cv2_data = [];
            for i = 1:length(isitemp)-1
                cv2_data(i) = 2*abs(isitemp(i)-isitemp(i+1))/(isitemp(i)+isitemp(i+1));
            end
            cv2_beforeinj.(experiment)(end+1) = nanmean(cv2_data);
            
        end
    end
    
    if strcmp(experiment, "dtx")
        freq_beforeinj.(experiment) = freq_beforeinj.(experiment)(freq_beforeinj.(experiment)<10);
        freq_afterinj.(experiment) = freq_afterinj.(experiment)(freq_beforeinj.(experiment)<10);
        cv2_beforeinj.(experiment) = cv2_beforeinj.(experiment)(freq_beforeinj.(experiment)<10);
        cv2_afterinj.(experiment) = cv2_afterinj.(experiment)(freq_beforeinj.(experiment)<10);
    end
    
    figure;hold;
    for i_sua = 1:size(freq_afterinj.(experiment),2)
        plot([1 2], [freq_beforeinj.(experiment)(i_sua)  freq_afterinj.(experiment)(i_sua)], '-sk', 'MarkerFaceColor','k');
    end
    errorbar([1 2], [mean(freq_beforeinj.(experiment))  mean(freq_afterinj.(experiment))], [std(freq_beforeinj.(experiment))  std(freq_afterinj.(experiment))],'--r','LineWidth',2);
    xlim([0.8 2.2]);
    title(experiment);
    
    figure;hold;
    for i_sua = 1:size(cv2_afterinj.(experiment),2)
        plot([1 2], [cv2_beforeinj.(experiment)(i_sua)  cv2_afterinj.(experiment)(i_sua)], '-sk', 'MarkerFaceColor','k');
    end
    errorbar([1 2], [nanmean(cv2_beforeinj.(experiment))  nanmean(cv2_afterinj.(experiment))], [nanstd(cv2_beforeinj.(experiment))  nanstd(cv2_afterinj.(experiment))],'--r','LineWidth',2);
    xlim([0.8 2.2]);
    title(experiment)
    
    %% compute stats
    p_freq.(experiment) = signrank(freq_beforeinj.(experiment),freq_afterinj.(experiment));
    p_cv2.(experiment)  = signrank(cv2_beforeinj.(experiment),cv2_afterinj.(experiment));
end
        