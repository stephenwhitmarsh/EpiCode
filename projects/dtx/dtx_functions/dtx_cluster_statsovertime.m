function dtx_cluster_statsovertime
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

rat_list    = 1:5;
ipart       = 1;
baseline_index      = 3;

%load precomputed spike stats (dtx_spike_project)
for irat = rat_list
    fprintf('loading data for irat = %d\n',irat);
    %      config{irat}.datasavedir = fullfile(config{irat}.datasavedir, 'test_alignXCorr');
    %     stats{irat}                         = spikeratestats_Events_Baseline(config{irat},false);
    config{irat}.statstime.suffix       = [];
    statsovertime{irat}                 = spikestatsOverTime(config{irat},[],false);
    config{irat}.statstime.suffix       = '_withoutbursts';
    statsovertime_withoutbursts{irat}   = spikestatsOverTime(config{irat},[],false);
end
%% compute avg over trials of each stat for the last 60 seconds, for each unit
% statsovertime{ipart}{ilabel}.(param){i_unit}{i_trial}

%PAUL : fenetre glissante de 2 secondes, donc 60s = 30 fenetres

for irat = rat_list
    for i_unit = 1:size(statsovertime{irat}{ipart}{baseline_index}.freq, 2)
        for itrial = 1:size(statsovertime{irat}{ipart}{baseline_index}.freq{i_unit},2)
            t = statsovertime{irat}{ipart}{baseline_index}.time{i_unit}{itrial};
            %ignore trial if it is too small
            if t(1) > t(end)-60
                continue
            end
            %select the last 60s
            sel = (length(t)-30):length(t);
            stats_concat.freq{ipart}{irat}{i_unit}(itrial,:)                        = statsovertime{irat}{ipart}{baseline_index}.freq{i_unit}{itrial}(sel);
            stats_concat.cv2{ipart}{irat}{i_unit}(itrial,:)                         = statsovertime{irat}{ipart}{baseline_index}.cv2{i_unit}{itrial}(sel);
            stats_concat.cv{ipart}{irat}{i_unit}(itrial,:)                          = statsovertime{irat}{ipart}{baseline_index}.cv{i_unit}{itrial}(sel);
            stats_concat.fanofactor{ipart}{irat}{i_unit}(itrial,:)                  = statsovertime{irat}{ipart}{baseline_index}.fanofactor{i_unit}{itrial}(sel);
            stats_concat.burstindex{ipart}{irat}{i_unit}(itrial,:)                  = statsovertime{irat}{ipart}{baseline_index}.burstindex{i_unit}{itrial}(sel);
            stats_concat.amplitude{ipart}{irat}{i_unit}(itrial,:)                   = statsovertime{irat}{ipart}{baseline_index}.amplitude{i_unit}{itrial}(sel);
            stats_concat.cv2_withoutbursts{ipart}{irat}{i_unit}(itrial,:)           = statsovertime_withoutbursts{irat}{ipart}{baseline_index}.cv2{i_unit}{itrial}(sel);
            stats_concat.cv_withoutbursts{ipart}{irat}{i_unit}(itrial,:)            = statsovertime_withoutbursts{irat}{ipart}{baseline_index}.cv{i_unit}{itrial}(sel);
            stats_concat.fanofactor_withoutbursts{ipart}{irat}{i_unit}(itrial,:)    = statsovertime_withoutbursts{irat}{ipart}{baseline_index}.fanofactor{i_unit}{itrial}(sel);
        end
    end
end
stats_concat.time = t(sel) - t(end); %timelock to the end;

for irat = rat_list
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
        %select outliers
        outliers = isoutlier(stats_concat.(method){ipart}{irat}{i_unit});
        %select empty trials
        trials_sel = logical.empty;
        for itrial = 1:size(stats_concat.(method){ipart}{irat}{i_unit},1)
            trials_sel(itrial,1) = sum(stats_concat.freq{ipart}{irat}{i_unit}(itrial,:))>0;
        end
        
        %remove those
        for method = ["cv2_withoutbursts","freq","cv2","burstindex","amplitude", "cv","fanofactor","cv_withoutbursts", "fanofactor_withoutbursts"]%"cv2";"cv2_withoutbursts";
            stats_concat.(method){ipart}{irat}{i_unit}(outliers) = NaN;
            stats_concat.(method){ipart}{irat}{i_unit} = stats_concat.(method){ipart}{irat}{i_unit}(trials_sel,:);
        end
    end
end

save(fullfile(config{1}.datasavedir, 'allrats-statsovertime_concat.mat'),'stats_concat', '-v7.3');
return

%% plot stats over time : freq, cv2, cv2 without bursts
load(fullfile(config{1}.datasavedir, 'allrats-statsovertime_concat'),'stats_concat');
rat_list = 1:5;
% stats_concat.(param){ipart}{irat}{i_unit}(trials,values)
% stats_concat.time

non_norm_extr = [];
last_norm_value = [];
for method = ["cv2_withoutbursts","freq"]%,"cv2","burstindex"]%"cv2";"cv2_withoutbursts";
    figure;hold;setfig();
    for irat = rat_list
        for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
                dofill = false;
%                 continue
            end
            %             if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') %sua
            %                 dofill = true;
            %                 %             continue
            %             end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                color = 'b';
                %             continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                color = 'r';
                %             continue
            end
            
            if isempty(stats_concat.(method){ipart}{irat}{i_unit})
                continue
            end
            
            y = nanmean(stats_concat.(method){ipart}{irat}{i_unit},1);
%             y = nanmedian(stats_concat.(method){ipart}{irat}{i_unit},1);
            x = stats_concat.time;
            
            
%             if max(y(x>-60)) <12 %remove one big outlier
                
                non_norm_extr.(method){irat}{i_unit} = [y(1), y(end)];
%                 y = y./y(x==-60); %normalize
%                 last_norm_value.(method){irat}{i_unit} = y(end);
                plot(x,y,'Color', color);
                %             end
%             end
        end
        %print to file
        ylabel(method, 'Interpreter','none');
        ax = axis;
        %     plot([ax(1) ax(2)], [1 1], 'g', 'LineWidth',2);
    end
    % xlim([-60 0]);
    % ylim([0.5 1.8]);
    % ylim([0 3]);
end


%% plot scatter norm and non norm
figure;hold;setfig();
count_pos = 0;
count_neg = 0;
for method = ["cv2_withoutbursts","freq"]
    for irat = rat_list
        for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
                markerfacecolor = 'none';
                continue
            end
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') %sua
                markerfacecolor = 'k';
                %             continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                plottype = '-^k';
                %             continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                plottype = '-ok';
                %             continue
            end
            if length(non_norm_extr.(method){irat}{i_unit}) == 2
                                            plot([-60 0], non_norm_extr.(method){irat}{i_unit},plottype,'MarkerFaceColor',markerfacecolor); %without normalization
%                 plot([-60 0], [1 last_norm_value.(method){irat}{i_unit}],plottype,'MarkerFaceColor',markerfacecolor); %with normalization
                if last_norm_value.(method){irat}{i_unit} >=1
                    count_pos = count_pos+1;
                else
                    count_neg = count_neg+1;
                end
                
            end
        end
    end
end
xlim([-70 10]);


plot([-70 10], [1 1], '--r');
set(gca,'YScale','log');

% test apparié pour cv2
cv2_withoutbursts_begin_mua = [];
cv2_withoutbursts_end_mua   = [];
cv2_withoutbursts_begin_sua = [];
cv2_withoutbursts_end_sua   = [];
for irat = rat_list
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            continue
        end
        if isempty(non_norm_extr.cv2_withoutbursts{irat}{i_unit})
            continue
        end
        cv2_withoutbursts_begin_mua(end+1) = non_norm_extr.cv2_withoutbursts{irat}{i_unit}(1);
        cv2_withoutbursts_end_mua(end+1) = non_norm_extr.cv2_withoutbursts{irat}{i_unit}(2);
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
            continue
        end
        cv2_withoutbursts_begin_sua(end+1) = non_norm_extr.cv2_withoutbursts{irat}{i_unit}(1);
        cv2_withoutbursts_end_sua(end+1)   = non_norm_extr.cv2_withoutbursts{irat}{i_unit}(2);
    end
end
p_mua_cv2 = signrank(cv2_withoutbursts_begin_mua,cv2_withoutbursts_end_mua);
p_sua_cv2 = signrank(cv2_withoutbursts_begin_sua,cv2_withoutbursts_end_sua);

% test apparié pour freq
freq_begin_mua = [];
freq_end_mua   = [];
freq_begin_sua = [];
freq_end_sua   = [];
for irat = rat_list
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            continue
        end
        if isempty(non_norm_extr.freq{irat}{i_unit})
            continue
        end
        freq_begin_mua(end+1) = non_norm_extr.freq{irat}{i_unit}(1);
        freq_end_mua(end+1) = non_norm_extr.freq{irat}{i_unit}(2);
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
            continue
        end
        freq_begin_sua(end+1) = non_norm_extr.freq{irat}{i_unit}(1);
        freq_end_sua(end+1)   = non_norm_extr.freq{irat}{i_unit}(2);
    end
end
% [h_mua,p_mua] = ttest(freq_begin_mua,freq_end_mua);
% [h_mua,p_mua] = ttest2(freq_begin_mua,freq_end_mua); %non apparié
% [h_sua,p_sua] = ttest(freq_begin_sua,freq_end_sua);
p_mua_freq = signrank(freq_begin_mua,freq_end_mua);
p_sua_freq = signrank(freq_begin_sua,freq_end_sua);
% ranksum
figure;hold;
plot(ones(size(freq_begin_mua))*2+rand, freq_end_mua,'o');
end