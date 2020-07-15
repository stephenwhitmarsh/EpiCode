function dtx_cluster_statsovertime
%run it on the cluster because stats over time is too big
% statsovertime{ipart}{ilabel}.(param){i_unit}{i_trial}
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

save(fullfile(config{1}.datasavedir, 'allrats-statsovertime_concat.mat'),'stats_concat', '-v7.3');
return

%% plot stats over time : freq, cv2, cv2 without bursts
load(fullfile(config{1}.datasavedir, 'allrats-statsovertime_concat'),'stats_concat');
% stats_concat.(param){ipart}{irat}{i_unit}(trials,values)
% stats_concat.time

non_norm_extr = [];
last_norm_value = [];
for method = ["cv2_withoutbursts","freq","cv2","burstindex"]%"cv2";"cv2_withoutbursts";
    figure;hold;setfig();
    for irat = rat_list
        for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).label, 2)
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                continue
            end
            if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
                dofill = false;
                                            continue
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
            
            %remove outliers
            outliers = isoutlier(stats_concat.(method){ipart}{irat}{i_unit});
            stats_concat.(method){ipart}{irat}{i_unit}(outliers) = NaN;
            
            y = nanmean(stats_concat.(method){ipart}{irat}{i_unit});
            x = stats_concat.time;
            
            y = y./y(x==-60); %normalize
            if max(y(x>-60)) <12 %remove one big outlier    
%                 last_norm_value{irat}{i_unit} = y(end);
                non_norm_extr{irat}{i_unit} = [nanmean(stats_concat.(method){ipart}{irat}{i_unit}(:,1)), nanmean(stats_concat.(method){ipart}{irat}{i_unit}(:,end))];
                plot(x,y,'Color', color);
%                 plot(x,nanmean(stats_concat.(method){ipart}{irat}{i_unit}),'Color', color);
%             end
        end
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
            if length(non_norm_extr{irat}{i_unit}) == 2
%                 plot([-60 0], non_norm_extr{irat}{i_unit},plottype,'MarkerFaceColor',markerfacecolor); %without normalization
                plot([-60 0], [1 last_norm_value{irat}{i_unit}],plottype,'MarkerFaceColor',markerfacecolor); %with normalization
                if last_norm_value{irat}{i_unit} >=1
                    count_pos = count_pos+1;
                else
                    count_neg = count_neg+1;
                end
                
            end
        end
end
xlim([-70 10]);


plot([-70 10], [1 1], '--r');
set(gca,'YScale','log');

end