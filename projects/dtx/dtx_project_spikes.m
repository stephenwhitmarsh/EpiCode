function dtx_project_spikes(slurm_task_id)

%vérifier rejet artefacts
%vérifier autant de trials que de trials non retirés
%vérifier que es images sont identiques

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


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx


config = dtx_setparams_probe_spikes([]);

%% prepare data for spyking-circus analysis
% for irat = slurm_task_id
%
%     [MuseStruct]                     = readMuseMarkers(config{irat}, false);
% 
%     %remove seizures for non-control experiments
%     if strcmp(config{irat}.type, 'dtx')
%         [MuseStruct]                 = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);
%         %save deadfile without removing seizures
%         writeSpykingCircusDeadFile(config{irat},MuseStruct);
%         %remove all seizures
%         [MuseStruct]                 = addMuseBAD(config{irat},MuseStruct);
%     end
%
%     %create multifile and dead file, and slurm file
%     writeSpykingCircus(config{irat}, MuseStruct, true, true);
%     %create param and prb file
%     writeSpykingCircusParameters(config{irat})
%
%     return
%
% end

%FIXME : dtx_remove_wrong_seizure : output markers to remove, clean it with
%editMuseMarkers


%% analyse spyking circus output
% separation of 'dtx' experiments from 'ctrl' experiments for some
% analysis. Done by filtering with cfg.type

for irat = slurm_task_id
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TEMPORARY_REMOVEME
    %config{irat}.imagesavedir = fullfile(config{irat}.imagesavedir, 'larger_statstime_window');
    %config{irat}.datasavedir = fullfile(config{irat}.datasavedir, 'test_new_orga');
    if ~isfolder(config{irat}.imagesavedir), mkdir (config{irat}.imagesavedir); end 
    if ~isfolder(config{irat}.datasavedir), mkdir (config{irat}.datasavedir); end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %align markers and remove wrong seizures
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false); %align to the begin of the bigger channel

    if strcmp(config{irat}.type, 'dtx') 
        MuseStruct = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true); 
    end
       
    %read LFP data
    dat_LFP                         = readLFP(config{irat}, MuseStruct, false);
    [config{irat},dat_LFP]          = dtx_correctDTX2name(config{irat},dat_LFP); %correct an error in channel name during acquisition, for rat DTX2
    
    %read spike data
    SpikeRaw                        = readSpikeRaw_Phy(config{irat},false);

    if strcmp(config{irat}.type, 'dtx')
        %make trials based on Muse Markers
        SpikeTrials                     = readSpikeTrials_MuseMarkers(config{irat}, MuseStruct,SpikeRaw, false);
    elseif strcmp(config{irat}.type, 'ctrl')
        %make trials of continuous length on all the data
        SpikeTrials                     = readSpikeTrials_continuous(config{irat},SpikeRaw, false);
    end

    clear SpikeRaw
    
    % remove artefacted trials    
    % Plot the trials to see the artefacted Slow Waves (no LFP loaded for ctrl experiments)
    [dat_LFP, ~]                      = removetrials_MuseMarkers(config{irat}, dat_LFP, MuseStruct,true);
    [SpikeTrials, ~]                  = removetrials_MuseMarkers(config{irat}, SpikeTrials, MuseStruct,true);

    %read spike waveforms
    SpikeWaveforms                  = readSpikeWaveforms(config{irat}, SpikeTrials, false);
    
    %create a separated config to avoid useless increase of memory use, if loop over patients
    cfgtemp                 = config{irat};
    cfgtemp.dataLFP         = dat_LFP;          %clear dat_LFP
    cfgtemp.SpikeTrials     = SpikeTrials;      %clear SpikeTrials
    cfgtemp.SpikeWaveforms  = SpikeWaveforms;   %clear SpikeWaveforms
    spikeratestats_Events_Baseline(cfgtemp,true); %no output. Load later for analysis over rats.
    %FIXME voir les try/end dans plot_morpho
    
    %get group info (put by hand)
    if ispc
        unit_table = readtable('Z:\analyses\lgi1\DTX-PROBE\classification_units.xlsx');
    elseif isunix
        unit_table = readtable('/network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/classification_units.xlsx');
    end
    unit_table = table2struct(unit_table);
    rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
    rat_table = unit_table(rat_idx,:);
    for ipart = 1:size(SpikeTrials,2)
        for markername = string(fieldnames(SpikeTrials{ipart}))'
            for i_unit = 1:size(SpikeTrials{ipart}.(markername).label, 2)
                unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), SpikeTrials{ipart}.(markername).label{i_unit});
                % find each element of the unit
                if sum(unit_idx) == 1
                    SpikeTrials{ipart}.(markername).cluster_group{i_unit}             = rat_table(unit_idx).group;
                else
                    SpikeTrials{ipart}.(markername).cluster_group{i_unit}             = 'noise';
                end
                if isempty(rat_table(unit_idx).group)
                    SpikeTrials{ipart}.(markername).cluster_group{i_unit} = 'noise';
                end
            end
        end
    end

    %plot a summary of all spike data
    plot_spike_quality(config{irat},SpikeTrials, SpikeWaveforms,true);
    cfgtemp                     = config{irat};
    cfgtemp.spikequal.suffix    = '_nostd';
    cfgtemp.spikequal.plot_std  = 'no';
    plot_spike_quality(cfgtemp,SpikeTrials, SpikeWaveforms,true);
    
end
return

%% Load precomputed data
%FIXME : à corriger avec la nouvelle organisation
ipart = 1;
rat_list = 1:7;
for irat = rat_list
%      config{irat}.datasavedir = fullfile(config{irat}.datasavedir, 'test_alignXCorr');
    stats{irat} = spikeratestats_Events_Baseline(config{irat},true);
end

%stats_concat : see dtx_cluster_statsovertime : stats des 60 dernières
%secondes pré-crise, seulement pour les trials où la période interictale
%est de plus de 60 secondes
load(fullfile(config{1}.datasavedir, 'allrats-statsovertime_concat.mat'),'stats_concat');


unit_table = readtable('Z:\analyses\lgi1\DTX-PROBE\classification_units.xlsx');
unit_table = table2struct(unit_table);
for irat = rat_list
    rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
    rat_table = unit_table(rat_idx,:);
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit});
        % find each element of the unit
        if sum(unit_idx) == 1
            stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}             = rat_table(unit_idx).group;
            if strcmp(config{irat}.type, 'dtx'), stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).code_slowwave{i_unit} = rat_table(unit_idx).code_slowwave_spikerate; end
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_spikerate{i_unit}    = rat_table(unit_idx).code_interictal_spikerate; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_cv2{i_unit}          = rat_table(unit_idx).code_interictal_cv2; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit}           = rat_table(unit_idx).Emax;
        end
    end
end

%replace empty cells by nans
for irat = rat_list
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit})
            stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit} = NaN;
        end
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit} = 'noise';
        end
    end
end

%% IN PN classification
figure;hold;
check_empty = [];
for irat = rat_list
    ipart = 1;
    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype = cell(1,size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2));
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                x = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.halfwidth{i_unit}];
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.peaktrough{i_unit}];
                z = [stats{irat}{ipart}.(config{irat}.spike.baselinename).template.troughpeak{i_unit}];
                %best clustering with halfwidth and throughpeak
                
                %values to separate in and pn determined empirically
%                 if x>2.25*10^-4 && z>6*10^-4 %PN
                if x>3*10^-4 || (x>2*10^-4&& z>4.5*10^-4 ) || z>6*10^-4%PN
                    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit} = 'pn';
                    if strcmp(config{irat}.type, 'dtx'), plottype = '^b';end
                    if strcmp(config{irat}.type, 'ctrl'), plottype = 'ob';end
                else %IN
                    stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit} = 'in';
                    if strcmp(config{irat}.type, 'dtx'), plottype = '^r';end
                    if strcmp(config{irat}.type, 'ctrl'), plottype = 'or';end
                end
                
                if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
                    scatter(x,z,plottype, 'filled');
                else
                     scatter(x,z,plottype);
                end
            end
        end
    end
end
ax = axis;
leg1 = scatter(1,1,'^k'); %mua
leg2 = scatter(1,1,'^k','filled'); %sua
leg3 = scatter(1,1,'^b', 'filled'); %pn
leg4 = scatter(1,1,'^r', 'filled'); %in
leg5 = scatter(1,1,'^k','filled'); %dtx
leg6 = scatter(1,1,'ok','filled'); %ctrl
legend([leg1 leg2 leg3 leg4 leg5 leg6],'MUA','SUA', 'PN','IN','DTX','CTRL');
% xlim([ax(1) ax(2)]); ylim([ax(3) ax(4)]);
xlim([0.1 0.6].*10^-3); ylim([0.25 1.55].*10^-3); %same as Bartho 2004
xticklabels(xticks.*1000); %convert to ms
yticklabels(yticks.*1000); %convert to ms
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('halfwidth(ms)');
ylabel('peak-trough (ms)');

%% norm morpho selon type
figure;hold;
check_empty = [];
for irat = rat_list
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
%             if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                    clear temp
                    for i=1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit},1)
                        temp(i,:) = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit}{i,:}; 
                    end
                    plot(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.time{1}{1}, nanmean(temp,1)/max(nanmean(temp,1)), 'b');
                elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                    clear temp
                    for i=1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit},1)
                        temp(i,:) = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.values{i_unit}{i,:}; 
                    end
                    x = stats{irat}{ipart}.(config{irat}.spike.baselinename).template.time{1}{1};
%                     plot(x+x(end)+0.0005, nanmean(temp,1)/max(nanmean(temp,1)), 'r');
                    plot(x, nanmean(temp,1)/max(nanmean(temp,1)), 'r');
                end

            end
        end
    end
end

%% FREQ
%comparer DTX toute la période intercritique, DTX les 10 dernières
%secondes, et ctrl
figure;hold;
% ilabel = 3; %interictal
markername = 'Interictal';
y = cell(1,3);
for irat = rat_list
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                plottype = 'ob';
            elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                plottype = 'or';
            end
            if strcmp(config{irat}.type, 'dtx'),   idx = 1;
            elseif strcmp(config{irat}.type, 'ctrl'), idx=3;
            end
            
            %plot data of all periods
            y{idx}(end+1) = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meanfreq{i_unit}];
            x = idx+(rand-0.5)*0.2;
            
            scatter(x,y{idx}(end),plottype, 'filled');
            
%             plot data of the last 10 seconds before the seizure
            if strcmp(config{irat}.type, 'dtx')
                t = stats_concat.time;
                idx = 2;
                y{idx}(end+1) = nanmean(nanmean(stats_concat.freq{ipart}{irat}{i_unit}(:,t>=-1)));
                x = idx+(rand-0.5)*0.2;
                scatter(x,y{idx}(end),plottype, 'filled');
            end
        end
    end
end
xlim([0 4]);

%stats
pval_freq = ranksum(y{1}, y{3});

sem = @(data) std(data)/sqrt(length(data));
errorbar([1 2 3], cellfun(@nanmean,y),cellfun(@sem,y));


%% CV2
figure;hold;
y = cell(1,2);
check_empty = [];
for irat = rat_list
    ipart = 1;
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
            check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
        else
            if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                    plottype = 'ob';
                elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                    plottype = 'or';
                end
                
                if strcmp (config{irat}.type, 'dtx'),   idx = 1; else, idx=2; end
                x = idx+(rand-0.5)*0.2;
                y{idx}(end+1) = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meancv2{i_unit}];
                
                scatter(x,y{idx}(end),plottype, 'filled');
            end
        end
    end
end
xlim([0 3]);
ylim([0 2]);

%stats
pval_cv2 = ranksum(y{1}, y{2});

%% plot sdf of each neuron during slow wave
% imagesc(log10(stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).sdf.avg));
% xlim([80 120])

figure;hold;
emax_list = nan;
sdf_avg = nan(1,length(stats{2}{ipart}.(config{1}.spike.eventsname{1}).sdfavg.avg));
sdf_time = stats{2}{ipart}.(config{1}.spike.eventsname{1}).sdfavg.time;
celltype = nan;
group = nan;
for irat = 1:5
    %for only sua
%     idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'mua') & ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    %for sua and mua
    idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    emax_list   = [emax_list, stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{idx}];
    sdf_avg     = vertcat(sdf_avg, stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).sdfavg.avg(idx,:));
    celltype    = [celltype, stats{irat}{ipart}.(config{1}.spike.baselinename).celltype(idx)];
    group       = [group, stats{irat}{ipart}.(config{1}.spike.baselinename).group(idx)];
end
[emax_list_sorted, deep_idx] = sort(emax_list, 'descend');

deep_idx            = deep_idx(~isnan(emax_list_sorted));
emax_list_sorted    = emax_list_sorted(~isnan(emax_list_sorted));
sdf_avg             = sdf_avg(deep_idx,:);
celltype            = celltype(deep_idx);
group               = group(deep_idx);

for i=1:size(sdf_avg,1)
%     sdf_avg_filt(i,:) = imgaussfilt(sdf_avg(i,:),1); %smooth the spike density function values
    sdf_blcorrect(i,:) = (sdf_avg(i,:)) - nanmean((sdf_avg(i,sdf_time>-1.9 & sdf_time<-1.4))); %correct baseline
    sdf_norm(i,:) = normalize(sdf_avg(i,:),'zscore');
end

% plot pn and in separately
figure;
for i_celltype = ["pn", "in"]
    
    if strcmp(i_celltype, "pn"), subplot(2,1,1);hold; else, subplot(2,1,2);hold; end
    
%     imagesc(log10(abs(sdf_norm(strcmp(celltype,i_celltype),:))));
%     imagesc(log10(abs(sdf_avg(strcmp(celltype,i_celltype),:))));
    imagesc(log10(abs(sdf_blcorrect(strcmp(celltype,i_celltype),:))));
    
    c = caxis;
    caxis([0 c(2)]);
    clb = colorbar;
    axis tight
    xticklabels(sdf_time(xticks+1));
    func_temp = @(v) sprintf('10\\^%s',v);
%       func_temp = @(v) eval(sprintf('round(10^%s)',v));
%     clb.TickLabels = cellfun(func_temp,clb.TickLabels,'UniformOutput',false);
   
    % plot baseline patch and zero line
    x1 = find(sdf_time > -1.9, 1, 'first');
    x2 = find(sdf_time < -1.4, 1, 'last');
    x = [x1 x2 x2 x1];
    ax = axis;
    y = [ax(3) ax(3) ax(4) ax(4)];
    patch('XData',x,'YData',y,'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);
    x = find(sdf_time == 0);
    plot([x x], [ax(3:4)], '--r');
    
    %plot electrode infos
    C = linspecer(sum(~isnan(unique(emax_list_sorted(strcmp(celltype,i_celltype))))),'qualitative');
    emax_pn = emax_list_sorted(strcmp(celltype,i_celltype));
    group_temp = group(strcmp(celltype,i_celltype));
    for ineuron = 1:size(celltype(strcmp(celltype,i_celltype)),2)
        chan_idx=find(unique(emax_pn)==emax_pn(ineuron));
        scatter(size(sdf_avg,2)+2,ineuron,'s','filled','MarkerEdgeColor', C(chan_idx,:), 'MarkerFaceColor', C(chan_idx,:));
        if ~contains(group_temp(ineuron), 'mua')%noise already removed
            scatter(size(sdf_avg,2)-1,ineuron,'<r','filled');
        end
    end
    set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
    xlabel('Time from begin of SlowWave (s)');
end

%% plot of cv2 with bursts, without bursts and freq over time : avg for each neuron
%Over time. OLD.
%normalized at t=-60
for method = ["cv2_withoutbursts","freq"]%"cv2";"cv2_withoutbursts";
    figure;hold;
    for irat = rat_list
        if strcmp(config{irat}.type, 'dtx')
            ipart = 1;
            for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
                if isempty(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit})
                    check_empty = [check_empty, [num2str(irat), ';',stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit},';']]
                else
                    if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
%                         if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
                        
                        x = stats{irat}{ipart}.(config{irat}.spike.baselinename).stats_over_time.(method){i_unit}.time;
                        y = stats{irat}{ipart}.(config{irat}.spike.baselinename).stats_over_time.(method){i_unit}.avg;
                        
                        if strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
                            color = 'b';
                        elseif strcmp(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
                            color = 'r';
                        end
                        y = y./y(x==-60); %normalize
                        if max(y(x>-60)) <10 %avoir outliers
                            plot(x,y,'Color', color);
                        end
                    end
                end
            end
        end
    end
end
xlim([-60 0]);
ylim([0.5 1.8]);
ylim([0 3]);

%% compter code SlowWave
code_slowwave_mua = [];
code_slowwave_sua = [];
freq_slowwave_mua = [];
freq_slowwave_sua = [];
celltype_slowwave_sua = {};
celltype_slowwave_mua = {};
maxchan_slowwave_mua = [];
maxchan_slowwave_sua = [];
maxfreq_slowwave_mua = [];
maxfreq_slowwave_sua = [];
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            continue
        end
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
            code_slowwave_mua(end+1) = stats{irat}{ipart}.SlowWave.code_slowwave{i_unit};
            freq_slowwave_mua(end+1) = stats{irat}{ipart}.Interictal.discharge.meanfreq{i_unit};
            celltype_slowwave_mua{end+1} = stats{irat}{ipart}.Interictal.celltype{i_unit};
            maxchan_slowwave_mua(end+1) = stats{irat}{ipart}.Interictal.maxchan{i_unit};
            maxfreq_slowwave_mua(end+1) = stats{irat}{ipart}.SlowWave.maxfreq.max{i_unit};
        elseif contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'sua')
%             code_slowwave_mua(end+1) = stats{irat}{ipart}.SlowWave.code_slowwave{i_unit};
%             freq_slowwave_mua(end+1) = stats{irat}{ipart}.Interictal.discharge.meanfreq{i_unit};
            code_slowwave_sua(end+1) = stats{irat}{ipart}.SlowWave.code_slowwave{i_unit};
            freq_slowwave_sua(end+1) = stats{irat}{ipart}.Interictal.discharge.meanfreq{i_unit};
            maxfreq_slowwave_sua(end+1) = stats{irat}{ipart}.SlowWave.maxfreq.max{i_unit};
            celltype_slowwave_sua{end+1} = stats{irat}{ipart}.Interictal.celltype{i_unit};
            maxchan_slowwave_sua(end+1) = stats{irat}{ipart}.Interictal.maxchan{i_unit};
        end
    end
end
bar(histcounts(code_slowwave_mua, 1:6));
bar(histcounts(code_slowwave_sua, 1:6));

%% groupe SW en fonction de la fréquence pré SW, du type cellulaire, de la localisation, du max de décharge pendant SW
figure;hold;
for i_unit = 1:size(code_slowwave_mua,2)
    if strcmp(celltype_slowwave_mua{i_unit}, 'pn')
        plottype = '^k';
        code_celltype(i_unit) = 1;
    else
        code_celltype(i_unit) = 2;
        plottype = 'ok';
    end
%     scatter(code_slowwave_mua(i_unit)+rand*0.3, freq_slowwave_mua(i_unit),plottype);
%     scatter(code_slowwave_mua(i_unit)+rand*0.3, freq_slowwave_mua(i_unit),plottype);
    scatter(code_slowwave_mua(i_unit)+rand*0.3, maxfreq_slowwave_mua(i_unit),plottype);
end

%count pn and in for each type of sw
nb_pn = sum(code_celltype==1);
nb_in = sum(code_celltype==2);
y=[];
for i_sw = 1:5
    y(end+1:end+2) = histcounts(code_celltype(code_slowwave_mua==i_sw), 1:3); %une valeur sur 2 : pn, in, pn, in etc
end
%normalize pn and in compared to there proportion in the whole recorded cells
for i = 1:2:length(y)
    y(i) = y(i)/nb_pn;
end
for i=2:2:length(y)
    y(i) = y(i)/nb_in;
end
bar(y);
%plot the ratio in/pn
ratio_in_pn = [];
for i=1:2:length(y)
    ratio_in_pn(end+1) = y(i+1)/y(i); %nb in divisé par nb pn
end
bar(ratio_in_pn);
ratio_in_pn(end+1) = nb_in/nb_pn;

for i_unit = 1:size(code_slowwave_sua,2)
    if strcmp(celltype_slowwave_sua{i_unit}, 'pn')
        plottype = '^k';
    else
        plottype = 'ok';
    end
%     scatter(code_slowwave_sua(i_unit)+rand*0.3, freq_slowwave_sua(i_unit), plottype,'filled');
%     scatter(code_slowwave_sua(i_unit)+rand*0.3, maxchan_slowwave_sua(i_unit), plottype,'filled');
scatter(code_slowwave_mua(i_unit)+rand*0.3, maxfreq_slowwave_mua(i_unit),plottype, 'filled');
end
xlim([0.5 5.5]);
setfig()

%% stats over time
% voir dtx_cluster_statsovertime

%% distrib of rpv
% ilabel = 3; %interictal
markername = 'Interictal';
for irat = 1:length(config)
    rpv{irat} = plot_spike_quality(config{irat},[], [],false);
end

figure;hold;setfig();
for irat = 1:5%length(config)
    for i_unit = 1:size(stats{irat}{ipart}.Interictal.label,2)
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'noise')
            continue
        end
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
            dofill = false;
        end
        if ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua') %sua
            dofill = true;
        end
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'pn')
            plottype = '^k';
        end
        if contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).celltype{i_unit}, 'in')
            plottype = 'ok';
        end
        
        x = rand;
        y = rpv{irat}{ipart}.(markername)(i_unit);
        
%         if y>10
%             toremove{irat}.label{i_unit} = stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit};
%             continue
%         end
%         
%         if y>1 && ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}, 'mua')
%             toremove{irat}.label{i_unit} = stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit};
%             continue
%         end
        
        if dofill
            scatter(x,y,plottype,'filled');
        else
            scatter(x,y,plottype);
        end
    end
end
xlim([-0.5 1.5]);
