function dtx_project_spikes(slurm_task_id)

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


%% analyse spyking circus output
% separation of 'dtx' experiments from 'ctrl' experiments for some
% analysis. Done by filtering with cfg.type

for irat = slurm_task_id
%     %%%%%%%%%%%%%%%% REMOVE
%     [MuseStruct]                    = readMuseMarkers(config{irat}, false);
%     SpikeRaw                        = readSpikeRaw_Phy(config{irat},false);
%     SpikeTrials                     = readSpikeTrials_MuseMarkers(config{irat}, MuseStruct,SpikeRaw, true);
%     return
%     %%%%%%%%%%%%%%%%%%%%%%%%
   
    %TEMPORARY_REMOVEME
%     config{irat}.imagesavedir = fullfile(config{irat}.imagesavedir, 'test_alignXCorr');
%     if ~isfolder(config{irat}.imagesavedir), mkdir (config{irat}.imagesavedir); end 

    %align markers and remove wrong seizures
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);
%     MuseStruct = alignMuseMarkersXcorr(config{irat}, MuseStruct, true);
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
    cfgtemp                         = [];
    cfgtemp                         = config{irat}; %info of where to save images
    cfgtemp.remove.plotdata         = 'yes';
    cfgtemp.remove.electrodetoplot  = ft_getopt(config{irat}.align,'channel', []);
    [dat_LFP, ~]                    = removetrials_MuseMarkers(cfgtemp, dat_LFP, MuseStruct);
    [SpikeTrials, ~]                = removetrials_MuseMarkers(cfgtemp, SpikeTrials, MuseStruct);

    %read spike waveforms
    SpikeWaveforms                  = readSpikeWaveforms_trials(config{irat}, SpikeTrials, true);
    
    %create a separated config to avoid useless increase of memory use, if loop over patients
    cfgtemp                 = [];
    cfgtemp                 = config{irat};
    cfgtemp.dataLFP         = dat_LFP;          clear dat_LFP
    cfgtemp.SpikeTrials     = SpikeTrials;      clear SpikeTrials
    cfgtemp.SpikeWaveforms  = SpikeWaveforms;   clear SpikeWaveforms
    spikeratestats_Events_Baseline(cfgtemp,true); %no output. Load later for analysis over rats.
    %FIXME voir les try/end dans plot_morpho
              
end %irat
return


%% Load precomputed data
ipart = 1;
for irat = 1:7
%      config{irat}.datasavedir = fullfile(config{irat}.datasavedir, 'test_alignXCorr');
    stats{irat} = spikeratestats_Events_Baseline(config{irat},false);
end

%do it after filtering and keeping only sua
unit_table = readtable('Z:\analyses\lgi1\DTX-PROBE\classification_units.xlsx');
unit_table = table2struct(unit_table);
for irat = 1:7
    rat_idx = strcmp({unit_table.ratID}, config{irat}.prefix(1:end-1))';
    rat_table = unit_table(rat_idx,:);
    for i_unit = 1:size(stats{irat}{ipart}.(config{irat}.spike.baselinename).spikewaveform.halfwidth, 2)
        unit_idx = strcmp(split(sprintf('cluster_%d,', rat_table.clusterID), ','), stats{irat}{ipart}.(config{irat}.spike.baselinename).label{i_unit});
        % find each element of the unit
        if sum(unit_idx) == 1
            stats{irat}{ipart}.(config{irat}.spike.baselinename).group{i_unit}              = rat_table(unit_idx).group;
            if strcmp(config{irat}.type, 'dtx'), stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).code_slowwave{i_unit} = rat_table(unit_idx).code_slowwave_spikerate; end
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_spikerate{i_unit}    = rat_table(unit_idx).code_interictal_spikerate; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).code_cv2{i_unit}          = rat_table(unit_idx).code_interictal_cv2; 
            stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{i_unit}            = rat_table(unit_idx).Emax;
        end
    end
end

%replace empty cells by nans
for irat = 1:7
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
for irat = 1:7
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
xlim([ax(1) ax(2)]); ylim([ax(3) ax(4)]);
xticklabels(xticks.*1000); %convert to ms
yticklabels(yticks.*1000); %convert to ms
set(gca, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 10);
xlabel('halfwidth(ms)');
ylabel('peak-trough (ms)');

%% norm morpho selon type
figure;hold;
check_empty = [];
for irat = 1:7
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
figure;hold;
check_empty = [];
for irat = 1:7
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
                if strcmp(config{irat}.type, 'dtx'),   x = 1; 
                elseif strcmp(config{irat}.type, 'ctrl'), x=2; 
                end
                
                x = x+(rand-0.5)*0.2;
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meanfreq{i_unit}];
                
                scatter(x,y,plottype, 'filled');
            end
        end
    end
end
xlim([0 3]);

%% CV2
figure;hold;
check_empty = [];
for irat = 1:7
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
                
                if strcmp (config{irat}.type, 'dtx'),   x = 1; else, x=2; end
                x = x+(rand-0.5)*0.2;
                y = [stats{irat}{ipart}.(config{irat}.spike.baselinename).discharge.meancv2{i_unit}];
                
                scatter(x,y,plottype, 'filled');
            end
        end
    end
end
xlim([0 3]);
ylim([0 2]);

%% plot sdf of each neuron during slow wave
% imagesc(log10(stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).sdf.avg));
% xlim([80 120])

figure;hold;
emax_list = nan;
sdf_avg = nan(1,length(stats{1}{ipart}.(config{1}.spike.eventsname{1}).sdf.avg));
sdf_time = stats{1}{ipart}.(config{1}.spike.eventsname{1}).sdf.time;
celltype = nan;
group = nan;
for irat = 1:5
    %for only sua
%     idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'mua') & ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    %for sua and mua
    idx = ~contains(stats{irat}{ipart}.(config{irat}.spike.baselinename).group, 'noise');
    emax_list   = [emax_list, stats{irat}{ipart}.(config{irat}.spike.baselinename).maxchan{idx}];
    sdf_avg     = vertcat(sdf_avg, stats{irat}{ipart}.(config{irat}.spike.eventsname{1}).sdf.avg(idx,:));
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
    sdf_avg_filt(i,:) = imgaussfilt(sdf_avg(i,:),1); %smooth the spike density function values
    sdf_filt_blcorrect(i,:) = (sdf_avg_filt(i,:)) - nanmean((sdf_avg_filt(i,1:20))); %correct baseline
end

% plot pn and in separately
figure;
for i_celltype = ["pn", "in"]
    
    if strcmp(i_celltype, "pn"), subplot(2,1,1);hold; else, subplot(2,1,2);hold; end
    imagesc(log10(abs(sdf_filt_blcorrect(strcmp(celltype,i_celltype),:))));
    c = caxis;
    caxis([0 1.5]);
    clb = colorbar;
    axis tight
    xticklabels(sdf_time(xticks+1));
    func_temp = @(v) sprintf('10\\^%s',v);
%       func_temp = @(v) eval(sprintf('round(10^%s)',v));
    clb.TickLabels = cellfun(func_temp,clb.TickLabels,'UniformOutput',false);
    % plot baseline patch and zero line
    x = [0 20 20 0];
    ax = axis;
    y = [ax(3) ax(3) ax(4) ax(4)];
    patch('XData',x,'YData',y,'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);
    plot([100 100], [ax(3:4)], '--r');
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

%normalized at t=-60
for method = ["cv2_withoutbursts","freq"]%"cv2";"cv2_withoutbursts";
    figure;hold;
    for irat = 1:7
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

end
%
