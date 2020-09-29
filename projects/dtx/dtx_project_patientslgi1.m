function dtx_project_patientslgi1(slurm_task_id)


%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end


ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = dtx_setparams_patients_lgi1;

for ipatient = slurm_task_id%1:size(config,2)
    MuseStruct = readMuseMarkers(config{ipatient},true);
    MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,true); 
    LFP = readLFP(config{ipatient},MuseStruct,true);
return
end

% config_origin = config;

% -	EMG contraction : duration of contraction :                         OK
% -	EEG-EMG delay :                                                     OK
% - frequency and regularity of seizures : 
%   * convertir temps des donn�es coup�es, grace au fichier txt         OK
%   * distrib des inter-seizure interval et des cv2                     OK
% Voir sur quel channel, et quels features, extraire
% -	Slow wave : 
%   * morphology (half width C4, half width better channel, autre ?) 
%   * topography : 1 topoplot par patient et par side (L, R)
%   * propagation : xcorr
% - traces brutes repr�sentatives : 1 crise pat 008, + toutes les crises
% superpos�es. Vid�o peccoud avec le m�me eeg
% Alignement pour topoplot
% - plot de chaque crise non artefact�e raw EEG + topoplot : avec dir et temps pour la
% retrouver. Pour choisir les mieux pour les measures. Plots flipped et
% non-flipped.  
% - pr�voir une vid�o de l'�v�nement
% comparer avg ol CSD et LFP



%% eeg emg delay, and emg duration
pat_list = 1:size(config,2);
for ipatient = pat_list
    
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, false);
    %take the last part
    ipart = size(MuseStruct,2); 
    
    %% compute eeg-emg delays and emg duration
    slowwave_begin_R      = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_R_begin');
    emg_begin_R           = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_R_EMG__START__');
    emg_end_R             = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_R_EMG__END__');
    slowwave_begin_L      = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_L_begin');
    emg_begin_L           = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_L_EMG__START__');
    emg_end_L             = concatenateMuseMarker(config{ipatient},MuseStruct,ipart,'SlowWave_L_EMG__END__');
    
    %put right and left seizures together
    slowwave_begin  = [slowwave_begin_R.synctime,slowwave_begin_L.synctime];
    emg_begin       = [emg_begin_R.synctime,emg_begin_L.synctime];
    emg_end         = [emg_end_R.synctime,emg_end_L.synctime];
    
    if isempty(slowwave_begin)
        delays{ipatient} = NaN;
        emg_duration{ipatient} = NaN;
        continue
    end
    
    delays{ipatient}        = emg_begin - slowwave_begin;
    emg_duration{ipatient}  = emg_end - emg_begin;
    eeg_dir{ipatient}           = [slowwave_begin_R.dir,slowwave_begin_L.dir];
    %     [bins, edges]       = histcounts(delays{ipatient},'BinWidth',0.01);
    %     bins_centers        = (edges(1:end-1)+edges(2:end))/2; %limiteinf + limitesup / 2
    %     bar(bins_centers,bins);
    %
end

%eeg emg delay
figure;hold
for ipatient = pat_list
    scatter(rand(size(delays{ipatient}))*0.2+ipatient-0.1, delays{ipatient}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    meandelay(ipatient) = nanmean(delays{ipatient});
    stddelay(ipatient)  = nanstd(delays{ipatient});
    errorbar(ipatient, meandelay(ipatient), stddelay(ipatient),'--rx');
    %     errorbar(ipatient, nanmean(delays{ipatient}), nanmean(delays{ipatient})/sqrt(size(delays{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 pat_list(end)+1]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('eeg-emg delay (s)');

figure;hold
scatter(rand(size(meandelay))*0.2, meandelay,'sk','MarkerFaceColor','k');
errorbar(0.1,nanmean(meandelay),nanstd(meandelay),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
xlim([-0.5 0.7]);
xticks([]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('mean eeg-emg delay, per patient (s)');
ax = axis; ylim([0 ax(4)]);

%emg duration
figure;hold
for ipatient = pat_list
    meanemgduration(ipatient) = nanmean(emg_duration{ipatient});
    stdemgduration(ipatient)  = nanstd(emg_duration{ipatient});
    scatter(rand(size(emg_duration{ipatient}))*0.2+ipatient-0.1, emg_duration{ipatient}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, meanemgduration(ipatient), stdemgduration(ipatient),'--rx');
    %     errorbar(ipatient, nanmean(emg_duration{ipatient}), nanmean(emg_duration{ipatient})/sqrt(size(emg_duration{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 pat_list(end)+1]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('emg duration (s)');

figure;hold
scatter(rand(size(meanemgduration))*0.2, meanemgduration,'sk','MarkerFaceColor','k');
errorbar(0.1,nanmean(meanemgduration),nanstd(meanemgduration),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
xlim([-0.5 0.7]);
xticks([]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('mean emg duration, per patient (s)');
ax = axis; ylim([0 ax(4)]);

% %check dir data
% figure;
% subplot(2,1,1);hold;
% scatter(eeg_dir{ipatient},emg_duration{ipatient});
% title(sprintf('emg duration patient %d',ipatient));
% subplot(2,1,2);hold;
% scatter(eeg_dir{ipatient},delays{ipatient});
% title(sprintf('eeg-emg delays patient %d',ipatient));
%svg stats en output


%% plot seizure frequency and regularity
%recover real time in segmented data OK
%correct concatenateMuseMarkers : clocktime non chang�, synctime en
%fonction du nb de samples. A ne pas utiliser si les fichiers ne sont pas
%continus.
for ipatient = 1:size(config,2)
    MuseStruct = readMuseMarkers_discontinuousMicromed(config{ipatient},false);
    
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_L';
    SlowWave_L(ipatient) = dtx_stats_seizure_timings(config{ipatient},MuseStruct,1);
    
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_R';
    SlowWave_R(ipatient) = dtx_stats_seizure_timings(config{ipatient},MuseStruct,1);
end

data_temp_L = {SlowWave_L.time_start};
data_temp_R = {SlowWave_R.time_start};

%time between 2 seizures
figure;hold
for ipatient = 1:size(config,2)
    %sort data and remove times between 2 dirs
    [data, indexes] = sort([data_temp_R{ipatient}.clock, data_temp_L{ipatient}.clock]);
    dirtemp = [data_temp_R{ipatient}.dir, data_temp_L{ipatient}.dir];
    dir_list = dirtemp(indexes);
    data = minutes(diff(data));
    data_L_R{ipatient} = data(~logical(diff(dir_list)));
    %plot
    scatter(rand(size(data_L_R{ipatient}))*0.2+ipatient-0.1, data_L_R{ipatient}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, nanmean(data_L_R{ipatient}), nanstd(data_L_R{ipatient}),'--rx');
end
xlim([0 size(config,2)+1]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('time between 2 slow waves (right and left, minutes)');
ax = axis;
ylim([0 ax(4)]);

data_L = {SlowWave_L.timebetween2seizures};
data_R = {SlowWave_R.timebetween2seizures};
%time between 2 seizures right
figure;hold
for ipatient = 1:size(config,2)
    data = minutes([data_R{ipatient}]);
    scatter(rand(size(data))*0.2+ipatient-0.1, data, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, nanmean(data), nanstd(data),'--rx');
end
xlim([0 size(config,2)+1]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('time between 2 slow waves (right only, minutes)');
ax = axis;
ylim([0 ax(4)]);

%time between 2 seizures left
figure;hold
for ipatient = 1:size(config,2)
    data = minutes([data_L{ipatient}]);
    scatter(rand(size(data))*0.2+ipatient-0.1, data, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(ipatient, nanmean(data), nanstd(data),'--rx');
end
xlim([0 size(config,2)+1]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('time between 2 slow waves (left only, minutes)');
ax = axis;
ylim([0 ax(4)]);

%r, l, r+l in the same plot
meandata_R_L = cellfun(@nanmean,data_L_R);
meandata_R = cellfun(@(v) nanmean(minutes(v)),data_R);
meandata_L = cellfun(@(v) nanmean(minutes(v)),data_L);
figure;hold
count_plot = 0;
for idata = {meandata_R_L,meandata_R,meandata_L}
    count_plot = count_plot+1;
    scatter(rand(size(idata{1}))*0.2+count_plot-0.1, idata{1},'sk','MarkerFaceColor','k');
    errorbar(count_plot,nanmean(idata{1}),nanstd(idata{1}),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
end
xlim([0.5 3.5]);
xticks(1:3);
xticklabels({'R + L', 'R', 'L'});
set(gca,'TickDir','out','FontWeight','bold');
ylabel('time between 2 slow waves, mean per patient');
ax = axis; ylim([0 ax(4)]);

%cv2 of right and left, only right, only left
for ipatient = 1:size(config,2)
    data = data_L_R{ipatient};
    data_L_R_cv2(ipatient) = nanmean(cv2(data));
    data = minutes([data_L{ipatient}]);
    data_L_cv2(ipatient) = nanmean(cv2(data));
    data = minutes([data_R{ipatient}]);
    data_R_cv2(ipatient) = nanmean(cv2(data));
end
figure;hold
scatter_position = 0;
for idata = {data_L_R_cv2,data_R_cv2,data_L_cv2}
    scatter_position = scatter_position +1;
    scatter(rand(size(idata{1}))*0.2+scatter_position-0.1, idata{1}, '.', 'MarkerEdgeColor', [0.6 0.6 0.6]);
    errorbar(scatter_position, nanmean(idata{1}), nanstd(idata{1}),'-rx');
end
xlim([0.5 3.5]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('Patients'' seizures cv2');
xticks(1:3);
xticklabels({'L + R', 'R only', 'L only'});
ylim([0 2]);

%line plot cv2
figure;hold
for ipatient = 1:size(config,2)%{data_L_R_cv2,data_L_cv2,data_R_cv2}
    plot([1 2 3], [data_L_R_cv2(ipatient),data_R_cv2(ipatient),data_L_cv2(ipatient)],'-sk','MarkerFaceColor','k');
end
errorbar([1 2 3], nanmean([data_L_R_cv2;data_R_cv2;data_L_cv2],2)',nanstd([data_L_R_cv2;data_R_cv2;data_L_cv2],0,2)','--rx','LineWidth',2);
xlim([0.5 3.5]);
set(gca,'TickDir','out','FontWeight','bold');
ylabel('Patients'' seizures cv2');
xticks(1:3);
xticklabels({'L + R', 'R only', 'L only'});
ylim([0 2]);

%signrank : paired test. ranksum(x,y) : x et y ind�pendants
p_all_vs_l = signrank(data_L_R_cv2,data_L_cv2);
p_all_vs_r = signrank(data_L_R_cv2,data_R_cv2);
p_r_vs_l   = signrank(data_R_cv2,data_L_cv2);
p_all_vs_l_indpt = ranksum(data_L_R_cv2,data_L_cv2);
p_all_vs_r_indpt = ranksum(data_L_R_cv2,data_R_cv2);
p_r_vs_l_indpt   = ranksum(data_R_cv2,data_L_cv2);

%writematrix([data_L_R_cv2;data_R_cv2;data_L_cv2]','\\lexport\iss01.charpier\echanges\aurelie.hanin\all_r_l.csv','Delimiter',';');


%% read precomputed LFP and compute CSD + visual selection of trials
pat_list = [1 3 5];
for ipatient = pat_list%1:size(config,2)
    MuseStruct = readMuseMarkers(config{ipatient},false);
    MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,false); 
    LFP = readLFP(config{ipatient},MuseStruct,false);
    
    load('elec1020_neighb.mat','neighbours');
    clear data* %one data file per patient
    for ipart = 1:size(LFP,2)
        for markername = string(config{ipatient}.LFP.name)
            if ~isfield(LFP{ipart},markername)
                continue
            end
            if isempty(LFP{ipart}.(markername))
                continue
            end
                
                % correct baseline
                cfgtemp                 = [];
                cfgtemp.demean          = 'yes';
                cfgtemp.baselinewindow  = [-2 -1];
                LFP_temp                = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
                               
                 % Compute CSD
                %     neigbours found here :https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
                %     .elec was already in the template folder
                cfgtemp                     = [];
                cfgtemp.method              = 'spline';
                cfgtemp.elec                = ft_read_sens('standard_1020.elc');
                cfgtemp.neighbours          = neighbours;
                CSD_temp                    = ft_scalpcurrentdensity(cfgtemp,LFP_temp);
                
                % select only time of interest for vizualization, remove emg channel
                cfgtemp         = [];
                cfgtemp.latency = [-1 0.5];
                cfgtemp.channel = {'all', '-EMG*'};
                data_sel_LFP    = ft_selectdata(cfgtemp, LFP_temp);
                
                %reject artefacted trials and/or channels
                waitfor(msgbox(sprintf('%s, p%d, %s : reject trials and/or channels',config{ipatient}.prefix(1:end-1),ipart,markername)));
                cfgtemp             = [];
                cfgtemp.method      = 'trial';
                cfgtemp.latency     = [-1 0.5];%only for vizualization (not selection of period)
                cfgtemp.keepchannel = 'no';
                cfgtemp.box         = 'yes';
                cfgtemp.axis        = 'yes';%do not work, so I modified it in the defaults of ft_plot_vector
                data_cleaned        = ft_rejectvisual(cfgtemp,data_sel_LFP);                %apply tyhis selection to CSD
                
                %reject same trials/channels in LFP and CSD whole data
                cfgtemp = [];
                cfgtemp.trials = data_cleaned.cfg.trials;
                cfgtemp.channel = data_cleaned.label;
                data.LFP{ipart}.(markername) = ft_selectdata(cfgtemp,LFP_temp);
                data.CSD{ipart}.(markername) = ft_selectdata(cfgtemp,CSD_temp);
                
                %temp_lfp = ft_timelockanalysis([],data.LFP{ipart}.(markername));
                %cfgtopo.xlim = config{ipatient}.morpho.toi.(markername);
                %ft_topoplotER(cfgtopo,temp_lfp);
                %ft_multiplotER(cfgtopo,data.LFP{ipart}.(markername));
                %ft_topoplotER(cfgtopo,data.CSD{ipart}.(markername));
                %ft_multiplotER(cfgtopo,data.CSD{ipart}.(markername));
        end
    end
    fname = fullfile(config{ipatient}.datasavedir, sprintf('%sSeizures_Visual_Selection.mat',config{ipatient}.prefix));
    fprintf('saving to %s\n', fname);
    save(fname,'data', '-v7.3');
end

%load aligned and selected lfp and csd data
for ipatient = 1:size(config,2)
    fname = fullfile(config{ipatient}.datasavedir, sprintf('%sSeizures_Visual_Selection.mat',config{ipatient}.prefix));
    fprintf('reading %s\n', fname);    
    temp = load(fname);
    data{ipatient} = temp.data;
    clear temp
end

% plot selected trials for each data : LFP/CSD, SlowWave_R/L, topoplot/slowwavemorpho/propagation
%not done yet. 

%% plot average data for each patient (to select toi morpho et topoplot)
%on plot with overdraw trial by trial, and avg of each channel, h = 120.
%one plot with all avg superpos�s

for ipatient = 1:size(config,2)
    for markername = string(config{ipatient}.LFP.name)
        iplot = 0;
        fig = figure('visible','off');
        for idata = ["LFP", "CSD"]
            iplot = iplot+1;
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            data_plot = data{ipatient}.(idata){ipart}.(markername);
            data_plot_avg = ft_timelockanalysis([],data_plot);
            color = [0.6 0.6 0.6];%color of trials when plotted
            if strcmp(idata, 'LFP')
                h = 60;
            elseif strcmp(idata, 'CSD')
                h = 0.01;
            end
            toi_temp = [-2 2];
            
            subplot(1,2,iplot);hold;
            
            %plot each trial
            for itrial = 1:size(data_plot.trialinfo,1)
                
                %plot EEG
                ichan = 0;
                for channame = string(config{ipatient}.LFP.channel)
                    ichan = ichan+1;
                    %select channel
                    cfgtemp = [];
                    cfgtemp.channel = convertStringsToChars(channame);
                    data_1chan = ft_selectdata(cfgtemp,data_plot);
                    if isempty(data_1chan.label)
                        continue
                    end
                    plot(data_1chan.time{itrial},data_1chan.trial{itrial}+(numel(data_plot.label)+1)*h-h*ichan,'Color', color); %first on top
                    %plot average after the last trial
                    if itrial == size(data_plot.trialinfo,1)
                        if ismember(channame, config{ipatient}.(markername).channel)
                            c = 'b';
                        else
                            c = 'k';
                        end
                        chan_idx = strcmp(channame, data_plot_avg.label);
                        plot(data_plot_avg.time, data_plot_avg.avg(chan_idx,:)+(numel(data_plot.label)+1)*h-h*ichan, c, 'LineWidth', 2);
                    end
                    
                end
                
                %plot EMG if any
                cfgtemp = [];
                cfgtemp.channel = {'EMG*'};
                EMG = ft_selectdata(cfgtemp,data_plot);
                if ~isempty(EMG.label)
                    ichan = ichan+1;
                    rescale = 2 * h / (max(EMG.trial{itrial}) - min(EMG.trial{itrial}));
                    plot(EMG.time{itrial}, EMG.trial{itrial}.*rescale + (numel(data_plot.label)+1)*h-h*ichan, 'Color',color); %first on top
                    %plot average after the last trial
                    if itrial == size(data_plot.trialinfo,1)
                        chan_idx = strcmp(EMG.label, data_plot_avg.label);
                        plot(data_plot_avg.time, data_plot_avg.avg(chan_idx,:).*rescale+ (numel(data_plot.label)+1)*h-h*ichan, 'k', 'LineWidth', 2);
                    end
                elseif length(EMG.label)>1
                    error('several EMG channel. It should have only one');
                    
                end
            end
            
            
            %set figure display
            axis tight
            ylim([-h (numel(data_plot.label)+2)*h]);
            xlim(toi_temp);
            xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
            ylabel('Channel name', 'Fontsize',15);
            tick = h;
            yticks(h : tick : numel(data_plot.label)*h);
            chan_name_all = [string(config{ipatient}.LFP.channel), string(EMG.label)];
            set(gca, 'YTickLabel',[flip(chan_name_all')]);
            set(gca, 'FontWeight','bold', 'Fontsize',15);
            set(gca,'TickDir','out');
            
            title(sprintf('%s : %s %s',config{ipatient}.prefix(1:end-1),convertStringsToChars(markername),convertStringsToChars(idata)),'Interpreter','none','Fontsize',18);
            
            %plot vertival line at marker time
            %ax = axis;
            %plot([0 0], [ax(3) ax(4)], '--b');
        end%idata
        %print to file
        imagesavedir = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_overdraw');
        
        if ~(exist (imagesavedir)==7)
            fprintf('Creating dir %s\n',imagesavedir);
            mkdir(imagesavedir);
        end
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        set(fig,'Renderer','Painters');
        print(fig, '-dpdf', fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_overdraw_allseizures.pdf',config{ipatient}.prefix,markername)),'-r600');
        print(fig, '-dpng', fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_overdraw_allseizures.png',config{ipatient}.prefix,markername)),'-r600');
        close all
        
    end%markername
end

%plot average of each channel

for ipatient = 1:size(config,2)
    for markername = string(config{ipatient}.LFP.name)
        if ~isfield(data{ipatient}.(idata){ipart},markername)
            continue
        end
        if isempty(data{ipatient}.(idata){ipart}.(markername))
            continue
        end
        iplot = 0;
        fig = figure('visible','off');
        for idata = ["LFP", "CSD"]
            iplot = iplot+1;
            
            data_plot = data{ipatient}.(idata){ipart}.(markername);
            data_plot_avg = ft_timelockanalysis([],data_plot);
            toi = [-2 2];
            
            subplot(1,2,iplot);hold;
            
            %select only channels of the current amrker
            cfgtemp         = [];
            cfgtemp.channel = config{ipatient}.(markername).channel;
            data_avg_temp   = ft_selectdata(cfgtemp,data_plot_avg);
            
            clear leg
            for ichan = 1:size(data_avg_temp.label,1)
                if ismember(data_avg_temp.label{ichan}, config{ipatient}.morpho.channel.(idata).(markername))
                    leg{ichan} = plot(data_avg_temp.time, data_avg_temp.avg(ichan,:));
                else
                    leg{ichan} = plot(data_avg_temp.time, data_avg_temp.avg(ichan,:),'k');
                end
            end
            
            %add patch of selected period
            ax  = axis;
            toi_temp = config{ipatient}.morpho.toi.(markername);
            x   = [toi_temp(1) toi_temp(2) toi_temp(2) toi_temp(1)];
            y   = [ax(3) ax(3) ax(4) ax(4)];
            patch('XData',x,'YData',y,'facecolor',[0 0 0],'edgecolor','none','facealpha',0.1);
            
            %set figure display
            legend([leg{:}],data_avg_temp.label,'Location','northeastoutside');
            axis tight
            xlim(toi);
            xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
            %ylabel('�V', 'Fontsize',15);
            set(gca, 'FontWeight','bold', 'Fontsize',15);
            set(gca,'TickDir','out');
            
            title(sprintf('%s : %s %s',config{ipatient}.prefix(1:end-1),convertStringsToChars(markername),convertStringsToChars(idata)),'Interpreter','none','Fontsize',18);
        end%idata
        
        imagesavedir = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_avg');
        
        if ~(exist (imagesavedir)==7)
            fprintf('Creating dir %s\n',imagesavedir);
            mkdir(imagesavedir);
        end
        
        fprintf('Print image to %s\n', fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_avg',config{ipatient}.prefix,markername)));
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        set(fig,'Renderer','Painters');
        print(fig, '-dpdf', fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_avg.pdf',config{ipatient}.prefix,markername)),'-r600');
        print(fig, '-dpng', fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_avg.png',config{ipatient}.prefix,markername)),'-r600');
        savefig(fig,fullfile(imagesavedir,sprintf('%s%s_LFP_CSD_avg.fig',config{ipatient}.prefix,markername)));
        close all
    end
end


%% topoplot : one figure for all patients, one figure R and one figure L
%voir si mieux avec toi morpho ou toi topoplot

%topoplot
for markername = string(config{ipatient}.LFP.name)
    for idata = ["LFP", "CSD"]
        fig = figure('visible','off');
        iplot = 0;
        for ipatient = 1:size(config,2)
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            to_plot = ft_timelockanalysis([], data{ipatient}.(idata){ipart}.(markername));
            toi = config{ipatient}.morpho.toi.(markername);%[config{ipatient}.morpho.toi.(markername)(1) - 0.3, config{ipatient}.morpho.toi.(markername)(2) + 0.1];
            
            iplot = iplot+1;
            subplot(4,3,iplot);hold;
            cfgtemp = [];
            cfgtemp.layout        = 'EEG1020';
            cfgtemp.colorbar      = 'yes';
            cfgtemp.zlim          = 'maxabs';
            cfgtemp.xlim          = toi;
            cfgtemp.comment       = 'xlim';
            cfgtemp.fontsize      = 15;
            cfgtemp.renderer      = 'painters';
            cfgtemp.colormap      = flip(parula(5000));%flip so negative peak is in yellow
            cfgtemp.comment = 'no';
            ft_topoplotER(cfgtemp,to_plot);
            %ft_multiplotER(cfgtemp,to_plot);
            
            title(config{ipatient}.prefix(1:end-1), 'Interpreter','none');
            
        end
        %print to file
        fprintf('Print image to %s\n',fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_topoplot_%s_%s',markername,idata)));
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_topoplot_%s_%s.pdf',markername,idata)),'-r600');
        print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_topoplot_%s_%s.png',markername,idata)),'-r600');
        close all
    end
end


%% slow wave morphology
for ipatient = 1:size(config,2)
    for idata = ["LFP", "CSD"]
        for markername = string(config{ipatient}.LFP.name)
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            %pas besoin de s�lectionner channel, d�j� un seul channel dans data
            %cfgmorpho.morpho.channame = config{ipatient}.morpho.channame.(markername);
            
            fig = figure('visible','off');hold;
            %compute data and plot results, for each trial
            for itrial = 1:size(data{ipatient}.(idata){ipart}.(markername).trial,2)
                cfgtemp = [];
                cfgtemp.trials = itrial;
                data_temp = ft_selectdata(cfgtemp, data{ipatient}.(idata){ipart}.(markername));
                
                %plot_morpho([], data_temp);
                config{ipatient}.morpho.toiac = [config{ipatient}.morpho.toi.(markername)(1) - 0.5, config{ipatient}.morpho.toi.(markername)(2)+0.2];
                config{ipatient}.morpho.toibl = [-1.5 -0.5];
                config{ipatient}.morpho.channame = config{ipatient}.align.channel.(markername);
                try
                [hw, ~, ~,amp] = plot_morpho(config{ipatient}, data_temp);
                catch
                    warning('cannot find slowwave in trial %d', itrial);
                    hw = nan;
                    amp = nan;
                end
                
                morpho.(markername).(idata).halfwidth{ipatient}(itrial) = hw;
                morpho.(markername).(idata).amplitude{ipatient}(itrial) = amp;
            end
            %remove text to make the figure readable
            delete(findall(gcf,'type','text'));
            
            %save fig
            if ~isfolder(fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl'))
                fprintf('creating directory %s\n', fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl'));
                mkdir(fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl'))
            end
            fprintf('Print to %s\n', fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl',sprintf('%smorpho_%s_%s',config{ipatient}.prefix,markername,idata)));
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl',sprintf('%smorpho_%s_%s.pdf',config{ipatient}.prefix,markername,idata)),'-r600');
            print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'sw_morpho_shortbl',sprintf('%smorpho_%s_%s.png',config{ipatient}.prefix,markername,idata)),'-r600');
            close all
        end
    end
end
%save computed morpho
fprintf('save morpho values to %s\n',fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'));
save(fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'), 'morpho', '-v7.3');

%plot distrib des amplitudes, et des halfwidth, pour chaque patient
for iparam = ["halfwidth", "amplitude"]
    for idata = ["LFP", "CSD"]
        for markername = string(config{ipatient}.LFP.name)
            fig=figure;hold
            for ipatient = 1:size(config,2)
                param = abs(morpho.(markername).(idata).(iparam){ipatient});
                meanparam.(idata).(markername)(ipatient) = nanmean(param);
                stdparam(ipatient)  = nanstd(param);
                scatter(rand(size(param))*0.2+ipatient-0.1, param, '.', 'MarkerEdgeColor', 'k');
                errorbar(ipatient, meanparam.(idata).(markername)(ipatient), stdparam(ipatient),'--rx');
                %     errorbar(ipatient, nanmean(emg_duration{ipatient}), nanmean(emg_duration{ipatient})/sqrt(size(emg_duration{ipatient},2)),'--rx'); %errbar : sem
            end
            xlim([0 size(config,2)+1]);
            set(gca,'TickDir','out','FontWeight','bold');
            ylabel(iparam);
            
            fprintf('Print to %s\n', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_%s',markername,idata,iparam)));
            fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_%s.pdf',markername,idata,iparam)),'-r600');
            print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_%s.png',markername,idata,iparam)),'-r600');
            close all
        end
    end
end

%plot mean of all patients in the same plot
for iparam = ["halfwidth", "amplitude"]
    for idata = ["LFP", "CSD"]
        fig=figure;hold;
        count_marker = 0;
        for markername = string(config{ipatient}.LFP.name)
            count_marker = count_marker+1;
            scatter(rand(size(meanparam.(idata).(markername)))*0.2+count_marker-0.1, meanparam.(idata).(markername),'sk','MarkerFaceColor','k');
            errorbar(count_marker,nanmean(meanparam.(idata).(markername)),nanstd(meanparam.(idata).(markername)),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
            xlim([0.5 2.5]);
            xticks(1:2);
            xticklabels(strrep(string(config{ipatient}.LFP.name),'_',' '));
            set(gca,'TickDir','out','FontWeight','bold');
            ylabel(sprintf('%s %s',iparam,idata));
            ax = axis; ylim([0 ax(4)]);
        end
        fprintf('Print to %s\n', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_%s_summary',markername,idata,iparam)));
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_summary.pdf',idata,iparam)),'-r600');
        print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_summary.png',idata,iparam)),'-r600');
        close all
    end
end



%% propagation 
%� revoir avec la m�thode de Havenith et al., 2011
%data.(idata){ipart}.(markername)
analysis = "propagation";
ipart = 1;
% datatest = select_trials_visual(config{ipatient},[],ipart,markername, false);

for ipatient = 1:size(config,2)
    for idata = ["LFP", "CSD"]
        for markername = string(config{ipatient}.LFP.name)
            if isfield(data{ipatient}.(idata){ipart},markername)
                
                data_temp = data{ipatient}.(idata){ipart}.(markername);
                                
                %compute average
                data_temp_avg = ft_timelockanalysis([], data_temp);
                %plot(data_temp_avg.time, data_temp_avg.avg)
                %plot(data_temp_avg.time, normalize(data_temp_avg.avg,2,'range'))

                % one plot per markername patient
                cfgtemp                     = config{ipatient};
                cfgtemp.LFP.xcorr.suffix    = sprintf('_%s_%s',markername,idata);
                cfgtemp.LFP.xcorr.xchan     = config{ipatient}.align.channel.(markername);
                cfgtemp.LFP.xcorr.toi       = 'all';%FIXME : voir si ne pas s�lectionner, et choisir une p�riode par d�faut
                cfgtemp.LFP.xcorr.plotdata  = 'yes';
                data_xcorr{ipatient}.(idata).(markername)        = dtx_xcorr_LFP(cfgtemp,data_temp_avg);
                
                % plot avg data 
                fig = figure;
                plot(data_temp_avg.time,data_temp_avg.avg);
                legend(data_temp_avg.label');
                %print to file
                fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'..',sprintf('%s_lfp_for_xcorr_%s_%s.pdf',config{ipatient}.prefix,markername,idata)),'-r600');
                print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'..',sprintf('%s_lfp_for_xcorr_%s_%s.png',config{ipatient}.prefix,markername,idata)),'-r600');
                close all
            end
        end
    end
end

save(fullfile(config{ipatient}.datasavedir,'allpatients_xcorr_lfp.mat'), 'data_xcorr', '-v7.3');
load(fullfile(config{ipatient}.datasavedir,'allpatients_xcorr_lfp.mat'), 'data_xcorr');

%plot distrib des lags pour chaque �lectrode (1 point = 1 patient)
latency_xcorr_all = [];
for idata = ["LFP", "CSD"]
    for markername = string(config{ipatient}.LFP.name)
        
        %gather data from all patients
        for ipatient = 1%:size(config,2)
            if isfield(data{ipatient}.(idata){ipart},markername)
                for ichan =1:size(data_xcorr{ipatient}.(idata).(markername).ylabel,2)
                    channame                            = data_xcorr{ipatient}.(idata).(markername).ylabel{ichan};
                    %initialize structure
                    latency_xcorr_all.(idata)                         = ft_getopt(latency_xcorr_all, convertStringsToChars(idata), []);
                    latency_xcorr_all.(idata).(markername)            = ft_getopt(latency_xcorr_all.(idata),convertStringsToChars(markername), []);
                    latency_xcorr_all.(idata).(markername).(channame) = ft_getopt(latency_xcorr_all.(idata).(markername),channame, []);
                    %get value
                    latency_xcorr_all.(idata).(markername).(channame)(end+1) = data_xcorr{ipatient}.(idata).(markername).lag(ichan);
                end
            end
        end
        
        %one figure per datatype
        fig = figure;hold;
        x_count = 0;
        for channame = string(config{ipatient}.(markername).channel)  
            x_count = x_count +1;
            toplot = latency_xcorr_all.(idata).(markername).(channame);
            scatter(rand(size(toplot))*0.2+x_count-0.1, toplot,'sk','MarkerFaceColor','k');
            errorbar(x_count,nanmean(toplot),nanstd(toplot),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
        end
        
        %figure display
        xlim([0.5 x_count+0.5]);
        xticks(1:x_count);
        set(gca,'TickDir','out','FontWeight','bold','FontSize',15);
        ylabel(sprintf('xcorr peak delay from %s',config{ipatient}.(markername).channel{2}));
        xticklabels(config{ipatient}.(markername).channel);
        ax=axis;
        plot([ax(1) ax(2)],[0 0],'--k');
        
        %print to file
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_xcorr_lfp_%s_%s.pdf',markername,idata)),'-r600');
        print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_xcorr_lfp_%s_%s.png',markername,idata)),'-r600');
        close all
        
    end%markername
end%idata

        
            
        

for ipatient = slurm_task_id
    
%     config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir, 'alignpeak');
    
    % read and align data
    [MuseStruct]                    = readMuseMarkers(config{ipatient}, true);
    [MuseStruct]                    = alignMuseMarkers(config{ipatient},MuseStruct, false);
%     [MuseStruct]                    = alignMuseMarkersXcorr(config{ipatient},MuseStruct, true);
    [LFP]                       = readLFP(config{ipatient}, MuseStruct, false); %dat_LFP{ipart}.(markername)
    
    % flip data
    if ft_getopt(config{ipatient}.LFP, 'flip', false) == true
        %inverse the data according to standard clinical visualisation
        for ipart = 1:size(LFP,2)
            for imarker = 1:size(LFP{ipart},2)
                if ~isempty(LFP{ipart}.(markername))
                    for itrial = 1:size(LFP{ipart}.(markername).trial,2)
                        LFP{ipart}.(markername).trial{itrial} = -LFP{ipart}.(markername).trial{itrial};
                    end
                end
            end
        end
    end
    
    % correct baseline
    for ipart = 1:size(LFP,2)
        for imarker = 1:size(LFP{ipart},2)
            if ~isempty(LFP{ipart}.(markername))
                cfgtemp                 = [];
                cfgtemp.demean          = ft_getopt(config{ipatient}.LFP, 'baseline', 'no');
                cfgtemp.baselinewindow  = ft_getopt(config{ipatient}.LFP, 'baselinewindow', []);
                LFP{ipart}.(markername) = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
            end
        end
    end
    
    % remove trials which intersest BAD markers
    
    [LFP, ~]           = removetrials_MuseMarkers(config{ipatient}, LFP, MuseStruct);
    
    % Compute CSD
    %     neigbours found here :https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
    %     .elec was already in the template folder
    load('elec1020_neighb.mat','neighbours');
    for ipart = 1:size(LFP,2)
        for imarker = 1:size(LFP{ipart},2)
            if ~isempty(LFP{ipart}.(markername))
                cfgtemp                     = [];
                cfgtemp.method              = 'spline';
                cfgtemp.elec                = ft_read_sens('standard_1020.elc');
                cfgtemp.neighbours          = neighbours;
                CSD{ipart}.(markername)     = ft_scalpcurrentdensity(cfgtemp, LFP{ipart}.(markername));
            end
        end
    end
    
   
    
    
    %% Topography
    % Pas de stats pour le moment
    config{ipatient}.topoplot.suffix = 'EEG';
    dtx_plot_SlowWaveTopography(config{ipatient}, LFP);
    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient}, LFP);
    
    config{ipatient}.topoplot.suffix = 'CSD';
    dtx_plot_SlowWaveTopography(config{ipatient}, CSD);
    dtx_plot_SlowWaveTopographyTimecouse(config{ipatient}, CSD);
    
    
    %% hw, amplitude, for both EEG and CSD
    ipart = 1;
    cfgmorpho = config{ipatient};
    
    
    %for those patients : one marker is SW_R, and the other is SW_L
    for imarker = 1:size(LFP{ipart},2)
        
        %select channel of interest
        cfgmorpho.morpho.channame = config{ipatient}.morpho.channame.(markername);
        
        %repeat the same analysis for eeg and csd
        i_analysis=1; %index for suffix
        analysis_type = {'eeg', 'csd'}; 
        for data = [LFP, CSD]
            
            fig = figure;hold;
            %compute data and plot results, for each trial
            for itrial = 1:size(data{ipart}.(markername).trial,2)
                cfgtemp = [];
                cfgtemp.trials = itrial;
                data_temp = ft_selectdata(cfgtemp, data{ipart}.(markername));
                
                [stats.(markername).(analysis_type{i_analysis}).halfwidth(itrial), ~, ~,...
                    stats.(markername).(analysis_type{i_analysis}).amplitude(itrial)] = ...
                    plot_morpho(cfgmorpho, data_temp);
            end
            %remove text to make the figure readable
            delete(findall(gcf,'type','text'));
            
            %save fig
            set(gca,'Fontsize',15);
            if ~(exist(config{ipatient}.imagesavedir)==7)
                mkdir(config{ipatient}.imagesavedir);
                fprintf('Create forlder %s',config{ipatient}.imagesavedir);
            end
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            
            fname = [config{ipatient}.prefix,'p',num2str(ipart),'-Morpho-',config{ipatient}.name.(markername),'_',cfgmorpho.morpho.channame,'_scale',strrep(num2str(config{ipatient}.morpho.toiplot),'  ','_')];
            print(fig, '-dpdf', fullfile(config{ipatient}.imagesavedir,[fname,'_',analysis_type{i_analysis},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(config{ipatient}.imagesavedir,[fname,'_',analysis_type{i_analysis},'.png']),'-r600');
            close all
            
            i_analysis=i_analysis+1; %index for suffix
        end
    end
    
    fname = fullfile(config{ipatient}.datasavedir, [config{ipatient}.prefix, 'stats_EEG.mat']);
    save(fname, stats, '-v7.3');
    return %STOP HERE FOR NOW
    %% count markers : time between seizures, eeg/emg delay, emg duration
    
    %% Propagation : plot canaux positifs, ou tous les canaux, normaliser, et compter d�lai � la main.
    
    %% scatter plot stats pooled entre les patients
    % hw, amplitudes, temps entre 2 crises, d�lai EEG EMG, EMG
    % duration. Faire une fonction pour scatter
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% OLD %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    %% plot morpho SW sans EMG
    %reset save path
    config = config_origin;
    config{ipatient}.imagesavedir = fullfile(config{ipatient}.imagesavedir,'..',sprintf('measurement_%s',config{ipatient}.LFP.name{1}), 'morpho_sw_maxchan');
    % check if images directory exists, if not create
    if ~isfolder(config{ipatient}.imagesavedir)
        ft_notice('creating directory %s', config{ipatient}.imagesavedir);
        mkdir(config{ipatient}.imagesavedir);
    end
    for imarker = 1:length(config{ipatient}.LFP.name)
        if ~isempty(LFP{ipart}.(markername))
            %             toi.toiplot = [-2 2];
            %             toi.toibl = cfg.align.toibaseline.(markername);
            %             toi.toiac = [-1 1];
            %             plot_morpho(config{ipatient},dat_LFP,ipart,imarker,config{ipatient}.LFP.electrodetoplot.(markername), toi,true, false,true);
            cfgtemp                     = [];
            cfgtemp.channame            = config{ipatient}.LFP.electrodetoplot.(markername);
            cfgtemp.plotstd             = 'yes';
            cfgtemp.removeoutliers      = 'no';
            cfgtemp.toiplot             = [-2 2];
            cfgtemp.toibl               = config{ipatient}.align.toibaseline.(markername);
            cfgtemp.toiac               = [-1 1];
            cfgtemp.measurehalfwidth     = 'yes';
            cfgtemp.halfwidthmethod     = 'bl';
            cfgtemp.measurepeaktrough    = 'yes';
            cfgtemp.name                = config{ipatient}.LFP.name.(markername);
            cfgtemp.saveplot            = 'yes';
            cfgtemp.imagesavedir        = config{ipatient}.imagesavedir;
            cfgtemp.prefix              = config{ipatient}.prefix;
            plot_morpho(cfgtemp,LFP{ipart}.(markername));
        end
    end
    
    
      
    %% seizures infos
    load(fullfile(config{1}.datasavedir,'All_patients_Seizure_Infos.mat'),'Seizure_Infos'); %previously computed
    
    %reset save path
    config = config_origin;
    config{1}.imagesavedir = fullfile(config{1}.imagesavedir,'..',sprintf('measurement_%s',config{1}.LFP.name{1}), 'seizures_infos');
    
    % check if images directory exists, if not create
    if ~isfolder(config{1}.imagesavedir)
        ft_notice('creating directory %s', config{1}.imagesavedir);
        mkdir(config{1}.imagesavedir);
    end
    
    dtx_stats_seizures_infos(config{1},Seizure_Infos);
    
end


end





