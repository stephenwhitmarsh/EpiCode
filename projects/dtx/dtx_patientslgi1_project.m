% analysis of LGI1 patients' EEG

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

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

feature('DefaultCharacterSet', 'CP1252'); % To fix bug for weird character problems in reading neurlynx

config = dtx_setparams_patients_lgi1;

setfig = @() set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);

% for ipatient = slurm_task_id%1:size(config,2)
%     MuseStruct = readMuseMarkers(config{ipatient},false);
%     MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,false); 
%     LFP = readLFP(config{ipatient},MuseStruct,true);
% return
% end

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\n************** Read marker data for %s **************\n',config{ipatient}.prefix(1:end-1));
    %read Muse marker and correct eeg file if it is micromed segmented data
    MuseStruct = readMuseMarkers_discontinuousMicromed(config{ipatient},false);
    %concatenate muse marekrs
    MuseStruct_concat = concatenateMuseMarkers(config{ipatient},MuseStruct,false);    %count seizure infos
    
    %analysis only on the last part 
    ipart = size(MuseStruct,2);
    
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_L';
    seizure_timing_L{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat,ipart);
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_R';
    seizure_timing_R{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat,ipart);
        
    %get data for each marker to have more readable script after
    marker.synctime = [];
    marker.clock    = datetime.empty;
    marker.dir      = [];
    slowwave_begin_R(ipatient)      = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_R_begin',        marker);  
    emg_begin_R(ipatient)           = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_R_EMG__START__', marker);  
    emg_end_R(ipatient)             = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_R_EMG__END__',   marker);  
    slowwave_begin_L(ipatient)      = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_L_begin',        marker);  
    emg_begin_L(ipatient)           = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_L_EMG__START__', marker);  
    emg_end_L(ipatient)             = ft_getopt(MuseStruct_concat{ipart}.markers,'SlowWave_L_EMG__END__',   marker); 
    
end

%% plot eeg emg delay, and emg duration, and count emg number
pat_list = 1:size(config,2);
for ipatient = pat_list
        %put right and left seizures together
    slowwave_begin  = [slowwave_begin_R(ipatient).synctime,slowwave_begin_L(ipatient).synctime];
    emg_begin       = [emg_begin_R(ipatient).synctime,emg_begin_L(ipatient).synctime];
    emg_end         = [emg_end_R(ipatient).synctime,emg_end_L(ipatient).synctime];
    
    if isempty(slowwave_begin)
        delays{ipatient} = NaN;
        emg_duration{ipatient} = NaN;
        continue
    end
    
    delays{ipatient}        = emg_begin - slowwave_begin;
    emg_duration{ipatient}  = emg_end - emg_begin;
    eeg_dir{ipatient}           = [slowwave_begin_R(ipatient).dir,slowwave_begin_L(ipatient).dir];
    %     [bins, edges]       = histcounts(delays{ipatient},'BinWidth',0.01);
    %     bins_centers        = (edges(1:end-1)+edges(2:end))/2; %limiteinf + limitesup / 2
    %     bar(bins_centers,bins);
    
    %count emg
    count.n_emg.left(ipatient)  = size(~isnan(emg_begin_L(ipatient).synctime),2);
    count.n_emg.right(ipatient) = size(~isnan(emg_begin_R(ipatient).synctime),2);
    count.n_emg.all(ipatient)   = count.n_emg.left(ipatient) + count.n_emg.right(ipatient);

end
count.n_emg.total = sum(count.n_emg.all);
count.n_emg.med = median(count.n_emg.all);
count.n_emg.med_non_zero = median(count.n_emg.all(count.n_emg.all>0));

%eeg emg delay
fig = figure;hold on 
for ipatient = pat_list
    meandelay(ipatient) = nanmean(delays{ipatient});
    stddelay(ipatient)  = nanstd(delays{ipatient});
    bar(ipatient, meandelay(ipatient),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(delays{ipatient}))*0.2+ipatient-0.1, delays{ipatient}, '.', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, meandelay(ipatient),0, stddelay(ipatient),'k','CapSize',10);
    %     errorbar(ipatient, nanmean(delays{ipatient}), nanmean(delays{ipatient})/sqrt(size(delays{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 pat_list(end)+1]);
ylabel('eeg-emg delay (s)');
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_eeg-emg_delay');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%eeg emg delay summary
fig = figure;hold on
x = (1:size(config,2))./size(config,2);
errorbar(x, meandelay,stddelay,'sk','MarkerFaceColor','k');
errorbar(0.5,nanmean(meandelay),nanstd(meandelay),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
xlim([-1 2]);
xticks([]);
ylabel('mean eeg-emg delay, per patient (s)');
ax = axis; ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_eeg-emg_delay_summary');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%emg duration
fig = figure;hold
for ipatient = pat_list
    meanemgduration(ipatient) = nanmean(emg_duration{ipatient});
    stdemgduration(ipatient)  = nanstd(emg_duration{ipatient});
    bar(ipatient, meanemgduration(ipatient),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(emg_duration{ipatient}))*0.2+ipatient-0.1, emg_duration{ipatient}, '.', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, meanemgduration(ipatient),0, stdemgduration(ipatient),'k','CapSize',10);
    %     errorbar(ipatient, nanmean(emg_duration{ipatient}), nanmean(emg_duration{ipatient})/sqrt(size(emg_duration{ipatient},2)),'--rx'); %errbar : sem
end
xlim([0 pat_list(end)+1]);
ylabel('emg duration (s)');
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_emg_duration');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%emg duration summary
fig = figure;hold
x = (1:size(config,2))./size(config,2);
errorbar(x, meanemgduration,stdemgduration,'sk','MarkerFaceColor','k');
errorbar(0.5,nanmean(meanemgduration),nanstd(meanemgduration),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
xlim([-1 2]);
xticks([]);
ylabel('mean emg duration, per patient (s)');
ax = axis; ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_emg_duration_summary');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

% %check dir data
% figure;
% subplot(2,1,1);hold;
% scatter(eeg_dir{ipatient},emg_duration{ipatient});
% title(sprintf('emg duration patient %d',ipatient));
% subplot(2,1,2);hold;
% scatter(eeg_dir{ipatient},delays{ipatient});
% title(sprintf('eeg-emg delays patient %d',ipatient));
%svg stats en output

%% plot seizure frequency and regularity, and count seizure number
%recover real time in segmented data OK
%correct concatenateMuseMarkers : clocktime non changï¿½, synctime en
%fonction du nb de samples. A ne pas utiliser si les fichiers ne sont pas
%continus.

%count seizures
for ipatient = 1:size(config,2)
    count.n_seizures.left(ipatient)  = size(~isnan(seizure_timing_L{ipatient}.time_start.synctime),2);
    count.n_seizures.right(ipatient) = size(~isnan(seizure_timing_R{ipatient}.time_start.synctime),2);
    count.n_seizures.all(ipatient)   = count.n_seizures.left(ipatient) + count.n_seizures.right(ipatient);
end
count.n_seizures.total = sum(count.n_seizures.all);
count.n_seizures.med = median(count.n_seizures.all);

%time between 2 seizures
fig=figure;hold
for ipatient = 1:size(config,2)
    %sort data and remove times between 2 dirs
    [data, indexes] = sort([seizure_timing_R{ipatient}.time_start.clock, seizure_timing_L{ipatient}.time_start.clock]);
    dirtemp = [seizure_timing_R{ipatient}.time_start.dir, seizure_timing_L{ipatient}.time_start.dir];
    dir_list = dirtemp(indexes);
    data = minutes(diff(data));
    data_L_R{ipatient} = data(~logical(diff(dir_list)));
    %plot
    bar(ipatient, nanmean(data_L_R{ipatient}),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(data_L_R{ipatient}))*0.2+ipatient-0.1, data_L_R{ipatient}, '.', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, nanmean(data_L_R{ipatient}), 0,nanstd(data_L_R{ipatient}),'k','CapSize',10);
end
xlim([0 size(config,2)+1]);
ylabel('time between 2 slow waves (right and left, minutes)');
ax = axis;
ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_R_L');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%time between 2 seizures right
fig=figure;hold
for ipatient = 1:size(config,2)
    clear data
    %sort data and remove times between 2 dirs
    data_R{ipatient} = seizure_timing_R{ipatient}.time_start.clock;
    dir_list = seizure_timing_R{ipatient}.time_start.dir;
    data_R{ipatient} = minutes(diff(data_R{ipatient}));
    data_R{ipatient} = data_R{ipatient}(~logical(diff(dir_list)));
    %plot
    bar(ipatient, nanmean(data_R{ipatient}),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(data_R{ipatient}))*0.2+ipatient-0.1, data_R{ipatient}, '.', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, nanmean(data_R{ipatient}), 0,nanstd(data_R{ipatient}),'k','CapSize',10);
end
% for ipatient = 1:size(config,2)
%     data = minutes([seizure_timing_R{ipatient}.time_start.synctime]);
%     scatter(rand(size(data))*0.2+ipatient-0.1, data, '.', 'MarkerEdgeColor', 'k');
%     errorbar(ipatient, nanmean(data), nanstd(data),'--rx');
% end
xlim([0 size(config,2)+1]);
ylabel('time between 2 slow waves (right only, minutes)');
ax = axis;
ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_R');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%time between 2 seizures left
fig=figure;hold
for ipatient = 1:size(config,2)
    clear data
    %sort data and remove times between 2 dirs
    data_L{ipatient} = seizure_timing_L{ipatient}.time_start.clock;
    dir_list = seizure_timing_L{ipatient}.time_start.dir;
    data_L{ipatient} = minutes(diff(data_L{ipatient}));
    data_L{ipatient} = data_L{ipatient}(~logical(diff(dir_list)));
    %plot
    bar(ipatient, nanmean(data_L{ipatient}),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(data_L{ipatient}))*0.2+ipatient-0.1, data_L{ipatient}, '.', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, nanmean(data_L{ipatient}), 0,nanstd(data_L{ipatient}),'k','CapSize',10);
end
xlim([0 size(config,2)+1]);
ylabel('time between 2 slow waves (left only, minutes)');
ax = axis;
ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_L');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%compare right and left
clear meandata stddata
meandata{1} = cellfun(@(v) nanmean((v)),data_R);
meandata{2} = cellfun(@(v) nanmean((v)),data_L);
fig=figure;hold on
for ipatient = 1:size(config,2)
    %errorbar([1 2], [meandata{1}(ipatient), meandata{2}(ipatient)],[stddata{1}(ipatient), stddata{2}(ipatient)],'-sk','MarkerFaceColor','k');
    plot([1 2], [meandata{2}(ipatient), meandata{1}(ipatient)],'-ok','MarkerEdgeColor','k','MarkerSize',10);
end
% errorbar([1 2],nanmean([meandata{1};meandata{2}],2)',nanstd([meandata{1};meandata{2}],0,2)','--rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Left', 'Right'});
ylabel('time between 2 slow waves, avg per patient');
ax = axis; ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_comparison_R_L');
dtx_savefigure(fig,fname,'pdf','png','close');


%time between 2 seizures summary
clear meandata stddata
meandata{1} = cellfun(@(v) nanmean((v)),data_R);
meandata{2} = cellfun(@(v) nanmean((v)),data_L);
meandata{3} = cellfun(@nanmean,data_L_R);
stddata{1} = cellfun(@(v) nanstd((v)),data_R);
stddata{2} = cellfun(@(v) nanstd((v)),data_L);
stddata{3} = cellfun(@nanstd,data_L_R);
fig=figure;hold
count_plot = 0;
for idata = 1:3
    count_plot = count_plot+2;
    x = (1:size(config,2))./size(config,2) + count_plot;
    errorbar(x, meandata{idata},stddata{idata},'sk','MarkerFaceColor','k');
    errorbar(count_plot+0.5,nanmean(meandata{idata}),nanstd(meandata{idata}),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
end
xlim([1 8]);
xticks(2.5:2:6.5);
xticklabels({'R', 'L','R + L'});
ylabel('time between 2 slow waves, mean per patient');
ax = axis; ylim([0 ax(4)]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_summary');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%cv2 of right and left, only right, only left
for ipatient = 1:size(config,2)
    data_L_R_cv2.mean(ipatient) = nanmean(cv2(data_L_R{ipatient}));
    data_L_R_cv2.std(ipatient) = nanstd(cv2(data_L_R{ipatient}));
    data_L_cv2.mean(ipatient) = nanmean(cv2(data_L{ipatient}));
    data_L_cv2.std(ipatient) = nanstd(cv2(data_L{ipatient}));
    data_R_cv2.mean(ipatient) = nanmean(cv2(data_R{ipatient}));
    data_R_cv2.std(ipatient) = nanstd(cv2(data_R{ipatient}));
end
fig=figure;hold
scatter_position = 0;
for idata = [data_R_cv2,data_L_cv2]% [data_L_R_cv2,data_R_cv2,data_L_cv2]
    scatter_position = scatter_position +1;
    scatter(rand(size(idata.mean))*0.2+scatter_position-0.1, idata.mean, 's', 'filled','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    errorbar(scatter_position, nanmean(idata.mean), nanstd(idata.mean),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
end
% plot([1 2 3], [nanmean(data_L_R_cv2.mean),nanmean(data_R_cv2.mean),nanmean(data_L_cv2.mean)], '--r', 'LineWidth',2);
% xlim([0.5 3.5]);
% xticks(1:3);
% xticklabels({'L + R', 'R only', 'L only'});
xlim([0.5 2.5]);
xticks(1:2);
xticklabels({'R', 'L'});
ylabel('Patients'' seizures cv2');
ylim([0 2]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_cv2_slowwaves_r_l');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%line plot cv2
fig=figure;hold
for ipatient = 1:size(config,2)%{data_L_R_cv2,data_L_cv2,data_R_cv2}
    plot([1 2 3], [data_L_R_cv2.mean(ipatient),data_R_cv2.mean(ipatient),data_L_cv2.mean(ipatient)],'-sk','MarkerFaceColor','k');
end
errorbar([1 2 3], nanmean([data_L_R_cv2.mean;data_R_cv2.mean;data_L_cv2.mean],2)',nanstd([data_L_R_cv2.mean;data_R_cv2.mean;data_L_cv2.mean],0,2)','--rx','LineWidth',2);
xlim([0.5 3.5]);
ylabel('Patients'' seizures cv2');
xticks(1:3);
xticklabels({'L + R', 'R only', 'L only'});
ylim([0 2]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_cv2_slowwaves_lineplot');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%summary cv2 right and left
clear meandata stddata
fig=figure;hold
count_plot = 0;
for idata = [data_R_cv2,data_L_cv2]
    count_plot = count_plot+2;
    x = (1:size(config,2))./size(config,2) + count_plot;
    errorbar(x, idata.mean,idata.std,'sk','MarkerFaceColor','k');
    errorbar(count_plot+0.5,nanmean(idata.mean),nanstd(idata.mean),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
end
xlim([1 6]);
xticks([2.5 4.5]);
xticklabels({'R', 'L'});
ylabel('cv2 of interseizure interval, mean per patient');
ylim([0 2]);
setfig();
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_CV2_summary');
dtx_savefigure(fig,fname,'pdf','png','fig','close');



%signrank : paired test. ranksum(x,y) : x et y indï¿½pendants
p_all_vs_l = signrank(data_L_R_cv2,data_L_cv2);
p_all_vs_r = signrank(data_L_R_cv2,data_R_cv2);
p_r_vs_l   = signrank(data_R_cv2,data_L_cv2);
p_all_vs_l_indpt = ranksum(data_L_R_cv2,data_L_cv2);
p_all_vs_r_indpt = ranksum(data_L_R_cv2,data_R_cv2);
p_r_vs_l_indpt   = ranksum(data_R_cv2,data_L_cv2);

% %writematrix([data_L_R_cv2;data_R_cv2;data_L_cv2]','\\lexport\iss01.charpier\echanges\aurelie.hanin\all_r_l.csv','Delimiter',';');

%% plot slow waves over time
fig = figure;
for ipatient = 1:size(config,2)
    subplot(size(config,2),1,ipatient); hold on;
  
    if seizure_timing_R{ipatient}.nrseizures >0
        length_previous_dirs(1) = 0;
        for idir = 1:size(MuseStruct{ipatient}{1},2)
            start_dir = MuseStruct{ipatient}{1}{idir}.starttime;
            dir_length = minutes(MuseStruct{ipatient}{1}{idir}.endtime - MuseStruct{ipatient}{1}{idir}.starttime);
            if ismember(idir,unique(seizure_timing_R{ipatient}.time_start.dir))
                for i_seizure= find(seizure_timing_R{ipatient}.time_start.dir == idir)
                    marker_time = seizure_timing_R{ipatient}.time_start.clock(i_seizure);
                    t = minutes(marker_time-start_dir) + length_previous_dirs(idir);
                    plot([t,t],[-1,1],'-b');
                end
            end
            length_previous_dirs(idir+1) = length_previous_dirs(idir) + dir_length;
        end
    end
    
    %left slow waves
    if seizure_timing_L{ipatient}.nrseizures >0
        length_previous_dirs(1) = 0;
        for idir = 1:size(MuseStruct{ipatient}{1},2)
            start_dir = MuseStruct{ipatient}{1}{idir}.starttime;
            dir_length = minutes(MuseStruct{ipatient}{1}{idir}.endtime - MuseStruct{ipatient}{1}{idir}.starttime);
            if ismember(idir,unique(seizure_timing_L{ipatient}.time_start.dir))
                for i_seizure= find(seizure_timing_L{ipatient}.time_start.dir == idir)
                    marker_time = seizure_timing_L{ipatient}.time_start.clock(i_seizure);
                    t = minutes(marker_time-start_dir) + length_previous_dirs(idir);
                    plot([t,t],[-1,1],'-r');
                end
            end
            length_previous_dirs(idir+1) = length_previous_dirs(idir) + dir_length;
        end
    end
    
    %plot limit betweeen 2 files
    for idir = 1:length(MuseStruct{ipatient}{ipart})
        if idir >1
            t = length_previous_dirs(idir);
            plot_dirlim = plot([t,t],[-2,2],'Color',[0.5 0.5 0.5],'Linewidth', 2, 'LineStyle', '-');
        end
    end
    
    ax = axis;
    xlim([0 ax(2)]);
    yticks([]);
    ylim([-2 2]);    
    set(gca,'TickDir','out','FontWeight','bold')
    if ipatient == size(config,2)
        xlabel('Time (minute)');
    end
end

fname = fullfile(config{ipatient}.imagesavedir,'..','allpatients_slowwaves_over_time');
dtx_savefigure(fig,fname,'pdf','png','close');

%plot only patient 2 and 19
pat_list = [2 4];
dir_list = [2 1];
for i = 1:2
    ipatient = pat_list(i);
    idir = dir_list(i);
    fig = figure;
    subplot(5,1,1); hold on;
    
    for i_seizure= find(seizure_timing_R{ipatient}.time_start.dir == idir)
        t = seizure_timing_R{ipatient}.time_start.clock(i_seizure);
        leg_R = plot([t,t],[-1,1],'-b','LineWidth',2);
    end
    
    for i_seizure= find(seizure_timing_L{ipatient}.time_start.dir == idir)
        t = seizure_timing_L{ipatient}.time_start.clock(i_seizure);
        leg_L = plot([t,t],[-1,1],'-r','LineWidth',2);
    end
    
    % xlim([0 ax(2)]);
    yticks([]);
    ylim([-2 2]);
    setfig();
    xlabel('Time');
    
    legend([leg_R, leg_L], 'SlowWave R', 'SlowWave L','location', 'eastoutside');
    
    fname = fullfile(config{ipatient}.imagesavedir,sprintf('%sslowwaves_over_time',config{ipatient}.prefix));
    dtx_savefigure(fig,fname,'pdf','png','close');
end

% figure;hold;
% for ipatient = 1:size(config,2)
%     if ~isempty(seizure_timing_L{ipatient}.time_start.synctime) && ~isempty(seizure_timing_R{ipatient}.time_start.synctime)
%         [X2, lag] = xcorr(seizure_timing_L{ipatient}.time_start.synctime,seizure_timing_R{ipatient}.time_start.synctime);
%         plot(lag,X2);
%     end
% end

%% #DATA# read and align LFP, compute CSD, visual rejection of artefacts (load precomputed data)
if false
    for ipatient = 1:size(config,2)
        MuseStruct = readMuseMarkers(config{ipatient},false);
        MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,true);
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
else
    %load aligned and selected lfp and csd data
    clear data
    for ipatient = 1:size(config,2)
        fname = fullfile(config{ipatient}.datasavedir, sprintf('%sSeizures_Visual_Selection.mat',config{ipatient}.prefix));
        fprintf('reading %s\n', fname);
        temp = load(fname);
        data{ipatient} = temp.data;
        clear temp
    end
end

%% plot average data for each patient (to select toi morpho and topoplot)
%on plot with overdraw trial by trial, and avg of each channel, h = 120.
%one plot with all avg superposï¿½s

for ipatient = 13%:size(config,2)
    for markername = string(config{ipatient}.LFP.name(1:2))
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
        fname = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_overdraw',sprintf('%s%s_LFP_CSD_overdraw_allseizures',config{ipatient}.prefix,markername));
        dtx_savefigure(fig,fname,'pdf','png','close');
        
    end%markername
end

%plot average of each channel

for ipatient = 13%1:size(config,2)
    for markername = string(config{ipatient}.LFP.name(1:2))
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
            %ylabel('ï¿½V', 'Fontsize',15);
            set(gca, 'FontWeight','bold', 'Fontsize',15);
            set(gca,'TickDir','out');
            
            title(sprintf('%s : %s %s',config{ipatient}.prefix(1:end-1),convertStringsToChars(markername),convertStringsToChars(idata)),'Interpreter','none','Fontsize',18);
        end%idata
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir, '..', 'plot_each_seizure_avg',sprintf('%s%s_LFP_CSD_avg',config{ipatient}.prefix,markername));
        dtx_savefigure(fig,fname,'pdf','png','fig','close');
        
    end
end

%% plot topo : one figure for all patients, one figure R and one figure L
%ok pour les 13 patients

%topoplot
for markername = string(config{ipatient}.LFP.name(1:2))
    for idata = ["LFP", "CSD"]
        fig = figure('visible','off');
        iplot = 0;
        for ipatient = 1:size(config,2)
            fprintf('topoplot %s : %s, patient %d\n', idata, markername, ipatient);
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            to_plot = ft_timelockanalysis([], data{ipatient}.(idata){ipart}.(markername));
            toi = config{ipatient}.morpho.toi.(markername);%[config{ipatient}.morpho.toi.(markername)(1) - 0.3, config{ipatient}.morpho.toi.(markername)(2) + 0.1];
            
            iplot = iplot+1;
            subplot(4,3,iplot);hold on;
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
            
            title(sprintf('pat %d (n=%d)',ipatient, size(data{ipatient}.(idata){ipart}.(markername).trial,2)), 'Interpreter','none');
            
        end
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir, '..',sprintf('allpatients_topoplot_%s_%s',markername,idata));
        dtx_savefigure(fig,fname,'pdf','png','close');
        
    end
end

%% #DATA# compute slow wave morphology (load precomputed data)
%ok pour les 13 patients
if false
    for ipatient = 1:size(config,2)
        for idata = ["LFP", "CSD"]
            for markername = string(config{ipatient}.LFP.name)
                if ~isfield(data{ipatient}.(idata){ipart},markername)
                    continue
                end
                if isempty(data{ipatient}.(idata){ipart}.(markername))
                    continue
                end
                
                %pas besoin de sï¿½lectionner channel, dï¿½jï¿½ un seul channel dans data
                %cfgmorpho.morpho.channame = config{ipatient}.morpho.channame.(markername);
                
                fig = figure('visible','off');hold;
                %compute data and plot results, for each trial
                for itrial = 1:size(data{ipatient}.(idata){ipart}.(markername).trial,2)
                    cfgtemp = [];
                    cfgtemp.trials = itrial;
                    data_temp = ft_selectdata(cfgtemp, data{ipatient}.(idata){ipart}.(markername));
                    
                    %plot_morpho([], data_temp);
                    config{ipatient}.morpho.toiac = [config{ipatient}.morpho.toi.(markername)(1) - 0.5, config{ipatient}.morpho.toi.(markername)(2)+0.2];
                    %config{ipatient}.morpho.toibl = [-1.5 -0.5];
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
                
                %print to file
                fname = fullfile(config{ipatient}.imagesavedir,'sw_morpho',sprintf('%smorpho_%s_%s',config{ipatient}.prefix,markername,idata));
                dtx_savefigure(fig,fname,'pdf','png','fig','close');
                
            end
        end
    end
    %save computed morpho
    fprintf('save morpho values to %s\n',fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'));
    save(fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'), 'morpho', '-v7.3');
else
    %plot distrib des amplitudes, et des halfwidth, pour chaque patient
    load(fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'), 'morpho');
end

%% plot slow wave morphology
for iparam = ["halfwidth", "amplitude"]
    for idata = ["LFP", "CSD"]
        count_marker = 0;
        fig=figure;hold on
        for markername = ["SlowWave_L", "SlowWave_R"]
            if strcmp(markername,"SlowWave_R")
            	count_marker = count_marker + size(config,2) + 5;
            end
            for ipatient = 1:size(config,2)
                if size(morpho.(markername).(idata).(iparam),2) < ipatient
                    meanparam.(idata).(markername).(iparam)(ipatient) = nan;
                    stdparam.(idata).(markername).(iparam)(ipatient) = nan;
                    continue
                end
                if isempty(morpho.(markername).(idata).(iparam){ipatient})
                    meanparam.(idata).(markername).(iparam)(ipatient) = nan;
                    stdparam.(idata).(markername).(iparam)(ipatient) = nan;
                    continue
                end
                
                param = abs(morpho.(markername).(idata).(iparam){ipatient});
                meanparam.(idata).(markername).(iparam)(ipatient) = nanmean(param);
                stdparam.(idata).(markername).(iparam)(ipatient)  = nanstd(param);
                
                
                bar(ipatient+ count_marker, meanparam.(idata).(markername).(iparam)(ipatient),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
                scatter(rand(size(param))*0.2+ipatient-0.1 + count_marker, param, '.', 'MarkerEdgeColor', 'k');
                errorbar(ipatient+ count_marker, meanparam.(idata).(markername).(iparam)(ipatient),0, stdparam.(idata).(markername).(iparam)(ipatient),'k');
                
            end
            
        end
        
        ax = axis;
        ylim([0 ax(4)]);
        xlim([-4, size(config,2)+count_marker+5]);
        xticks([1:size(config,2), count_marker+1:count_marker+size(config,2)]);
        xticklabels([1:size(config,2), 1:size(config,2)]);
        setfig();
        ylabel(iparam);
        set(gca, 'Fontsize', 11);
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s',idata,iparam));
        dtx_savefigure(fig,fname,'pdf','png','close');
        
    end
end

%compare right and left
for iparam = ["halfwidth", "amplitude"]
    for idata = ["LFP", "CSD"]
        fig=figure;hold on
        for ipatient = 1:size(config,2)
            plot([1 2], [meanparam.(idata).SlowWave_L.(iparam)(ipatient), meanparam.(idata).SlowWave_R.(iparam)(ipatient)],'-ok','MarkerEdgeColor','k','MarkerSize',10);
        end
        % errorbar([1 2],nanmean([meandata{1};meandata{2}],2)',nanstd([meandata{1};meandata{2}],0,2)','--rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels({'Left', 'Right'});
        ylabel(sprintf('%s : avg of each patient', iparam));
        ax = axis; ylim([0 ax(4)]);
        setfig();
        fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_compare_R_L_%s_%s',idata,iparam));
        dtx_savefigure(fig,fname,'pdf','png','close');
    end
end

%plot mean of all patients in the same plot
for iparam = ["halfwidth", "amplitude"]
    for idata = ["LFP", "CSD"]
        fig=figure;hold;
        count_marker = 0;
        for markername = string(config{ipatient}.LFP.name(1:2))
            count_marker = count_marker+1;
            %scatter(rand(size(meanparam.(idata).(markername).(iparam)))*0.2+count_marker-0.1, meanparam.(idata).(markername).(iparam),'sk','MarkerFaceColor','k');
            nb_points = size(meanparam.(idata).(markername).(iparam),2);
            x = (1:nb_points)./ nb_points .* 0.6 + 0.7 + (count_marker-1);
            errorbar(x, meanparam.(idata).(markername).(iparam),stdparam.(idata).(markername).(iparam),'sk','MarkerFaceColor','k');
            errorbar(count_marker,nanmean(meanparam.(idata).(markername).(iparam)),nanstd(meanparam.(idata).(markername).(iparam)),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
            xlim([0.5 2.5]);
            xticks(1:2);
            xticklabels(strrep(string(config{ipatient}.LFP.name),'_',' '));
            set(gca,'TickDir','out','FontWeight','bold');
            ylabel(sprintf('%s %s',iparam,idata));
            ax = axis; ylim([0 ax(4)]);
        end
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_%s_%s_summary',idata,iparam));
        dtx_savefigure(fig,fname,'pdf','png','close');
            
    end
end

%% #DATA# read and align EMG data (load precomputed data)
markers_emg_begin = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align"];

if false
    %get EMG data
    clear data_EMG
    for ipatient = 1:size(config,2)
        MuseStruct = readMuseMarkers(config{ipatient},true);
        MuseStruct = alignMuseMarkers(config{ipatient},MuseStruct,true);
        LFP = readLFP(config{ipatient},MuseStruct,true);
        %remove EMG with BAD markers
        LFP = removetrials_MuseMarkers([],LFP,MuseStruct,true);
        for markername = markers_emg_begin
            if isempty(LFP{1}.(markername))
                continue
            end
            cfgtemp                         = [];
            cfgtemp.channel                 = config{ipatient}.EMG.(markername);
            data_EMG{ipatient}.(markername) = ft_selectdata(cfgtemp, LFP{1}.(markername));
        end
    end
    save(fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'), 'data_EMG', '-v7.3');
else
    fprintf('Reading %s\n',fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'));
    load(fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'), 'data_EMG');
end

%% EMG : compute envelopes
for markername = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align"]
    
    for ipatient = 1:size(config,2)
        if ipatient > size(data_EMG,2)
                continue
            end
        if ~isfield(data_EMG{ipatient}, markername)
            continue
        end
        if isempty(data_EMG{ipatient}.(markername))
            continue
        end
        if isempty(data_EMG{ipatient}.(markername).label)
            continue
        end
        
        cfgtemp = [];
        cfgtemp.channel = config{ipatient}.EMG.(markername);
        EMG = ft_selectdata(cfgtemp,data_EMG{ipatient}.(markername));
        
        if isempty(EMG.label)
            continue
        end
        
        t = EMG.time{1};
        
        %compute all envelopes, create a fieldtrip structure to average
        for itrial = 1 : size(EMG.trial,2)
            rect_emg = abs(EMG.trial{itrial}(1,:));
            [env{ipatient}.trial{itrial}, ~] = envelope(rect_emg,config{ipatient}.EMG.envparam,config{ipatient}.EMG.envmethod);
            env{ipatient}.time{itrial} = t;
        end
        env{ipatient}.label = {'dummy'};
        
        env_avg{ipatient} = ft_timelockanalysis([], env{ipatient});

        %value to do normalize later
        bl_avg(ipatient) = nanmean(env_avg{ipatient}.avg(t>-2&t<-0.5));
        norm_factor(ipatient) = max(env_avg{ipatient}.avg - bl_avg(ipatient));

    end
end


%% plot average EMG morpho
markers_emg_begin = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align"];

%assume all trials from all patients have the same length
t = data_EMG{1}.SlowWave_R_EMG_align.time{1};

for plot_std = [true false]
    
    count_plot = 0;
    fig = figure;hold;
    for markername = markers_emg_begin
        %create nan array so if patients are missing they are put as nan
        %env_avg.(markername).avg = nan(size(config,2), size(t,2));
        
        for ipatient = 1:size(config,2)
            
            if ipatient > size(data_EMG,2)
                continue
            end
            if ~isfield(data_EMG{ipatient}, markername)
                continue
            end
            if isempty(data_EMG{ipatient}.(markername).label)
                continue
            end
            if size(data_EMG{ipatient}.(markername).label,1) > 1
                error('supposed to have only 1 EMG channel');
            end
            t = data_EMG{ipatient}.(markername).time{1};
            
            %plot trial by trial : rect and envelope
            %         for itrial = 1 : size(data_EMG{ipatient}.(markername).trial,2)
            %             t = data_EMG{ipatient}.(markername).time{itrial};
            %             rect_emg = abs(data_EMG{ipatient}.(markername).trial{itrial}(1,:));
            %             plot(t,rect_emg,'color',[0.6 0.6 0.6]); %first on top
            %         end
            
            %compute all envelopes
%             env = [];
%             for itrial = 1 : size(data_EMG{ipatient}.(markername).trial,2)
%                 rect_emg = abs(data_EMG{ipatient}.(markername).trial{itrial}(1,:));
%                 [env{itrial}, ~] = envelope(rect_emg,config{ipatient}.EMG.envparam,config{ipatient}.EMG.envmethod);
%                 %plot(t,env{itrial},'color','k');
%             end
            
%             %average envelopes
%             env_avg.(markername).time = t; %suppose that all trials have the same length
%             for ioffset = 1:length(env{1}) %all trials must have the same size
%                 for itrial = 1:length(env)
%                     env_by_offset(itrial) = env{itrial}(ioffset);
%                 end
%                 env_avg.(markername).avg(ipatient,ioffset) = mean(env_by_offset);
%             end
%             env_avg.(markername).avg(ipatient,:) = normalize(env_avg.(markername).avg(ipatient,:), 'range');
            
            %normalize
            env_avg_norm = (env_avg{ipatient}.avg - bl_avg(ipatient)) ./ norm_factor(ipatient);
            env_std_norm = sqrt(env_avg{ipatient}.var) ./ norm_factor(ipatient);
            
            count_plot = count_plot+1;
            
            if plot_std
                h=2;
            else
                h=1;
            end
            
            %plot avg
            p = plot(env_avg{ipatient}.time,env_avg_norm+count_plot*h, 'k');
            p.ZData = ones(size(p.YData));
            
            if plot_std
                %plot std
                std_data = env_std_norm;
                x = env_avg{ipatient}.time;
                y = [env_avg_norm+count_plot*h - std_data; std_data; std_data]';
                filled_SD = area(x,y);
                filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
                filled_SD(1).ShowBaseLine = 'off';
            end
            
            
            %plot(env_avg.(markername).time,env_avg.(markername).avg(ipatient,:)+count_plot*h, 'k');
            y_labels{count_plot} = sprintf('pat %d (n=%d)',ipatient, size(data_EMG{ipatient}.(markername).trial,2));
            patient{count_plot} = config{ipatient}.prefix(1:end-1);
        end
    end
    
    %set figure display
    ax = axis;
    xlim([-2 2])
    ylim([h-h/4 ax(4)]);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    set(gca,'TickDir','out');
    yticks(h:h:count_plot*h);
    yticklabels(y_labels);
    
    if plot_std
        title('EMG morphology : average+/-std of envelopes of each trial','Fontsize',18);
        fname = fullfile(config{ipatient}.imagesavedir,'..', 'allpatients_emg_morphology_withstd2');
    else
        title('EMG morphology : average of envelopes of each trial','Fontsize',18);
        fname = fullfile(config{ipatient}.imagesavedir,'..', 'allpatients_emg_morphology');
    end
        
    dtx_savefigure(fig,fname, 'png','pdf','close');
end

%% plot EMG method 
markers_emg_begin = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align"];
count_plot = 0;
for markername = markers_emg_begin
    for ipatient = 1:size(config,2)
        if ipatient > size(data_EMG,2)
            continue
        end
        if ~isfield(data_EMG{ipatient}, markername)
            continue
        end
        if isempty(data_EMG{ipatient}.(markername).label)
            continue
        end
        if size(data_EMG{ipatient}.(markername).label,1) > 1
            error('supposed to have only 1 EMG channel');
        end
        
        fig = figure;hold;
        t = data_EMG{1}.(markername).time{1};
        
        %compute all envelopes
        env = [];
        for itrial = 1 : size(data_EMG{ipatient}.(markername).trial,2)
            rect_emg = abs(data_EMG{ipatient}.(markername).trial{itrial}(1,:));
            [env{itrial}, ~] = envelope(rect_emg,config{ipatient}.EMG.envparam,config{ipatient}.EMG.envmethod);
            %plot(t,env{itrial},'color','k');
        end
        
        nb_trials = size(data_EMG{ipatient}.(markername).trial,2);
        
        %h automatic setting :
        for itrial = 1 : nb_trials
            h_temp_max = max(data_EMG{ipatient}.(markername).trial{itrial}(1,:));
            h_temp_min = min(data_EMG{ipatient}.(markername).trial{itrial}(1,:));
            h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
        end
        h = mean(h_temp_amplitude);
        
        % plot raw emg
        subplot(1,3,1); hold on
        for itrial = 1 : nb_trials
            plot(t,data_EMG{ipatient}.(markername).trial{itrial}(1,:)+(nb_trials+1)*h- itrial*h,'k'); %first on top
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        %xlabel(sprintf('Time from \n%s (s)', config{ipatient}.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('Number of seizures');
        title(sprintf('Timecourse of \n%s (%d trials)', data_EMG{ipatient}.(markername).label{1}, size(data_EMG{ipatient}.(markername).trial,2)),'Interpreter','none','Fontsize',15);
        setfig();
        tick = h;
        yticks(tick : tick*10 : nb_trials*h);
        yticklabels(nb_trials : -10 : 0);
        axis tight
        xlim([-2 2]);
        
        % Plot rectified EMG
        h = h/2; %because no more negative values
        
        subplot(1,3,2); hold on
        
        for itrial = 1 : nb_trials
            plot(t,abs(data_EMG{ipatient}.(markername).trial{itrial}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        xlabel(sprintf('Time from %s (s)', markername),'Interpreter','none');
        title(sprintf('Timecourse of \nrectified %s (%d trials)', data_EMG{ipatient}.(markername).label{1}, size(data_EMG{ipatient}.(markername).trial,2)),'Interpreter','none','Fontsize',15);
        setfig();
        tick = h;
        yticks(tick : tick*10 : nb_trials*h);
        yticklabels(nb_trials : -10 : 0);
        set(gca,'TickDir','out');
        axis tight
        xlim([-2 2]);
        
        % Plot envelope of rectified EMG
        subplot(1,3,3);
        hold on
        
        for itrial = 1 : nb_trials
            plot(t,abs(data_EMG{ipatient}.(markername).trial{itrial}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
            plot(t,env{itrial}+ (nb_trials+1)*h - itrial*h,'c','LineWidth',2);
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        
        title(sprintf('Timecourse of envelope \nof rectified %s (%d trials)', data_EMG{ipatient}.(markername).label{1}, size(data_EMG{ipatient}.(markername).trial,2)),'Interpreter','none','Fontsize',15);
        setfig();
        tick = h;
        yticks(tick : tick*10 : nb_trials*h);
        yticklabels(nb_trials : -10 : 0);
        set(gca,'TickDir','out');
        axis tight
        xlim([-2 2]);
        
        %print figure
        fname = fullfile(config{ipatient}.imagesavedir,'emg_method',[config{ipatient}.prefix,convertStringsToChars(markername),'_emg_method_',data_EMG{ipatient}.(markername).label{1}]);
        dtx_savefigure(fig,fname,'pdf','png','close');
    end
end

%% count seizure without artefacts
for ipatient = 1:size(config,2)
    for markername = string(config{ipatient}.LFP.name(1:2))
        idata = "LFP";
        if ~isfield(data{ipatient}.(idata){ipart},markername)
            continue
        end
        if isempty(data{ipatient}.(idata){ipart}.(markername))
            continue
        end
        switch markername
            case 'SlowWave_R'
                count.without_artefacts.n_seizures.right(ipatient)    = size(data{ipatient}.(idata){ipart}.(markername).trial,2);
            case 'SlowWave_L'
                count.without_artefacts.n_seizures.left(ipatient)     = size(data{ipatient}.(idata){ipart}.(markername).trial,2);
        end
    end
    try
        count.without_artefacts.n_seizures.all(ipatient) = count.without_artefacts.n_seizures.right(ipatient) + count.without_artefacts.n_seizures.left(ipatient);
    catch
        count.without_artefacts.n_seizures.all(ipatient) = count.without_artefacts.n_seizures.left(ipatient);
    end
end
count.without_artefacts.n_seizures.total    = sum(count.without_artefacts.n_seizures.all);
count.without_artefacts.n_seizures.med      = median(count.without_artefacts.n_seizures.all);

% plot selected trials for each data : LFP/CSD, SlowWave_R/L, topoplot/slowwavemorpho/propagation
%not done yet. 

%% count eeg duration per patient and save count data


for ipatient = 1:size(config,2)
    ipart = 1;
    for idir = 1:size(config{ipatient}.directorylist{ipart},2)
        [isNeuralynx, isMicromed, isBrainvision] = get_data_format(config{ipatient});
        % find data file
        if isNeuralynx
            temp        = dir(fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir},'*.ncs'));
            fname    = fullfile(config{ipatient}.rawdir,config{ipatient}.directorylist{ipart}{idir},temp.name);
        elseif isMicromed
            fname    = fullfile(config{ipatient}.rawdir, [config{ipatient}.directorylist{ipart}{idir} '.TRC']);
        elseif isBrainvision
            fname    = fullfile(config{ipatient}.rawdir, [config{ipatient}.directorylist{ipart}{idir} '.eeg']);
        end
        fprintf('Reading header from  %s \n',fname);
        hdr = ft_read_header(fname);
        
        count.data_length.dir{ipatient}(idir) = hdr.nSamples / hdr.Fs;
    end
    count.data_length.all(ipatient) = sum(count.data_length.dir{ipatient});
end

count.data_length.total = sum(data_length.all);
count.data_length.min = min(data_length.all);
count.data_length.med = median(data_length.all);
count.data_length.max = max(data_length.all);

%% count emg without artefacts
for ipatient = 1:size(config,2)
    has_R = false;
    has_L = false;
    for markername = string(config{ipatient}.LFP.name(3:4))
        if ipatient > size(data_EMG,2)
            continue
        end
        if ~isfield(data_EMG{ipatient},markername)
            continue
        end
        if isempty(data_EMG{ipatient}.(markername))
            continue
        end
        switch markername
            case 'SlowWave_R_EMG_align'
                has_R = true;
                count.without_artefacts.n_emg.right(ipatient)    = size(data_EMG{ipatient}.(markername).trial,2);
            case 'SlowWave_L_EMG_align'
                has_L = true;
                count.without_artefacts.n_emg.left(ipatient)     = size(data_EMG{ipatient}.(markername).trial,2);
        end
    end
    count.without_artefacts.n_emg.all(ipatient) = 0;
    if has_R
        count.without_artefacts.n_emg.all(ipatient) = count.without_artefacts.n_emg.right(ipatient) + count.without_artefacts.n_emg.all(ipatient);
    end
    if has_L
        count.without_artefacts.n_emg.all(ipatient) = count.without_artefacts.n_emg.left(ipatient) + count.without_artefacts.n_emg.all(ipatient);
    end
end
count.without_artefacts.n_emg.total    = sum(count.without_artefacts.n_emg.all);
count.without_artefacts.n_emg.med      = median(count.without_artefacts.n_emg.all);
count.without_artefacts.n_emg.med_non_zero      = median(count.without_artefacts.n_emg.all(count.without_artefacts.n_emg.all>0));

%% save count variable
%computed in several separated parts of this script :
% - count.data_length
% - count.n_seizures
% - count.n_emg
% - count.without_artefacts.n_seizures
% - count.without_artefacts.n_emg
save(fullfile(config{ipatient}.datasavedir, 'allpatients_data_length.mat'), 'count'); 



%% STOP ICI
%% Old : propagation : abandonné

toremove.SlowWave_R = [2 3 6];
toremove.SlowWave_L = [];

for idata = ["LFP", "CSD"]
    for markername = string(config{ipatient}.LFP.name)
        for ipatient = 1:size(config,2)
            if ismember(ipatient, toremove.(markername))
                continue
            end
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            data_xcorr        = data{ipatient}.(idata){1}.(markername);
            cfgtemp           = [];
            cfgtemp.lpfilter  = 'yes';
            cfgtemp.lpfreq    = 5;
            data_xcorr        = ft_preprocessing(cfgtemp, data_xcorr);
            cfgtemp           = [];
            cfgtemp.channel   = config{ipatient}.morpho.channel.(idata).(markername);
            cfgtemp.latency   = config{ipatient}.morpho.toi.(markername);
            data_xcorr        = ft_selectdata(cfgtemp, data_xcorr);
            data_xcorr        = ft_timelockanalysis([], data_xcorr);
            
            for chan_name = string(data_xcorr.label')
                ichan = strcmp(chan_name, data_xcorr.label);
                [p, l] = findpeaks(-data_xcorr.avg(ichan,:), data_xcorr.time,'NPeaks',1,'SortStr','descend');
                if isempty(p)
                    peak.(idata).(markername){ipatient}.(chan_name) = nan;
                    lag.(idata).(markername){ipatient}.(chan_name) = nan;
                else
                    peak.(idata).(markername){ipatient}.(chan_name) = p;
                    lag.(idata).(markername){ipatient}.(chan_name) = l;
                end
            end
            
            
            %plot xcorr. remove channel if correlation is too small
            [data_xcorr_aligned, nshifts{ipatient}.(idata).(markername), shifted_corr{ipatient}.(idata).(markername)] = alignchannels_Xcorr(data_xcorr,100);
            
            %fig = figure;
            %leg = plot(data_xcorr_aligned.time,data_xcorr_aligned.avg);
            %figure;hold
            %plot(data_xcorr.time,data_xcorr.avg,'k');
            %plot(data_xcorr.time,nanmean(data_xcorr.avg,1),'r','LineWidth',2);
        end
    end
end

% refchan.SlowWave_L = {'Fp1','F3','F7','Fz','C3','Cz'};
refchan.SlowWave_L = {'Cz','Fz','F7','F3','Fp1','C3'};
% refchan.SlowWave_R = {'Fp2','F4','F8','Fz','C4','Cz'};
refchan.SlowWave_R = {'Fz','F4','Cz','Fp2','F8','C4'};


%idata = "LFP"; markername = "SlowWave_R";

for idata = ["LFP", "CSD"]
    for markername = string(config{ipatient}.LFP.name)
        
        clear pat_list leg x timeshifts
        fig = figure; hold;
        count_plot = 0;
        for ipatient = 1:size(config,2)
            if ismember(ipatient, toremove.(markername))
                continue
            end
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            
            %plot xcorr results
            timeshifts{ipatient} = struct2array(nshifts{ipatient}.(idata).(markername))./ config{ipatient}.LFP.resamplefs;
            %re-order the channel to have the same order for each patient
            [x{ipatient}, idx] = match_str(refchan.(markername),fieldnames(nshifts{ipatient}.(idata).(markername))');
%             %plot findpeaks results
%             timeshifts{ipatient} = struct2array(lag.(idata).(markername){ipatient});
%             %re-order the channel to have the same order for each patient
%             [x{ipatient}, idx] = match_str(refchan.(markername),fieldnames(lag.(idata).(markername){ipatient})');

            timeshifts{ipatient} = timeshifts{ipatient}(idx');
            x{ipatient} = x{ipatient}';
            %test = fieldnames(nshifts{ipatient}.(idata).(markername))';
            %test = test(idx');
            
            count_plot = count_plot+1;
            scatter(x{ipatient}+rand*0.2-0.1,timeshifts{ipatient},'ok','filled');
            pat_list{count_plot} = config{ipatient}.prefix(1:end-1);
            
        end
        xlim([0 7]);
        xticks([1:6]);
        xticklabels(refchan.(markername));
        
        %gather patient's data
        for ipatient = 1:size(config,2)
            if ismember(ipatient, toremove.(markername))
                continue
            end
            if ~isfield(data{ipatient}.(idata){ipart},markername)
                continue
            end
            if isempty(data{ipatient}.(idata){ipart}.(markername))
                continue
            end
            propa.(idata).(markername).ipatient = ipatient;
            for ichan = 1:6
                channame = refchan.(markername){ichan};
                if ismember(ichan, x{ipatient})
                    idx = find(x{ipatient}==ichan);
                    propa.(idata).(markername).timeshifts.(channame)(ipatient) = timeshifts{ipatient}(idx);
                else
                    propa.(idata).(markername).timeshifts.(channame)(ipatient) = nan;
                end
            end
            for ichan = 1:6
                channame = refchan.(markername){ichan};
                propa.(idata).(markername).meandata(ichan) = nanmean(propa.(idata).(markername).timeshifts.(channame));
                propa.(idata).(markername).stddata(ichan)  = nanstd(propa.(idata).(markername).timeshifts.(channame));
            end
        end
                
        errorbar(1:6, propa.(idata).(markername).meandata, propa.(idata).(markername).stddata,'--xr', 'LineWidth', 2,'CapSize',10,'MarkerSize',10);
        
    end
end

%% OLD LFP xcorr


for idata = ["LFP", "CSD"]
    for markername = string(config{ipatient}.LFP.name)
            for ipatient = 1:size(config,2)
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
                cfgtemp.LFP.xcorr.toi       = 'all';%FIXME : voir si ne pas sï¿½lectionner, et choisir une pï¿½riode par dï¿½faut
                cfgtemp.LFP.xcorr.plotdata  = 'yes';
                data_xcorr{ipatient}.(idata).(markername)        = dtx_xcorr_LFP(cfgtemp,data_temp_avg);
                
%                 % plot avg data
%                 fig = figure;
%                 plot(data_temp_avg.time,data_temp_avg.avg);
%                 legend(data_temp_avg.label');
%                 
%                 %print to file
%                 fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('%s_lfp_for_xcorr_%s_%s.pdf',config{ipatient}.prefix,markername,idata));
%                 dtx_savefigure(fig,fname,'pdf','png','close');
                
            end
        end
    end
end

save(fullfile(config{ipatient}.datasavedir,'allpatients_xcorr_lfp.mat'), 'data_xcorr', '-v7.3');
load(fullfile(config{ipatient}.datasavedir,'allpatients_xcorr_lfp.mat'), 'data_xcorr');

%plot distrib des lags pour chaque ï¿½lectrode (1 point = 1 patient)
latency_xcorr_all = [];
for idata = ["LFP", "CSD"]
    for markername = string(config{ipatient}.LFP.name)
        
        %gather data from all patients
        for ipatient = 1:size(config,2)
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
        
        fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_xcorr_lfp_%s_%s',markername,idata));
        dtx_savefigure(fig,fname,'pdf','png','close');
     
    end%markername
end%idata







