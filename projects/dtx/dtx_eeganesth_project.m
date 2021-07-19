% DTX sedated rats, EEG analysis

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\development
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/development
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

config = dtx_eeganesth_setparams;
ipart = 1;

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\nRead marker data for %s\n',config{ipatient}.prefix(1:end-1));
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient},false);
    MuseStruct_concat{ipatient} = concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},false);
end

%% compute seizures descriptive stats
for ipatient = 1:size(config,2)
    config{ipatient}.seizuretimings.injection_clock= config{ipatient}.injectiontime;
    seizure_timing{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
    count.seizures.all(ipatient)  = size(seizure_timing{ipatient}.time_start.synctime,2);
    count.slowwaves.all(ipatient) = size(MuseStruct_concat{ipatient}{1}.markers.SlowWave.synctime,2);
    count.data_length.all(ipatient) = hours(seizure_timing{ipatient}.recordduration);
end
count.seizures.total   = sum(count.seizures.all);
count.seizures.min     = min(count.seizures.all);
count.seizures.med     = median(count.seizures.all);
count.seizures.max     = max(count.seizures.all);
count.slowwaves.total  = sum(count.slowwaves.all);
count.slowwaves.min    = min(count.slowwaves.all);
count.slowwaves.med    = median(count.slowwaves.all);
count.slowwaves.max    = max(count.slowwaves.all);
count.data_length.time_unit = 'hours';
count.data_length.total = sum(count.data_length.all);
count.data_length.min = min(count.data_length.all);
count.data_length.med = median(count.data_length.all);
count.data_length.max = max(count.data_length.all);

%% plot number of seizures over time
iparam = "nb_seizures";
clear data
fig = figure;
sgtitle(iparam,'Interpreter', 'none','Fontsize', 18, 'FontWeight', 'bold');
subplot(2,1,1);hold on
data.label = {'dummy'};
for ipatient = 1:size(config,2)
    data.time{ipatient} = hours(seizure_timing{ipatient}.statsovertime.endtime);
    data.trial{ipatient}= seizure_timing{ipatient}.statsovertime.(iparam);
    plot(data.time{ipatient},data.trial{ipatient}, 'k');
    idxstart = find(~isnan(data.trial{ipatient}),1,'first');
    s = plot(data.time{ipatient}(idxstart),data.trial{ipatient}(idxstart),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
    s.ZData = ones(size(s.YData)).*2;
    idxend = find(~isnan(data.trial{ipatient}),1,'last');
    s = plot(data.time{ipatient}(idxend),data.trial{ipatient}(idxend),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
    s.ZData = ones(size(s.YData)).*2;
end
data_avg = ft_timelockanalysis([],data);
std_data = sqrt(data_avg.var);
x = data_avg.time;
y = [data_avg.avg - std_data; std_data; std_data]';
filled_SD = area(x,y);
filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
filled_SD(1).ShowBaseLine = 'off';

p = plot(data_avg.time, data_avg.avg, 'r', 'LineWidth', 2);
p.ZData = ones(size(p.YData));
axis tight
ylabel('nb. per hour');
xlabel('time (hours post injection)');
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);

%plot smoothed avg and std
subplot(2,1,2);hold on
avg_data_smooth = movmean(data_avg.avg,1,'omitnan');
std_data_smooth = movmean(std_data,1,'omitnan');
x = data_avg.time;
y = [avg_data_smooth - std_data_smooth; std_data_smooth; std_data_smooth]';
filled_SD = area(x,y);
filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
filled_SD(2).FaceColor = 'k'; filled_SD(3).FaceColor = 'k';
filled_SD(1).ShowBaseLine = 'off';
p = plot(data_avg.time, avg_data_smooth, 'Color','k', 'LineWidth', 2);
p.ZData = ones(size(p.YData));
axis tight
ax = axis;
xlim([0 ax(2)]);
ylabel('nb. per hour');
xlabel('time (hours post injection)');
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
axis tight;
xlim([0 11.2]);

fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_over_time',iparam));
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%% plot seizure frequency, cv2 and duration
for iparam = ["timebetween2seizures", "timebetween2seizures_cv2", "seizureduration"]
    fig=figure; hold on;
    for ipatient = 1:size(config,2)
        if contains(iparam, "timebetween2seizures")
            if ~contains(iparam, 'cv2')
                to_plot = minutes(seizure_timing{ipatient}.(iparam));
            else
                to_plot = seizure_timing{ipatient}.(iparam);
            end
        elseif iparam == "seizureduration"
            to_plot = seconds(seizure_timing{ipatient}.(iparam));
        end
        bar(ipatient, mean(to_plot, 'omitnan'),'FaceColor','white','EdgeColor','k', 'LineWidth', 2);
        scatter(rand(1,size(to_plot,2))*0.2+ipatient-0.1,to_plot, 'o', 'MarkerEdgeColor', 'k');
        errorbar(ipatient, mean(to_plot, 'omitnan'), 0,std(to_plot, 'omitnan'),'k','CapSize',10, 'LineWidth', 2);
        datamean.(iparam)(ipatient) = mean(to_plot, 'omitnan');
    end
    set(gca, 'LineWidth', 2);
    title(iparam,'Interpreter','none');
    ax = axis;
    ylim([0 ax(4)]);
    xticks(1:size(config,2));
    xlabel('rat nr');
    set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_distrib',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end

%% #DATA EEG, computed in dtx_eegrodents_cluster.m (load precomputed data)
for ipatient = 1:size(config,2)
    fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
    temp = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
    LFP{ipatient} = temp.LFP;
    clear temp
    LFP{ipatient} = removeArtefactedTrials(config{ipatient}, LFP{ipatient});
end

%% TDW : overdraw of averages from all rats
markername = "SlowWave";
toi = 2;
fig = figure; hold on;
for ipatient = 1:size(config,2)
    cfgtemp = [];
    cfgtemp.lpfilter = 'yes';
    cfgtemp.lpfreq = 30;
    data_plot = ft_preprocessing(cfgtemp, LFP{ipatient}{1}.(markername));
    
    cfgtemp         = [];
    cfgtemp.channel = config{ipatient}.align.channel.(markername);
    data_plot       = ft_timelockanalysis(cfgtemp, data_plot);
    plot(data_plot.time, normalize(data_plot.avg, 'range'), color);
end
xlim([-toi toi]);
xlabel('time (s)');
ylabel('uV');
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
figname = fullfile(config{ipatient}.imagesavedir, '..','average_M1G', sprintf('allrats_average_M1G_%s_%d', markername, toi));
dtx_savefigure(fig, figname, 'png', 'pdf', 'close');

%% TDW : overdraw of all trials, per rat
for ipatient = 1:size(config,2)
    cfgtemp         = [];
    cfgtemp.channel = config{ipatient}.align.channel.SlowWave;
    to_plot         = ft_selectdata(cfgtemp, LFP{ipatient}{1}.SlowWave);
    avg = ft_timelockanalysis([], to_plot);
    
    fig = figure; hold on;
    for itrial = 1:size(to_plot.trial, 2)
        y = to_plot.trial{itrial}/10;
        p = plot(to_plot.time{itrial}, y,'k');
    end
    plot(avg.time, avg.avg/10, 'white', 'LineWidth', 2);
    xlim([-2 2]);
    set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
    xlabel('Time (s)');
    ylabel('uV');
    fname = fullfile(config{ipatient}.imagesavedir, '..', 'slowwave_trials', sprintf('%sSlowWave', config{ipatient}.prefix));
    dtx_savefigure(fig, fname, 'png', 'pdf', 'close');
end

%% #DATA computed in dtx_cluster_eegvideo.m# slow wave morphology
clear morpho
for ipatient = 1:size(config, 2)
    temp = load(fullfile(config{ipatient}.datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)), 'morpho');
    morpho{ipatient} = temp.morpho;
end

%% TDW halfwidth
clear param meanparam stdparam
iparam = "halfwidth";
fig=figure;hold on
for ipatient = 1:size(config, 2)
    if isempty(morpho{ipatient})
        continue
    end
    param = abs(morpho{ipatient}.(iparam));
    meanparam.(iparam)(ipatient) = mean(param, 'omitnan');
    stdparam.(iparam)(ipatient)  = std(param, 'omitnan');
    bar(ipatient, meanparam.(iparam)(ipatient),'FaceColor','white','EdgeColor','k', 'LineWidth', 2);
    scatter(rand(size(param))*0.4+ipatient-0.2, param, 'o', 'MarkerEdgeColor', 'k');
    errorbar(ipatient, meanparam.(iparam)(ipatient),0, stdparam.(iparam)(ipatient),'k', 'LineWidth', 2);
    datamean.(iparam)(ipatient) = meanparam.(iparam)(ipatient);
end
set(gca, 'LineWidth', 2);
ax = axis;
ylim([0 ax(4)]);
xlim([0, size(config,2)+1]);
xticks(1:size(config,2));
xticklabels(1:size(config,2));
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
ylabel(iparam);
set(gca, 'Fontsize', 11);
fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_morpho_%s',iparam));
dtx_savefigure(fig,fname,'pdf','png','close');

%% count seizures without artefacts
for ipatient = 1:size(config,2)
    count.without_artefacts.seizures.all(ipatient)= size(LFP{ipatient}{ipart}.SlowWave.trial,2);
end
count.without_artefacts.seizures.total    = sum(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.min      = min(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.med      = median(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.max      = max(count.without_artefacts.seizures.all);

%% save count variable
save(fullfile(config{ipatient}.datasavedir, 'allpatients_count_datalength.mat'), 'count'); 