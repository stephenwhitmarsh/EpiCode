% DTX EEG video rodents analysis

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

config = dtx_eegvideo_setparams;
ipart = 1;

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\n Read marker data for %s \n',config{ipatient}.prefix(1:end-1));
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient},false);
    MuseStruct_concat{ipatient} = concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},false);  
end

%% compute seizures and emg descriptive stats
for ipatient = 1:size(config,2)
    config{ipatient}.seizuretimings.injection_clock = config{ipatient}.injectiontime;
    seizure_timing{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
    config{ipatient}.emgtimings.injection_clock= config{ipatient}.injectiontime;
    emg_timing{ipatient} = dtx_stats_emg_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
end

%% count nr of seizures 
for ipatient = 1:size(config,2)
    count.seizures.all(ipatient)  = sum(seizure_timing{ipatient}.time_start.clock - config{ipatient}.injectiontime < hours(10));
    count.slowwaves.all(ipatient) = sum(MuseStruct_concat{ipatient}{1}.markers.SlowWave.clock  - config{ipatient}.injectiontime < hours(10));
    count.emg.all(ipatient)       = sum(emg_timing{ipatient}.time_start.clock - config{ipatient}.injectiontime < hours(10));
    count.data_length.t_start(ipatient) = hours(MuseStruct_concat{ipatient}{1}.markers.Analysis_Start.clock - config{ipatient}.injectiontime);
    sel = find((MuseStruct_concat{ipatient}{1}.markers.Crise_Start.clock - config{ipatient}.injectiontime) < hours(10));
    count.data_length.t_end(ipatient) = hours(MuseStruct_concat{ipatient}{1}.markers.Crise_Start.clock(sel(end)) - config{ipatient}.injectiontime);
    count.data_length.all(ipatient) = count.data_length.t_end(ipatient) - count.data_length.t_start(ipatient);
end
count.seizures.total        = sum(count.seizures.all);
count.seizures.min          = min(count.seizures.all);
count.seizures.med          = median(count.seizures.all);
count.seizures.max          = max(count.seizures.all);
count.slowwaves.total       = sum(count.slowwaves.all);
count.slowwaves.min         = min(count.slowwaves.all);
count.slowwaves.med         = median(count.slowwaves.all);
count.slowwaves.max         = max(count.slowwaves.all);
count.emg.total             = sum(count.emg.all);
count.emg.min               = min(count.emg.all);
count.emg.med               = median(count.emg.all);
count.emg.max               = max(count.emg.all);
count.emg.min_non_zero      = min(count.emg.all(count.emg.all>0));
count.emg.med_non_zero      = median(count.emg.all(count.emg.all>0));
count.data_length.time_unit = 'hours';
count.data_length.total     = sum(count.data_length.all);
count.data_length.min       = min(count.data_length.all);
count.data_length.med       = median(count.data_length.all);
count.data_length.max       = max(count.data_length.all);

%% plot number of seizures over time
iparam = "nb_seizures";
clear data
fig = figure;
sgtitle(iparam,'Interpreter', 'none','Fontsize', 18, 'FontWeight', 'bold');
subplot(2,1,1);hold on
data.label = {'dummy'};
for ipatient = 1:size(config,2)
    data.time{ipatient} = hours(seizure_timing{ipatient}.statsovertime.starttime + seizure_timing{ipatient}.statsovertime.endtime) ./ 2; %middle of the window
    data.trial{ipatient}= seizure_timing{ipatient}.statsovertime.(iparam);
    p = plot(data.time{ipatient},data.trial{ipatient}, 'k');
    p.ZData = ones(size(p.YData));
    idxstart = find(~isnan(data.trial{ipatient}),1,'first');
    s = plot(data.time{ipatient}(idxstart),data.trial{ipatient}(idxstart),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
    s.ZData = ones(size(s.YData)).*2;
    idxend = find(~isnan(data.trial{ipatient}),1,'last');
    s = plot(data.time{ipatient}(idxend),data.trial{ipatient}(idxend),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
    s.ZData = ones(size(s.YData)).*2;
end

%plot avg and std
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

%set figure display
axis tight
ax = axis;
xlim([0 ax(2)]);
ylabel('nb. per hour');
xlabel('time (hours post injection)');
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);

%plot smoothed avg and std
subplot(2,1,2);hold on
avg_data_smooth = movmean(data_avg.avg(3:end),1);
std_data_smooth = movmean(std_data(3:end),1);

x = data_avg.time(3:end);
y = [avg_data_smooth - std_data_smooth; std_data_smooth; std_data_smooth]';
filled_SD = area(x,y);
filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
filled_SD(2).FaceColor = 'k'; filled_SD(3).FaceColor = 'k';
filled_SD(1).ShowBaseLine = 'off';

p = plot(x, avg_data_smooth, 'Color','k', 'LineWidth', 2);
p.ZData = ones(size(p.YData));

%set figure display
axis tight
ax = axis;
xlim([0 ax(2)]);
ylabel('nb. per hour');
xlabel('time (hours post injection)');
set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);

fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_over_time',iparam));
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%% plot seizure frequency, cv2 and duration
for iparam = ["timebetween2seizures", "timebetween2seizures_cv2", "seizureduration"]
    for ipatient = 1:size(config,2)
        if contains(iparam, "timebetween2seizures")
            sel = seizure_timing{ipatient}.x_timebetween2seizures - config{ipatient}.injectiontime < hours(10); %select the 10 first hours
            if ~contains(iparam, 'cv2')
                to_plot{ipatient} = minutes(seizure_timing{ipatient}.(iparam)(sel));
            else
                to_plot{ipatient} = seizure_timing{ipatient}.(iparam)(sel(1:end-1));
            end
        else
            sel = seizure_timing{ipatient}.x_timebetween2seizures - config{ipatient}.injectiontime < hours(10); %select the 10 first hours
            if ~contains(iparam, 'cv2')
                to_plot{ipatient} = seconds(seizure_timing{ipatient}.(iparam)(sel));
            else
                to_plot{ipatient} = seizure_timing{ipatient}.(iparam)(sel(1:end-1));
            end
        end
    end
    y = [to_plot{:}]';
    g = [];
    for i = 1:length(to_plot)
        g = [g; repmat(i, length(to_plot{i}), 1)];
    end
    g(isnan(y)) = [];
    y(isnan(y)) = [];
    fig = figure;hold on;
    boxplot(y, g, 'symbol', 'ok');
    set(gca, 'tickdir', 'out', 'fontsize', 25);
    set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
    if strcmp(iparam, "timebetween2seizures")
        ylim([-10 45]);
    end
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_boxplot',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig');
    
    %table
    h = findobj(fig,'Type','line');
    tableout = table.empty;
    for ifield = unique(string({h.Tag}))
        if ifield == "Outliers"
            continue
        end
        ipat = length(config) + 1;
        for i = find(strcmp({h.Tag}, ifield))
            ipat = ipat-1;
            if ifield == "Box"
                tableout.(sprintf('%s_25', ifield))(ipat) = min(h(i).YData);
                tableout.(sprintf('%s_75', ifield))(ipat) = max(h(i).YData);
            elseif ifield == "Lower Whisker"
                tableout.(ifield)(ipat) = min(h(i).YData);
            elseif ifield == "Upper Whisker"
                tableout.(ifield)(ipat) = max(h(i).YData);
            else
                tableout.(ifield)(ipat) = unique(h(i).YData);
            end
        end
    end
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_table.csv',iparam));
    writetable(tableout, fname, 'delimiter', ';');
    
    %average
    y = [];
    for ipatient = 1:size(to_plot, 2)
        y = [y, mean(to_plot{ipatient}, 'omitnan')];
    end
    fig = figure; hold on;
    scatter(rand(size(y))*0.1-0.05+1, y, 50, 'ok', 'filled');
    boxplot(y, 'symbol', 'k');
    set(gca, 'tickdir', 'out', 'fontsize', 25);
    set(findall(gca, 'type', 'line'), 'linewidth', 2);
    if contains(iparam, 'cv2')
        ylabel('CV2');
        x = xlim;
        ylim([0 2]);
        plot(x, [1 1], '--k');
    end
    
    if strcmp(iparam, "timebetween2seizures")
        ylim([-10 45]);
    end
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_boxplot_average',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end

%% plot emg duration, eeg-emg delay

%emg duration
to_plot = cellfun(@(c) c.emg_duration, emg_timing, 'uniformoutput', false);
y = [to_plot{:}]';
g = [];
for i = 1:length(to_plot)
    g = [g; repmat(i, length(to_plot{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
fig = figure;hold on;
boxplot(y, g, 'symbol', 'ok');
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
ylim([-1 3]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description','allpatients_emg_duration_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig', 'close');

%table
h = findobj(fig,'Type','line');
tableout = table.empty;
for ifield = unique(string({h.Tag}))
    if ifield == "Outliers"
        continue
    end
    ipat = length(config) + 1;
    for i = find(strcmp({h.Tag}, ifield))
        ipat = ipat-1;
        if ifield == "Box"
            tableout.(sprintf('%s_25', ifield))(ipat) = min(h(i).YData);
            tableout.(sprintf('%s_75', ifield))(ipat) = max(h(i).YData);
        elseif ifield == "Lower Whisker"
            tableout.(ifield)(ipat) = min(h(i).YData);
        elseif ifield == "Upper Whisker"
            tableout.(ifield)(ipat) = max(h(i).YData);
        else
            tableout.(ifield)(ipat) = unique(h(i).YData);
        end
    end
end
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description','allpatients_emg_duration_table.csv');
writetable(tableout, fname, 'delimiter', ';');

y = [];
for ipatient = 1:size(to_plot, 2)
    y = [y, mean(to_plot{ipatient}, 'omitnan')];
end
fig = figure; hold on;
scatter(rand(size(y))*0.1-0.05+1, y, 50, 'ok', 'filled');
boxplot(y, 'symbol', 'k');
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylim([-1 3]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description','allpatients_emg_duration_boxplot_average');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%eeg-emg delay
to_plot = cellfun(@(c) c.eeg_emg_delay, emg_timing, 'uniformoutput', false);
y = [to_plot{:}]';
g = [];
for i = 1:length(to_plot)
    g = [g; repmat(i, length(to_plot{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
fig = figure;hold on;
boxplot(y, g, 'symbol', 'ok');
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description','allpatients_eeg_emg_delay_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig', 'close');

%table
h = findobj(fig,'Type','line');
tableout = table.empty;

for ifield = unique(string({h.Tag}))
    if ifield == "Outliers"
        continue
    end
    ipat = length(config) + 1;
    for i = find(strcmp({h.Tag}, ifield))
        ipat = ipat-1;
        if ifield == "Box"
            tableout.(sprintf('%s_25', ifield))(ipat) = min(h(i).YData);
            tableout.(sprintf('%s_75', ifield))(ipat) = max(h(i).YData);
        elseif ifield == "Lower Whisker"
            tableout.(ifield)(ipat) = min(h(i).YData);
        elseif ifield == "Upper Whisker"
            tableout.(ifield)(ipat) = max(h(i).YData);
        else
            tableout.(ifield)(ipat) = unique(h(i).YData);
        end
    end
end
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description','allpatients_eeg_emg_delay_table.csv');
writetable(tableout, fname, 'delimiter', ';');
    
%% #DATA EEG, computed in dtx_eegrodents_cluster.m (load precomputed data)
for ipatient = 1:size(config,2)
    fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
    lfptemp = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
    lfptemp.LFP = removeArtefactedTrials(config{ipatient}, lfptemp.LFP);
    LFP{ipatient} = lfptemp.LFP;
    clear lfptemp
end

%% TDW : overdraw of averages from all rats
fig = figure; hold on;
markername = "SlowWave";
for ipatient = 1:size(config, 2)
    cfgtemp         = [];
    cfgtemp.trials  = LFP{ipatient}{1}.(markername).trialinfo.starttime - config{ipatient}.injectiontime < hours(10);
    cfgtemp.channel = config{ipatient}.align.channel.(markername);
    to_plot         = ft_selectdata(cfgtemp, LFP{ipatient}{1}.(markername));
    
    cfgtemp         = [];
    cfgtemp.lpfilter = 'yes';
    cfgtemp.lpfreq  = 100;
    to_plot         = ft_preprocessing(cfgtemp, to_plot);
    
    y   = cat(1,to_plot.trial{:});
    y   = mean(y, 1, 'omitnan');
    x   = to_plot.time{1};
    y   = y - mean(y(x>-2&x<-1));
    [ymin, loc] = min(y(x>-0.3&x<0.2));
    y   = y ./ -ymin;
    temp = x(x>-0.3&x<0.2);
    x   = x - temp(loc);
    
    plot(x, y, 'LineWidth', 2);
end

axis tight;
xlim([-2 2]);
set(gca, 'TickDir', 'out', 'FontSize', 25, 'FontWeight', 'bold');
xlabel('Time (s)');
ylabel('uV');
fname = fullfile(config{ipatient}.imagesavedir, '..', 'figure_morpho', 'Allpatients-SlowWave_morpho_averages');
dtx_savefigure(fig,fname, 'png', 'pdf', 'close');

%% TDW : overdraw of all trials, per rat
for ipatient = 1:size(config, 2)
    cfgtemp         = [];
    cfgtemp.channel = 'M1G';
    cfgtemp.trials  = LFP{ipatient}{1}.SlowWave.trialinfo.starttime - config{ipatient}.injectiontime < hours(10);
    to_plot         = ft_selectdata(cfgtemp, LFP{ipatient}{1}.SlowWave);
    
    cfgtemp                 = [];
    cfgtemp.demean          = 'yes';
    cfgtemp.baselinewindow  = [-2 -1];
    to_plot                 = ft_preprocessing(cfgtemp, to_plot);
    
    avg = ft_timelockanalysis([], to_plot);
    
    fig = figure; hold on;
    for itrial = 1:size(to_plot.trial, 2)
        y = to_plot.trial{itrial}/10;
        p = plot(to_plot.time{itrial}, y, 'k');
    end
    
    plot(avg.time, avg.avg/10, 'Color', [0.2 0.6 1], 'LineWidth', 2);
    
    xlim([-2 2]);
    set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
    xlabel('Time (s)');
    ylabel('uV');
    fname = fullfile(config{ipatient}.imagesavedir, sprintf('%sSlowWave_Figure', config{ipatient}.prefix));
    dtx_savefigure(fig, fname, 'png', 'close');
end

%% #DATA computed in dtx_cluster_eegvideo.m# slow wave morphology
clear morpho
for ipatient = 1:size(config,2)
    temp = load(fullfile(config{ipatient}.datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)), 'morpho');    
    morpho{ipatient} = temp.morpho;
end

%% TDW half width
iparam = "halfwidth";
y = [];
for ipatient = 1:size(config,2)
    if isempty(morpho{ipatient})
        continue
    end
    sel = morpho{ipatient}.time - config{ipatient}.injectiontime < hours(10);
    param{ipatient} = abs(morpho{ipatient}.(iparam)(sel));
end
y = [param{:}]';
g = [];
for i = 1:length(param)
    g = [g; repmat(i, length(param{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
fig = figure;hold on;
boxplot(y, g, 'symbol', 'ok');
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
ylim([-0.1 1.5]);
ylabel('ms');
yticklabels(yticks.*1000);

%print to file
fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_morpho_%s_boxplot',iparam));
dtx_savefigure(fig,fname,'pdf','png','close');

%table
h = findobj(fig,'Type','line');
tableout = table.empty;
for ifield = unique(string({h.Tag}))
    if ifield == "Outliers"
        continue
    end
    ipat = length(config) + 1;
    for i = find(strcmp({h.Tag}, ifield))
        ipat = ipat-1;
        if ifield == "Box"
            tableout.(sprintf('%s_25', ifield))(ipat) = min(h(i).YData);
            tableout.(sprintf('%s_75', ifield))(ipat) = max(h(i).YData);
        elseif ifield == "Lower Whisker"
            tableout.(ifield)(ipat) = min(h(i).YData);
        elseif ifield == "Upper Whisker"
            tableout.(ifield)(ipat) = max(h(i).YData);
        else
            tableout.(ifield)(ipat) = unique(h(i).YData);
        end
    end
end
fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_morpho_%s_table.csv',iparam));
writetable(tableout, fname, 'delimiter', ';');

%average
y = [];
for ipatient = 1:size(param, 2)
    y = [y, mean(param{ipatient}, 'omitnan')];
end
fig = figure; hold on;
scatter(rand(size(y))*0.1-0.05+1, y, 50, 'ok', 'filled');
boxplot(y);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylim([-0.1 1.5]);
ylabel('ms');
yticklabels(yticks.*1000);
fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_morpho_%s_boxplot_average',iparam));
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%% #DATA EMG : clean ekg and compute envelopes (load precomputed data)
if false
    clear env*
    for markername = ["SlowWave_EMG_begin", "SlowWave"]
        for ipatient = 1:size(config,2)
            if ~isfield(LFP{ipatient}{ipart}, "SlowWave_EMG_begin")
                continue
            end
            if isempty(LFP{ipatient}{ipart}.SlowWave_EMG_begin)
                continue
            end
            if isempty(LFP{ipatient}{ipart}.SlowWave_EMG_begin.label)
                continue
            end
            
            %find non-artefacted EMG trials
            if strcmp(markername, "SlowWave_EMG_begin")
                triallist = 1:size(LFP{ipatient}{ipart}.(markername).trial, 2);
            elseif strcmp(markername, "SlowWave")
                triallist = [];
                for itrial = 1:size(LFP{ipatient}{ipart}.(markername).trial, 2)
                    starttime = LFP{ipatient}{ipart}.SlowWave.trialinfo.starttime(itrial);
                    difftime = seconds(min(abs(LFP{ipatient}{ipart}.SlowWave_EMG_begin.trialinfo.starttime - starttime)));
                    if difftime < 5
                        triallist = [triallist itrial];
                    end
                end
            end
            
            cfgtemp         = [];
            cfgtemp.channel = 'EMG1';
            EMG             = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
            
            cfgtemp         = [];
            cfgtemp.channel = 'EMG2';
            ref             = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
            
            if isempty(EMG.label) %some rats do not have EMG 
                continue
            end

            %separate ekg artefacts from the EMG
            cfgtemp         = [];
            cfgtemp.channel = {'EMG1', 'EMG2'};
            emg_bipolar     = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(markername));
            emg_ica         = ft_componentanalysis([], emg_bipolar);
            
            %compute envelopes
            i = 0;
            tokeep = true(size(triallist));
            t = EMG.time{1};
            for itrial = triallist
                emg_cleaned = emg_ica.trial{itrial}(2,:);
                rect_emg = abs(emg_cleaned);
                e = envelope(rect_emg, config{ipatient}.EMG.envparam, config{ipatient}.EMG.envmethod);
                if max(e(t>-2&t<-0.5)) > 200 %remove EMG artefacted before TDW
                    tokeep(itrial) = false;
                    continue
                end
                i = i+1;
                env.(markername){ipatient}.trial{i} = e;
                env.(markername){ipatient}.time{i}  = t;
            end
            env.(markername){ipatient}.label = {'dummy'};
            env.(markername){ipatient}.trialinfo = EMG.trialinfo(triallist(tokeep),:);
            env_avg.(markername){ipatient} = ft_timelockanalysis([], env.(markername){ipatient});
            
            %value to do normalization later
            bl_avg.(markername)(ipatient) = mean(env_avg.(markername){ipatient}.avg(t>-2&t<-0.5), 'omitnan');
            norm_factor.(markername)(ipatient) = max(env_avg.(markername){ipatient}.avg - bl_avg.(markername)(ipatient));
            
            cfgtemp = [];
            cfgtemp.trials = triallist(tokeep);
            eeg.(markername){ipatient} = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(markername));
            
        end
    end
    fname = fullfile(config{1}.datasavedir, 'allrats_emg_envelopes.mat');
    save(fname, 'env', 'env_avg', 'bl_avg', 'norm_factor', 'eeg', '-v7.3');
else
    fname = fullfile(config{1}.datasavedir, 'allrats_emg_envelopes.mat');
    fprintf('reading %s\n', fname);
    load(fname, 'env', 'env_avg', 'bl_avg', 'norm_factor', 'eeg');
end

%% plot overdraw and avg of EMG envelopes for each rat
for markername = ["SlowWave_EMG_begin", "SlowWave"]
    for ipatient = 1:size(config,2)
        if ipatient > size(env.(markername), 2)
            continue
        end
        if isempty(env.(markername){ipatient})
            continue
        end

        fig = figure; hold on
        count = 0;
        sgtitle(sprintf('%s : %d trials', config{ipatient}.prefix(1:end-1), size(env.(markername){ipatient}.trial,2)),'Interpreter','none','FontSize',18,'FontWeight','bold');
        for itrial = 1:size(env.(markername){ipatient}.trial,2)
            t        = env.(markername){ipatient}.time{itrial};
            bl_trial = mean(env.(markername){ipatient}.trial{itrial}(t>-2&t<-0.5), 'omitnan');
            trial    = env.(markername){ipatient}.trial{itrial} - bl_trial;
            plot(t,trial,'color', [0.5 0.25 0]);
        end
        EMG_avg_blcorrected = env_avg.(markername){ipatient}.avg-bl_avg.(markername)(ipatient);
        p = plot(t,EMG_avg_blcorrected, 'Color', [0.2 0.6 1],'LineWidth',2);

        ax = axis;
        xlim([-2 2]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho', sprintf('%semg_morphology_figure_%s',config{ipatient}.prefix, markername));
        dtx_savefigure(fig,fname, 'png','pdf','close'); %do not save in pdf as it takes 6 hours per figure
    end
end

%% plot EMG and EEG averages, for all rats
for markername = ["SlowWave_EMG_begin", "SlowWave"]
    fig = figure; hold on;
    sgtitle(markername, 'interpreter', 'none');
    % plot eeg
    subplot(2,1,1); hold on; 
    for ipatient = 1:size(config, 2)
        if isempty(LFP{ipatient}{1}.(markername))
            continue
        end
        cfgtemp = [];
        cfgtemp.keeptrials = 'no';
        cfgtemp.channel = 'M1G';
        eegplot = ft_timelockanalysis(cfgtemp, LFP{ipatient}{1}.(markername));
        
        cfgtemp = [];
        cfgtemp.lpfilter = 'yes';
        cfgtemp.lpfreq = 70;
        eegplot = ft_preprocessing(cfgtemp, eegplot);
        
        y = eegplot.avg;
        x = eegplot.time;
        y = y - mean(y(x>-2&x<-1));
        [ymin, loc] = min(y(x>-0.3&x<1));
        y = y ./ -ymin;
        plot(x, y, 'k');
        
        ax = axis;
        xlim([-2 2]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
    end
    
    %plot emg
    for ipatient = 1:size(config,2)
        if ipatient > size(env.(markername), 2)
            continue
        end
        if isempty(env.(markername){ipatient})
            continue
        end
        subplot(2,1,2); hold on;
        for itrial = 1:size(env.(markername){ipatient}.trial,2)
            t = env.(markername){ipatient}.time{itrial};
            bl_trial   = mean(env.(markername){ipatient}.trial{itrial}(t>-2&t<-0.5), 'omitnan');
            trial = env.(markername){ipatient}.trial{itrial} - bl_trial;
        end
        avg_toplot = (env_avg.(markername){ipatient}.avg - bl_avg.(markername)(ipatient));  
        x = env_avg.(markername){ipatient}.time;
        y = avg_toplot ./ max(avg_toplot(x>-0.5&x<0.5));
        p = plot(x, y,'k');
        ax = axis;
        xlim([-2 2]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);
    end
    fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho_align_begin', sprintf('allpatients-emg_morphology_figure_%s', markername));
    dtx_savefigure(fig,fname, 'png','pdf','close');
end

%% TFR : plot one example for the figure 1
ipatient = 3;
itrial = 271;
toi = [-5 15];

fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
temp     = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
temp.LFP = removeArtefactedTrials(config{ipatient}, temp.LFP);
LFP      = temp.LFP;

LFP{1}.Seizure.trialinfo.startfile = LFP{1}.Seizure.trialinfo.begsample / LFP{1}.Seizure.fsample;

cfgtemp         = [];
cfgtemp.trials  = itrial;
cfgtemp.channel = config{ipatient}.LFP.motorcortex;
sel             = ft_selectdata(cfgtemp, LFP{1}.Seizure);

cfgtemp             = [];
cfgtemp.hpfilter    = 'no';
cfgtemp.hpfreq      = 1;
cfgtemp.hpfilttype  = 'fir';
sel_filt            = ft_preprocessing(cfgtemp,sel);

cfgtemp            = [];
cfgtemp.channel    = 'all';
cfgtemp.method     = 'mtmconvol';
cfgtemp.output     = 'pow';
cfgtemp.taper      = 'dpss';
cfgtemp.tapsmofrq  = 2.86;
cfgtemp.pad        = 'nextpow2';
cfgtemp.keeptrials = 'yes';
cfgtemp.foi        = 1:50;
cfgtemp.t_ftimwin  = ones(size(cfgtemp.foi)) .* config{ipatient}.TFR.timewinsize; 
cfgtemp.toi        = '50%';
TFR                = ft_freqanalysis(cfgtemp,sel_filt);

%plot TFR
fig = figure; subplot(4,1,2);
cfgtemp         = [];
cfgtemp.zlim    = 'maxmin';
cfgtemp.figure  = 'gcf';
cfgtemp.xlim    = toi;
cfgtemp.colormat = jet;
cfgtemp.colorbar = 'no';
cfgtemp.baseline = [-5 -2];
cfgtemp.baselinetype = 'relchange';
ft_singleplotTFR(cfgtemp, TFR);
ft_pimpplot(gcf, jet, true);
set(gca, 'XColor', 'white', 'TickDir', 'out');
title([]);
xticks([]);
ylabel('Frequency (Hz)');

caxis([0, 2500]);
c = colorbar;
set(c, 'Location', 'eastoutside', 'color', [0 0 0]);
pos = get(c, 'Position');
pos(1) = 0.92;
set(c, 'pos', pos);

%plot M1G
fig = figure;
cfgtemp         = [];
cfgtemp.trials  = itrial;
cfgtemp.channel = config{ipatient}.LFP.motorcortex;
sel             = ft_selectdata(cfgtemp, LFP{1}.Seizure);
cfgtemp         = [];
cfgtemp.lpfilter = 'yes';
cfgtemp.lpfreq  = 50;
sel = ft_preprocessing(cfgtemp, sel);

plot(sel.time{1}, sel.trial{1}, 'k', 'LineWidth', 0.3);
xlim(toi);
y = ylim;
set(gca, 'XColor', 'white', 'TickDir', 'out');
xticks([]);

%plot M1D
cfgtemp         = [];
cfgtemp.trials  = itrial;
cfgtemp.channel = 'M1D';
sel             = ft_selectdata(cfgtemp, LFP{1}.Seizure);

cfgtemp         = [];
cfgtemp.lpfilter = 'yes';
cfgtemp.lpfreq  = 50;
sel = ft_preprocessing(cfgtemp, sel);

subplot(4,1,3);
plot(sel.time{1}, sel.trial{1}, 'k', 'LineWidth', 0.3);
xlim(toi);
ylim(y);
xticks([]);
set(gca, 'XColor', 'white', 'TickDir', 'out');

% remove ekg
cfgtemp         = [];
cfgtemp.trials  = itrial;
cfgtemp.channel = {'EMG1', 'EMG2'};
data_ica        = ft_selectdata(cfgtemp, LFP{3}{1}.Seizure);
data_ica        = ft_componentanalysis([], data_ica);

subplot(4,1,4);hold on
plot(data_ica.time{1}, data_ica.trial{1}(2,:), 'LineWidth', 0.3);
xlim(toi);
ylim([-1500 1500]);
xlabel('Time (s)');

fname = fullfile(config{ipatient}.imagesavedir, '..', sprintf('TFR_example_for_figure_%strial%d',config{ipatient}.prefix, itrial));
dtx_savefigure(fig, fname, 'png', 'pdf', 'close');

%% count seizures and emg without artefacts
for ipatient = 1:size(config,2)
    for markername = string(config{ipatient}.LFP.name(1:2))
        if ~isfield(LFP{ipatient}{ipart},markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        switch markername
            case 'SlowWave'
                count.without_artefacts.seizures.all(ipatient)= sum(LFP{ipatient}{ipart}.(markername).trialinfo.starttime - config{ipatient}.injectiontime < hours(10));
            case 'SlowWave_EMG_begin'
                count.without_artefacts.emg.all(ipatient)     = sum(LFP{ipatient}{ipart}.(markername).trialinfo.starttime - config{ipatient}.injectiontime < hours(10));
        end
    end
    count.without_artefacts.emg_eeg_proportion.all(ipatient) = count.without_artefacts.emg.all(ipatient)/count.without_artefacts.seizures.all(ipatient);
end
count.without_artefacts.seizures.total                      = sum(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.min                        = min(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.med                        = median(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.max                        = max(count.without_artefacts.seizures.all);

count.without_artefacts.emg.total                           = sum(count.without_artefacts.emg.all);
count.without_artefacts.emg.min                             = min(count.without_artefacts.emg.all);
count.without_artefacts.emg.med                             = median(count.without_artefacts.emg.all);
count.without_artefacts.emg.max                             = max(count.without_artefacts.emg.all);
count.without_artefacts.emg.min_non_zero                    = min(count.without_artefacts.emg.all(count.without_artefacts.emg.all>0));
count.without_artefacts.emg.med_non_zero                    = median(count.without_artefacts.emg.all(count.without_artefacts.emg.all>0));

count.without_artefacts.emg_eeg_proportion.min              = min(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.med              = median(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.max              = max(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.mean             = mean(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.mean_non_zero    = mean(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));
count.without_artefacts.emg_eeg_proportion.min_non_zero     = min(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));
count.without_artefacts.emg_eeg_proportion.med_non_zero     = median(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));

%% save count variable
save(fullfile(config{ipatient}.datasavedir, 'allpatients_count_datalength.mat'), 'count'); 