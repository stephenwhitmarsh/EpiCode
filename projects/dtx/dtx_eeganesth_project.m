% analysis of DTX EEG video rodents

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

config = dtx_eeganesth_setparams;
ipart = 1;
setfig = @() set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 25);

pat_list = 1:size(config,2);

%% vérifier les marqueurs
% for ipatient=1:size(config,2)
%     MuseStruct{ipatient} = readMuseMarkers(config{ipatient},true);
%     config{ipatient}.prefix(1:end-1)
%     check_nr_crises_startend(config{ipatient},MuseStruct{ipatient},1);
%     concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},true);
% end

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\n************** Read marker data for %s **************\n',config{ipatient}.prefix(1:end-1));
    %read Muse marker and correct eeg file if it is micromed segmented data
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient},false);
%     if strcmp(config{ipatient}.type, 'dtx')
%         MuseStruct = dtx_remove_wrong_seizure(config{ipatient}, MuseStruct,true);
%     end
    %concatenate muse marekrs
    MuseStruct_concat{ipatient} = concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},false);    %count seizure infos
end

%% compute seizures descriptive stats
for ipatient = 1:size(config,2)
    %analysis only on the last part 
    ipart = size(MuseStruct{ipatient},2);

    config{ipatient}.seizuretimings.injection_clock= config{ipatient}.injectiontime;
    seizure_timing{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
    
    %count seizures, slowwaves and emg
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

%% plot seizures over time
%seizures occurence
%seizures duration
%seizures cv2 of occurences and duration

for iparam = ["nb_seizures", "timebetween2seizures_mean", "timebetween2seizures_cv2", "seizureduration_mean", "seizureduration_cv2"]
    for do_focus = ["all", "10h"]
    clear data
    %Crise_Start over time
    fig = figure;
    sgtitle(iparam,'Interpreter', 'none','Fontsize', 18, 'FontWeight', 'bold');
    subplot(2,1,1);hold on
    data.label = {'dummy'};
    for ipatient = 1:size(config,2)
        data.time{ipatient} = hours(seizure_timing{ipatient}.statsovertime.endtime); %start of the window
        data.trial{ipatient}= seizure_timing{ipatient}.statsovertime.(iparam);
        %     idx = y(x<hours(17)) == 0;
        %     y(idx) = nan;
        leg{ipatient} = plot(data.time{ipatient},data.trial{ipatient}, 'k');
        leg{ipatient}.ZData = ones(size(leg{ipatient}.YData));
        if strcmp(iparam,'nb_seizures')
            %scatter begining
            idxstart = find(~isnan(data.trial{ipatient}),1,'first');
            s = plot(data.time{ipatient}(idxstart),data.trial{ipatient}(idxstart),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
            s.ZData = ones(size(s.YData)).*2;
            %scatter end
            idxend = find(~isnan(data.trial{ipatient}),1,'last');
            s = plot(data.time{ipatient}(idxend),data.trial{ipatient}(idxend),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
            s.ZData = ones(size(s.YData)).*2;
            name{ipatient} = config{ipatient}.prefix(1:end-1);
        end
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
    %xlim([0 ax(2)]);
    ylabel('raw + mean + std');
    xlabel('time (hours post injection)');
    setfig();
    if strcmp(do_focus,"10h")
        xlim([0 10]);
    end
    
    %plot smoothed avg and std
    subplot(2,1,2);hold on
    avg_data_smooth = movmean(data_avg.avg,10,'omitnan');
    std_data_smooth = movmean(std_data,10,'omitnan');
    
    x = data_avg.time;
    y = [avg_data_smooth - std_data_smooth; std_data_smooth; std_data_smooth]';
    filled_SD = area(x,y);
    filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
    filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
    filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
    filled_SD(1).ShowBaseLine = 'off';
    
    p = plot(data_avg.time, avg_data_smooth, 'Color','r', 'LineWidth', 2);
    p.ZData = ones(size(p.YData));
        
    %set figure display
    axis tight
    ax = axis;
    xlim([0 ax(2)]);
    ylabel('smoothed mean +/- std');
    xlabel('time (hours post injection)');
    setfig();
    if strcmp(do_focus,"10h")
        xlim([0 10]);
    end
    
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_over_time_%s',iparam,do_focus));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
    end
end

%% plot distrib seizures stats for each patient
for iparam = ["timebetween2seizures_mean", "timebetween2seizures_cv2", "seizureduration_mean", "seizureduration_cv2"]
    fig=figure;hold
    for ipatient = 1:size(config,2)
        to_plot = seizure_timing{ipatient}.statsovertime.(iparam);
        bar(ipatient, nanmean(to_plot),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
        scatter(rand(1,size(to_plot,2))*0.2+ipatient-0.1,to_plot, '.', 'MarkerEdgeColor', 'k');
        errorbar(ipatient, nanmean(to_plot), 0,nanstd(to_plot),'k','CapSize',10);
    end
    % xlim([0 size(config,2)+1]);
    title(iparam,'Interpreter','none');
    ax = axis;
    ylim([0 ax(4)]);
    xticks(1:size(config,2));
    xlabel('rat nr');
    setfig();
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_distrib',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end


%% #DATA COMPUTED IN dtx_cluster_eegvideo.m# read and align LFP an remove artefacts (load precomputed data)
% load precomputed data
for ipatient = 1:size(config,2)
    try
    fprintf('Reading %s\n', fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']));
    temp = load(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix,'LFP.mat']), 'LFP');
    LFP{ipatient} = temp.LFP;
    clear temp
    catch
        LFP{ipatient}{1}.SlowWave = [];
        LFP{ipatient}{1}.SlowWave_EMG_begin = [];
    end
end

%% plot seizures for each patient (to select toi morpho and topoplot)
%on plot with overdraw trial by trial, and avg of each channel, h = 120.
%one plot with all avg superposï¿½s

%done with dtx_cluster_eegvideo_plotalldata

%% #DATA COMPUTED IN dtx_cluster_eegvideo.m# slow wave morphology + average over a sliding time window
clear morpho
%only computed for markername "SlowWave"
for ipatient = pat_list
    try
        %load precomputed data
        temp = load(fullfile(config{ipatient}.datasavedir,sprintf('%sslowwave_morpho.mat',config{ipatient}.prefix)), 'morpho');
        morpho{ipatient} = temp.morpho;
    catch
        morpho{ipatient} = [];
    end
end

%% plot slow wave morphology distribution
clear param meanparam stdparam
for iparam = ["halfwidth", "amplitude"]
    for i_filt = ["raw","hpfilt_0_15"]
        fig=figure;hold on
        for ipatient = pat_list
            if isempty(morpho{ipatient})
                continue
            end
            param = abs(morpho{ipatient}.(i_filt).(iparam));
            meanparam.(i_filt).(iparam)(ipatient) = nanmean(param);
            stdparam.(i_filt).(iparam)(ipatient)  = nanstd(param);
            
            bar(ipatient, meanparam.(i_filt).(iparam)(ipatient),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
            scatter(rand(size(param))*0.2+ipatient-0.1, param, '.', 'MarkerEdgeColor', 'k');
            errorbar(ipatient, meanparam.(i_filt).(iparam)(ipatient),0, stdparam.(i_filt).(iparam)(ipatient),'k');
            
        end
        
        ax = axis;
        ylim([0 ax(4)]);
        xlim([0, size(pat_list,2)+1]);
        xticks(1:size(pat_list,2));
        xticklabels(1:size(pat_list,2));
        setfig();
        ylabel(iparam);
        set(gca, 'Fontsize', 11);
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_morpho_%s_%s',iparam,i_filt));
        dtx_savefigure(fig,fname,'pdf','png','close');
    end
end

%plot mean of all patients in the same plot
for iparam = ["halfwidth", "amplitude"]
    for i_filt = ["raw","hpfilt_0_15"]
        fig=figure;hold;
        nb_points = size(meanparam.(i_filt).(iparam),2);
        x = (1:nb_points)./ nb_points .* 0.7 + 0.6;
        errorbar(x, meanparam.(i_filt).(iparam),stdparam.(i_filt).(iparam),'sk','MarkerFaceColor','k');
        errorbar(mean(x),nanmean(meanparam.(i_filt).(iparam)),nanstd(meanparam.(i_filt).(iparam)),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
        ax = axis;
        ylim([0 ax(4)]);
        xlim([0 2]);
        xticks([]);
        setfig();
        ylabel(iparam);
        set(gca,'TickDir','out','FontWeight','bold');
        ax = axis; ylim([0 ax(4)]);
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_morpho_%s_%s_summary',iparam,i_filt));
        dtx_savefigure(fig,fname,'pdf','png','close');
    end
end

%% plot slowwave morphology over time

for i_filt = ["raw","hpfilt_0_15"]
    for iparam = ["halfwidth_mean", "amplitude_mean"]
        clear data
        %Crise_Start over time
        fig = figure;
        sgtitle(sprintf('%s (%s)',iparam,i_filt),'Interpreter', 'none','Fontsize', 18, 'FontWeight', 'bold');
        subplot(2,1,1);hold on
        data.label = {'dummy'};
        for ipatient = 1:size(config,2)
            if isempty(morpho{ipatient})
                data.time{ipatient} = data.time{ipatient-1};
                data.trial{ipatient} = nan(size(data.time{ipatient}));
                continue
            end
            data.time{ipatient} = hours(morpho{ipatient}.(i_filt).statsovertime.endtime); %end of the window
            data.trial{ipatient}= morpho{ipatient}.(i_filt).statsovertime.(iparam);
            %     idx = y(x<hours(17)) == 0;
            %     y(idx) = nan;
            leg{ipatient} = plot(data.time{ipatient},data.trial{ipatient}, 'k');
            leg{ipatient}.ZData = ones(size(leg{ipatient}.YData));
            %scatter begining
            idxstart = find(~isnan(data.trial{ipatient}),1,'first');
            s = plot(data.time{ipatient}(idxstart),data.trial{ipatient}(idxstart),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
            s.ZData = ones(size(s.YData)).*2;
            %scatter end
            idxend = find(~isnan(data.trial{ipatient}),1,'last');
            s = plot(data.time{ipatient}(idxend),data.trial{ipatient}(idxend),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
            s.ZData = ones(size(s.YData)).*2;
            name{ipatient} = config{ipatient}.prefix(1:end-1);
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
        ylabel('raw + mean + std');
        xlabel('time (hours post injection)');
        setfig();
        
        %plot smoothed avg and std
        subplot(2,1,2);hold on
        avg_data_smooth = movmean(data_avg.avg,10,'omitnan');
        std_data_smooth = movmean(std_data,10,'omitnan');
        
        x = data_avg.time;
        y = [avg_data_smooth - std_data_smooth; std_data_smooth; std_data_smooth]';
        filled_SD = area(x,y);
        filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
        filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
        filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
        filled_SD(1).ShowBaseLine = 'off';
        
        p = plot(data_avg.time, avg_data_smooth, 'Color','r', 'LineWidth', 2);
        p.ZData = ones(size(p.YData));
        
        %set figure display
        axis tight
        ax = axis;
        xlim([0 ax(2)]);
        ylabel('smoothed mean +/- std');
        xlabel('time (hours post injection)');
        setfig();
        
        fname = fullfile(config{ipatient}.imagesavedir,'..','morpho_all',sprintf('allpatients_%s_over_time',iparam));
        dtx_savefigure(fig,fname,'pdf','png','fig','close');
    end
end

%% count seizures without artefacts

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
                count.without_artefacts.seizures.all(ipatient)    = size(LFP{ipatient}{ipart}.(markername).trial,2);
        end
    end
end
count.without_artefacts.seizures.total    = sum(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.min      = min(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.med      = median(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.max      = max(count.without_artefacts.seizures.all);


%% save count variable
%computed in several separated parts of this script :
% - count.data_length
% - count.seizures
% - count.emg
% - count.without_artefacts.seizures
% - count.without_artefacts.emg
save(fullfile(config{ipatient}.datasavedir, 'allpatients_count_datalength.mat'), 'count'); 

%% TFR : clean data
if false
    for ipatient = 1:size(config,2)
        i=0;
        for markername = ["SlowWave_begin", "SlowWave", "Crise_End"]
            i=i+1;
            if ~isempty(LFP{ipatient}{1}.(markername))
                test(ipatient,i) = size(LFP{ipatient}{1}.(markername).trial,2);
            else
                test(ipatient,i) = nan;
            end
            a = find(diff(LFP{ipatient}{1}.SlowWave.trialinfo.trialnr)>1);
            b = find(diff(LFP{ipatient}{1}.SlowWave_begin.trialinfo.trialnr)>1);
            
%             cfgtemp = [];
%             cfgtemp.channel     = config{ipatient}.align.channel.SlowWave;
%             cfgtemp.hpfilter    ='yes';
%             cfgtemp.hpfilttype  ='fir';
%             cfgtemp.hpfreq      = config{ipatient}.TFR.foi(1);
%             cfgtemp.lpfilter    ='yes';
%             cfgtemp.lpfilttype  ='fir';
%             cfgtemp.lpfreq      = config{ipatient}.TFR.foi(end);
%             data_sel            = ft_preprocessing(cfgtemp,LFP{ipatient}{1}.(markername));
            
            %reject artefacted trials and/or channels
            waitfor(msgbox(sprintf('%s, p%d, %s : reject trials and/or channels',config{ipatient}.prefix(1:end-1),ipart,markername)));
            cfgtemp             = [];
            cfgtemp.method      = 'summary';%trial
            cfgtemp.channel     = config{ipatient}.align.channel.SlowWave;
            cfgtemp.latency     = 'all';%only for vizualization (not selection of period)
            cfgtemp.box         = 'no';
            cfgtemp.axis        = 'no';%do not work, so I modified it in the defaults of ft_plot_vector
            data_cleaned.(markername) = ft_rejectvisual(cfgtemp,LFP{ipatient}{1}.(markername));                %apply tyhis selection to CSD
        end
        %save cleaned data
        fname = fullfile(config{ipatient}.datasavedir,sprintf('%sLFP_cleanedforTFR.mat',config{ipatient}.prefix));
        save(fname,'data_cleaned', '-v7.3');
        
    end
end








