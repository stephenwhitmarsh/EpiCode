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

config = dtx_eegvideo_setparams;
ipart = 1;
setfig = @() set(gca,'TickDir','out','FontWeight','bold', 'FontSize', 15);

pat_list = 1:size(config,2);

%% vérifier les marqueurs
% ipatient=9;
% MuseStruct = readMuseMarkers(config{ipatient},true);
% config{ipatient}.prefix(1:end-1)
% check_nr_crises_startend(config{ipatient},MuseStruct,1);

% concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},true);

%% #DATA# read Muse marker data
for ipatient = 1:size(config,2)
    fprintf('\n************** Read marker data for %s **************\n',config{ipatient}.prefix(1:end-1));
    %read Muse marker and correct eeg file if it is micromed segmented data
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient},false);
    %concatenate muse marekrs
    MuseStruct_concat{ipatient} = concatenateMuseMarkers(config{ipatient},MuseStruct{ipatient},false);    %count seizure infos
end

%% compute seizures and emg descriptive stats
for ipatient = 1:size(config,2)
    %analysis only on the last part 
    ipart = size(MuseStruct{ipatient},2);
    
    %FIXME à déplacer dans le setparams
    %seizure stats

    config{ipatient}.seizuretimings.injection_clock= config{ipatient}.injectiontime;
   
    seizure_timing{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
    
    %emg stats

    config{ipatient}.emgtimings.injection_clock= config{ipatient}.injectiontime;
    emg_timing{ipatient} = dtx_stats_emg_timings(config{ipatient},MuseStruct_concat{ipatient},ipart);
    
    %count seizures, slowwaves and emg
    count.seizures.all(ipatient)  = size(seizure_timing{ipatient}.time_start.synctime,2);
    count.slowwaves.all(ipatient) = size(MuseStruct_concat{ipatient}{1}.markers.SlowWave.synctime,2);
    count.emg.all(ipatient)       = size(emg_timing{ipatient}.time_start.synctime,2);
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
count.emg.total        = sum(count.emg.all);
count.emg.min          = min(count.emg.all);
count.emg.med          = median(count.emg.all);
count.emg.max          = max(count.emg.all);
count.emg.min_non_zero = min(count.emg.all(count.emg.all>0));
count.emg.med_non_zero = median(count.emg.all(count.emg.all>0));
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
            data.time{ipatient} = hours(seizure_timing{ipatient}.statsovertime.starttime + seizure_timing{ipatient}.statsovertime.endtime) ./ 2; %middle of the window
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
            end
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


%% plot emg duration and eeg-emg delay over time
for iparam = ["eeg_emg_delay_mean", "emg_duration_mean"]
    clear data
    %Crise_Start over time
    fig = figure;
    sgtitle(iparam,'Interpreter', 'none','Fontsize', 18, 'FontWeight', 'bold');
    subplot(2,1,1);hold on
    data.label = {'dummy'};
    i_emg = 0;
    for ipatient = 1:size(config,2)
        if emg_timing{ipatient}.nr_emg == 0
            continue
        end
        i_emg = i_emg+1;
        data.time{i_emg} = hours(emg_timing{ipatient}.statsovertime.starttime + emg_timing{ipatient}.statsovertime.endtime) ./ 2; %middle of the window
        data.trial{i_emg}= emg_timing{ipatient}.statsovertime.(iparam);
        %     idx = y(x<hours(17)) == 0;
        %     y(idx) = nan;
        leg{i_emg} = plot(data.time{i_emg},data.trial{i_emg}, 'k');
        leg{i_emg}.ZData = ones(size(leg{i_emg}.YData));
        %scatter begining
        idxstart = find(~isnan(data.trial{i_emg}),1,'first');
        s = plot(data.time{i_emg}(idxstart),data.trial{i_emg}(idxstart),'o','MarkerEdgeColor','g','MarkerFaceColor','g');
        s.ZData = ones(size(s.YData)).*2;
        %scatter end
        idxend = find(~isnan(data.trial{i_emg}),1,'last');
        s = plot(data.time{i_emg}(idxend),data.trial{i_emg}(idxend),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
        s.ZData = ones(size(s.YData)).*2;
        name{i_emg} = config{ipatient}.prefix(1:end-1);
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
    
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_over_time',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end


%% plot distrib of eeg emg delay, and emg duration
%voir si différencer emg artefacté ou absence d'emg
%eeg emg delay
for iparam = ["eeg_emg_delay", "emg_duration"]
    fig = figure;hold on
    for ipatient = 1:size(config,2)
        bar(ipatient, nanmean(emg_timing{ipatient}.(iparam)),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
        scatter(rand(size(emg_timing{ipatient}.(iparam)))*0.2+ipatient-0.1, emg_timing{ipatient}.(iparam), '.', 'MarkerEdgeColor', 'k');
        errorbar(ipatient, nanmean(emg_timing{ipatient}.(iparam)),0, nanstd(emg_timing{ipatient}.(iparam)),'k','CapSize',10);
    end
    xlim([0 size(config,2)+1]);
    ylabel(iparam, 'Interpreter','none');
    setfig();
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end

%eeg emg delay summary
for iparam = ["eeg_emg_delay", "emg_duration"]
    for ipatient = 1:size(config,2)
        meanparam(ipatient) = nanmean(emg_timing{ipatient}.(iparam));
        stdparam(ipatient)  = nanstd(emg_timing{ipatient}.(iparam));
    end
    fig = figure;hold on
    x = (1:size(config,2))./size(config,2);
    errorbar(x, meanparam,stdparam,'sk','MarkerFaceColor','k');
    errorbar(0.5,nanmean(meanparam),nanstd(meanparam),'-rx','LineWidth',2,'CapSize',10,'MarkerSize',10);
    xlim([-1 2]);
    xticks([]);
    ylabel(sprintf('mean %s, per patient',iparam),'Interpreter','none');
    ax = axis; ylim([0 ax(4)]);
    setfig();
    fname = fullfile(config{ipatient}.imagesavedir,'..','seizure_description',sprintf('allpatients_%s_summary',iparam));
    dtx_savefigure(fig,fname,'pdf','png','fig','close');
end


% for ipatient = 1:size(config,2)
%     fprintf('patient %d : last seizure %g hours post injection\n',ipatient, hours(seizure_timing{ipatient}.time_start.clock(end) - config{ipatient}.injectiontime));
% end

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
            data.time{ipatient} = hours(morpho{ipatient}.(i_filt).statsovertime.starttime + morpho{ipatient}.(i_filt).statsovertime.endtime) ./ 2; %middle of the window
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

%% EMG : compute envelopes
for markername = "SlowWave_EMG_begin"
    
    for ipatient = pat_list
        
        if ~isfield(LFP{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername).label)
            continue
        end
        
        cfgtemp = [];
        cfgtemp.channel = config{ipatient}.EMG.(markername);
        EMG = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
        
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

%% EMG : plot overdraw and avg of envelopes for each rat (no normalization)
for markername = "SlowWave_EMG_begin"
    for ipatient = pat_list
        if ~isfield(LFP{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername).label)
            continue
        end
        
        fig = figure; hold on
        sgtitle(sprintf('%s : %d trials', config{ipatient}.prefix(1:end-1), size(LFP{ipatient}{ipart}.(markername).trial,2)),'Interpreter','none','FontSize',18,'FontWeight','bold');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % emg raw + avg
%         subplot(3,1,1);hold on;
        %plot each trial with baseline substraction
        for itrial = 1:size(env{ipatient}.trial,2)
            bl_trial   = nanmean(env{ipatient}.trial{itrial}(t>-2&t<-0.5));
            trial = env{ipatient}.trial{itrial} - bl_trial;
            p = plot(t,trial,'k');
            p.Color(4) = 0.2;
            p.ZData = ones(size(p.YData));
        end
        
        %plot avg
        EMG_avg_blcorrected = env_avg{ipatient}.avg-bl_avg(ipatient);
        p = plot(t,EMG_avg_blcorrected, 'b','LineWidth',2);
        p.ZData = ones(size(p.YData)).*2;
        
        %plot std
        std_data = sqrt(env_avg{ipatient}.var);
        x = t;
        y = [EMG_avg_blcorrected - std_data; std_data; std_data]';
        filled_SD = area(x,y);
        filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
        filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
        filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
        filled_SD(1).ShowBaseLine = 'off';
    
        %set figure display
        ax = axis;
        xlim([-1 1]);
        ylim([min(EMG_avg_blcorrected(t>-2&t<2) - std_data(t>-2&t<2)) max(EMG_avg_blcorrected(t>-2&t<2) + std_data(t>-2&t<2))]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca, 'FontWeight','bold', 'Fontsize',15);
        set(gca,'TickDir','out');
%         title('EMG morphology : envelopes of each trial, + average +/- std','Fontsize',18);
%         set(gca,'ycolor','b');
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %eeg raw + avg
%          subplot(3,1,2);hold on;
%          cfgtemp         = [];
%          cfgtemp.channel = config{ipatient}.LFP.motorcortex;
%          LFP_motorcortex = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
%          
%          t = LFP_motorcortex.time{1};
%         %plot each trial with baseline substraction
%         for itrial = 1:size(LFP_motorcortex.trial,2)
%             bl_trial   = nanmean(LFP_motorcortex.trial{itrial}(t>-2&t<-0.5));
%             trial = LFP_motorcortex.trial{itrial} - bl_trial;
%             p = plot(t,trial,'k');
%             p.Color(4) = 0.2;
%             p.ZData = ones(size(p.YData));
%         end
%         
%         %plot avg
%         LFP_avg = ft_timelockanalysis([],LFP_motorcortex);
%         LFP_avg_blcorrected = LFP_avg.avg-nanmean(LFP_avg.avg(t>-2&t<-0.5));
%         p = plot(t,LFP_avg_blcorrected, 'r','LineWidth',2);
%         p.ZData = ones(size(p.YData)).*2;
%         
%         %plot std
%         std_data = sqrt(LFP_avg.var);
%         x = t;
%         y = [LFP_avg_blcorrected - std_data; std_data; std_data]';
%         filled_SD = area(x,y);
%         filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%         filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%         filled_SD(2).FaceColor = 'r'; filled_SD(3).FaceColor = 'r';
%         filled_SD(1).ShowBaseLine = 'off';
%     
%         %set figure display
%         ax = axis;
%         xlim([-2 2])
%         ylim([min(LFP_avg_blcorrected(t>-2&t<2) - std_data(t>-2&t<2)) max(LFP_avg_blcorrected(t>-2&t<2) + std_data(t>-2&t<2))]);
%         xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
%         ylabel('uV');
%         set(gca, 'FontWeight','bold', 'Fontsize',15);
%         set(gca,'TickDir','out');
%         title('EEG morphology (raw trials, + average +/- std)','Fontsize',18);
% %         set(gca,'ycolor','r');
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %eeg/emg comparison
%         subplot(3,1,3);hold on;
%         
%         %plot avg EEG mean + std
% %         yyaxis left
%         plot(t,-LFP_avg_blcorrected./max(-LFP_avg_blcorrected), 'r','LineWidth',2);
%         std_data = sqrt(LFP_avg.var)./max(-LFP_avg_blcorrected);
%         x = t;
%         y = [-LFP_avg_blcorrected./max(-LFP_avg_blcorrected) - std_data; std_data; std_data]';
%         filled_SD = area(x,y);
%         filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%         filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%         filled_SD(2).FaceColor = 'r'; filled_SD(3).FaceColor = 'r';
%         filled_SD(1).ShowBaseLine = 'off';
% %         
% %         axis tight
% %         xlim([-2 2])
% %         ylabel('uV');
% %         set(gca, 'FontWeight','bold', 'Fontsize',15);
% %         set(gca,'TickDir','out');
% %         set(gca,'ycolor','r');
%         
%         % plot EMG mean + std
% %         yyaxis right
%         plot(t,EMG_avg_blcorrected./max(EMG_avg_blcorrected), 'b','LineWidth',2);
%         std_data = sqrt(env_avg{ipatient}.var)./max(EMG_avg_blcorrected);
%         x = t;
%         y = [EMG_avg_blcorrected./max(EMG_avg_blcorrected) - std_data; std_data; std_data]';
%         filled_SD = area(x,y);
%         filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%         filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%         filled_SD(2).FaceColor = 'b'; filled_SD(3).FaceColor = 'b';
%         filled_SD(1).ShowBaseLine = 'off';
%         
%         axis tight
%         xlim([-2 2])
%         ylabel('uV');
%         set(gca, 'FontWeight','bold', 'Fontsize',15);
%         set(gca,'TickDir','out');
% %         set(gca,'ycolor','b');        
%         xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
%         ylim([-1 2]);
        
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho', sprintf('%semg_morphology',config{ipatient}.prefix));
        dtx_savefigure(fig,fname, 'png','close'); %do not save in pdf as it takes 6 hours per figure
    end
end

%% EMG : plot envelope for each dir
for markername = "SlowWave_EMG_begin"
    for ipatient = pat_list
        if ~isfield(LFP{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername).label)
            continue
        end
        
        for idir = unique(LFP{ipatient}{ipart}.(markername).trialinfo.idir)'
            fig = figure; hold on
            sel = find(LFP{ipatient}{ipart}.(markername).trialinfo.idir==idir)';
           
            %plot each trial with baseline substraction
            for itrial = 1:size(sel,2)
                bl_trial   = nanmean(env{ipatient}.trial{sel(itrial)}(t>-2&t<-0.5));
                trial = env{ipatient}.trial{sel(itrial)} - bl_trial;
                p = plot(t,trial,'k');
                p.Color(4) = 0.2;
                p.ZData = ones(size(p.YData));
            end
            
       
            %set figure display
            ax = axis;
            xlim([-2 2]);
            xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
            ylabel('uV');
            set(gca, 'FontWeight','bold', 'Fontsize',15);
            set(gca,'TickDir','out');
            title(sprintf('%s : %s (dir n° %d)', config{ipatient}.prefix(1:end-1),config{ipatient}.directorylist{1}{idir}, idir));
            
            
            fname = fullfile(config{ipatient}.imagesavedir,'emg_morpho_eachdir', sprintf('%semg_morphology_%d',config{ipatient}.prefix,idir));
            dtx_savefigure(fig,fname, 'png','close'); %do not save in pdf as it takes 6 hours per figure
        end%idir
    end
end


%% EMG : plot averages of envelopes of all rats

markername = "SlowWave_EMG_begin";

for plot_std = [true false]
    count_plot = 0;
    fig=figure;hold;
    for ipatient = pat_list
        
        if ~isfield(LFP{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername).label)
            continue
        end
        
        count_plot = count_plot+1;
        
        if plot_std
            h=3;
        else
            h=1;
        end
        
        %compute normalize avg and std
        env_avg_norm = (env_avg{ipatient}.avg - bl_avg(ipatient)) ./ norm_factor(ipatient);
        env_std_norm = sqrt(env_avg{ipatient}.var) ./ norm_factor(ipatient);
        
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
        
        y_labels{count_plot} = sprintf('pat %d (n=%d)',ipatient, size(env{ipatient}.trial,2));
        patient{count_plot} = config{ipatient}.prefix(1:end-1);
        
    end
    
    ax = axis;
    xlim([-1 1])
    %ylim([h-h/4 ax(4)]);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    set(gca,'TickDir','out');
    yticks(h:h:count_plot*h);
    yticklabels(y_labels);
    
    if plot_std
        title('EMG morphology : average+/-std of envelopes of each trial','Fontsize',18);
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho', 'allpatients_emg_morphology_withstd');
    else
        title('EMG morphology : average of envelopes of each trial','Fontsize',18);
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho', 'allpatients_emg_morphology');
    end
    
    dtx_savefigure(fig,fname, 'png','pdf','close');
end

%% plot EMG method for 10 randomly selected EMG
count_plot = 0;
for markername = "SlowWave_EMG_begin"
    for ipatient = pat_list
        if ~isfield(LFP{ipatient}{ipart}, markername)
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername))
            continue
        end
        if isempty(LFP{ipatient}{ipart}.(markername).label)
            continue
        end
        
        cfgtemp = [];
        cfgtemp.channel = config{ipatient}.EMG.(markername);
        EMG = ft_selectdata(cfgtemp,LFP{ipatient}{ipart}.(markername));
        
        if isempty(EMG.label)
            continue
        end
        
        fig = figure;hold on;
        t = EMG.time{1};
        %select 10 emg
        nb_trials = 10;
        trial_list = randperm(size(EMG.trial,2),nb_trials);
        
        %compute all envelopes
        env = [];
        for itrial = 1:nb_trials%trial_list
            rect_emg = abs(EMG.trial{trial_list(itrial)}(1,:));
            [env{itrial}, ~] = envelope(rect_emg,config{ipatient}.EMG.envparam,config{ipatient}.EMG.envmethod);
            %plot(t,env{itrial},'color','k');
        end
        
        %h automatic setting :
        for itrial = 1:nb_trials
            h_temp_max = max(EMG.trial{trial_list(itrial)}(1,:));
            h_temp_min = min(EMG.trial{trial_list(itrial)}(1,:));
            h_temp_amplitude(itrial) = h_temp_max - h_temp_min;
        end
        h = mean(h_temp_amplitude);
        
        % plot raw emg
        subplot(1,3,1); hold on
        for itrial = 1:nb_trials
            plot(t,EMG.trial{trial_list(itrial)}(1,:)+(nb_trials+1)*h- itrial*h,'k'); %first on top
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        %xlabel(sprintf('Time from \n%s (s)', config{ipatient}.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
        ylabel('Number of seizures');
        title(sprintf('Timecourse of \n%s (%d trials)', EMG.label{1}, size(EMG.trial,2)),'Interpreter','none','Fontsize',15);
        setfig();
        tick = h;
        yticks(tick : tick*10 : nb_trials*h);
        yticklabels(nb_trials : -10 : 0);
        axis tight
        xlim([-2 2]);
        
        % Plot rectified EMG
        h = h/2; %because no more negative values
        
        subplot(1,3,2); hold on
        
        for itrial = 1:nb_trials
            plot(t,abs(EMG.trial{trial_list(itrial)}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        xlabel(sprintf('Time from %s (s)', markername),'Interpreter','none');
        title(sprintf('Timecourse of \nrectified %s (%d trials)', EMG.label{1}, size(EMG.trial,2)),'Interpreter','none','Fontsize',15);
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
        
        for itrial = 1:nb_trials
            plot(t,abs(EMG.trial{trial_list(itrial)}(1,:))+ (nb_trials+1)*h - itrial*h,'k'); %first on top
            plot(t,env{itrial}+ (nb_trials+1)*h - itrial*h,'c','LineWidth',2);
        end
        plot([0 0],[h/2 (nb_trials+1)*h], '--r');
        
        title(sprintf('Timecourse of envelope \nof rectified %s (%d trials)', EMG.label{1}, size(EMG.trial,2)),'Interpreter','none','Fontsize',15);
        setfig();
        tick = h;
        yticks(tick : tick*10 : nb_trials*h);
        yticklabels(nb_trials : -10 : 0);
        set(gca,'TickDir','out');
        axis tight
        xlim([-2 2]);
        
        %print figure
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_method',[config{ipatient}.prefix,convertStringsToChars(markername),'_emg_method_',EMG.label{1}]);
        dtx_savefigure(fig,fname,'pdf','png','close');
    end
end

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

if false
    for ipatient = 7
        fname = fullfile(config{ipatient}.datasavedir,sprintf('%sLFP_cleanedforTFR.mat',config{ipatient}.prefix));
        load(fname,'data_cleaned');
        for markername = ["SlowWave_begin", "SlowWave", "Crise_End"]
            waitfor(msgbox(sprintf('%s, p%d, %s : reject trials and/or channels',config{ipatient}.prefix(1:end-1),ipart,markername)));
            cfgtemp             = [];
            cfgtemp.method      = 'trial';
            cfgtemp.channel     = config{ipatient}.align.channel.SlowWave;
            cfgtemp.latency     = 'all';%only for vizualization (not selection of period)
            cfgtemp.box         = 'no';
            cfgtemp.axis        = 'no';%do not work, so I modified it in the defaults of ft_plot_vector
            data_cleaned.(markername) = ft_rejectvisual(cfgtemp,data_cleaned.(markername));                %apply tyhis selection to CSD
        end
        %save cleaned data
        fname = fullfile(config{ipatient}.datasavedir,sprintf('%sLFP_cleanedforTFR.mat',config{ipatient}.prefix));
        save(fname,'data_cleaned', '-v7.3');
        
    end
end

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
                count.without_artefacts.seizures.all(ipatient)    = size(LFP{ipatient}{ipart}.(markername).trial,2);
            case 'SlowWave_EMG_begin'
                count.without_artefacts.emg.all(ipatient)     = size(LFP{ipatient}{ipart}.(markername).trial,2);
        end
    end
    count.without_artefacts.emg_eeg_proportion.all(ipatient) = count.without_artefacts.emg.all(ipatient)/count.without_artefacts.seizures.all(ipatient);
end
count.without_artefacts.seizures.total    = sum(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.min      = min(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.med      = median(count.without_artefacts.seizures.all);
count.without_artefacts.seizures.max      = max(count.without_artefacts.seizures.all);

count.without_artefacts.emg.total    = sum(count.without_artefacts.emg.all);
count.without_artefacts.emg.min      = min(count.without_artefacts.emg.all);
count.without_artefacts.emg.med      = median(count.without_artefacts.emg.all);
count.without_artefacts.emg.max      = max(count.without_artefacts.emg.all);
count.without_artefacts.emg.min_non_zero = min(count.without_artefacts.emg.all(count.without_artefacts.emg.all>0));
count.without_artefacts.emg.med_non_zero = median(count.without_artefacts.emg.all(count.without_artefacts.emg.all>0));

count.without_artefacts.emg_eeg_proportion.min = min(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.med = median(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.max = max(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.mean = mean(count.without_artefacts.emg_eeg_proportion.all);
count.without_artefacts.emg_eeg_proportion.mean_non_zero = mean(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));
count.without_artefacts.emg_eeg_proportion.min_non_zero = min(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));
count.without_artefacts.emg_eeg_proportion.med_non_zero = median(count.without_artefacts.emg_eeg_proportion.all(count.without_artefacts.emg_eeg_proportion.all>0));


%% save count variable
%computed in several separated parts of this script :
% - count.data_length
% - count.seizures
% - count.emg
% - count.without_artefacts.seizures
% - count.without_artefacts.emg
save(fullfile(config{ipatient}.datasavedir, 'allpatients_count_datalength.mat'), 'count'); 









