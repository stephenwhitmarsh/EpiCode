% LGI1 patients data analysis

%% Set parameters
if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\distributionPlot;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/distributionPlot/
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

config = dtx_patients_lgi1_setparams;
ipart=1;

%% #DATA# read Muse marker data

for ipatient = 1:size(config,2)
    fprintf('\n Read marker data for %s\n',config{ipatient}.prefix(1:end-1));
    %read Muse marker and correct eeg file if it is micromed segmented data
    MuseStruct = readMuseMarkers_discontinuousMicromed(config{ipatient},false);
    %concatenate muse markers
    MuseStruct_concat = concatenateMuseMarkers(config{ipatient},MuseStruct,false);    %count seizure infos
   
    % compute seizures descriptive stats
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_L';
    seizure_timing_L{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat,ipart);
    
    config{ipatient}.seizuretimings.marker_start = 'SlowWave_R';
    seizure_timing_R{ipatient} = dtx_stats_seizure_timings(config{ipatient},MuseStruct_concat,ipart);
        
    %symplify marker structures
    marker.synctime = [];
    marker.clock    = datetime.empty;
    marker.dir      = [];
    slowwave_begin_R(ipatient)  = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_R_begin',        marker);  
    emg_begin_R(ipatient)       = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_R_EMG__START__', marker);  
    emg_end_R(ipatient)         = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_R_EMG__END__',   marker);  
    slowwave_begin_L(ipatient)  = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_L_begin',        marker);  
    emg_begin_L(ipatient)       = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_L_EMG__START__', marker);  
    emg_end_L(ipatient)         = ft_getopt(MuseStruct_concat{ipart}.markers, 'SlowWave_L_EMG__END__',   marker); 
    
end

%% eeg emg delay, emg duration, and count emg number
for ipatient = 1:size(config,2)
    
    %put right and left seizures together
    slowwave_begin  = [slowwave_begin_R(ipatient).synctime,slowwave_begin_L(ipatient).synctime];
    emg_begin       = [emg_begin_R(ipatient).synctime,emg_begin_L(ipatient).synctime];
    emg_end         = [emg_end_R(ipatient).synctime,emg_end_L(ipatient).synctime];
    if isempty(slowwave_begin) %some patients do not have EMG
        delays{ipatient} = NaN;
        emg_duration{ipatient} = NaN;
        continue
    end
    
    delays{ipatient}        = emg_begin - slowwave_begin;
    emg_duration{ipatient}  = emg_end - emg_begin;
    eeg_dir{ipatient}       = [slowwave_begin_R(ipatient).dir,slowwave_begin_L(ipatient).dir];
    
    %count emg
    count.n_emg.left(ipatient)  = size(~isnan(emg_begin_L(ipatient).synctime),2);
    count.n_emg.right(ipatient) = size(~isnan(emg_begin_R(ipatient).synctime),2);
    count.n_emg.all(ipatient)   = count.n_emg.left(ipatient) + count.n_emg.right(ipatient);

end
count.n_emg.total = sum(count.n_emg.all);
count.n_emg.med = median(count.n_emg.all);
count.n_emg.med_non_zero = median(count.n_emg.all(count.n_emg.all>0));

%emg duration : boxplot of each patient
fig = figure;hold on;
y = [emg_duration{:}]';
g = [];
for i = 1:length(emg_duration)
    g = [g; repmat(i, length(emg_duration{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
boxplot(y, g);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markersize', 1);
ylim([-1 3]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_emg_duration_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig');

%emg duration : table
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
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_emg_duration_table.csv');
writetable(tableout, fname);

% emg duration distribution of patients' means
fig = figure; hold on;
y = [];
for ipatient = 1:size(config,2)
    y = [y, mean(emg_duration{ipatient}, 'omitnan')];
end
scatter(rand(size(y))*0.1-0.05+1, y, 50, 'ok', 'filled');
boxplot(y);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylim([-1 3]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_emg_duration_boxplot_average');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

% eeg-emg delay : boxplot of each patient
fig = figure;hold on;
y = [delays{:}]';
g = [];
for i = 1:length(delays)
    g = [g; repmat(i, length(delays{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
boxplot(y, g);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markersize', 1);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_eeg_emg_delays_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig');

% eeg-emg delay : table
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
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_eeg-emg-delays_table.csv');
writetable(tableout, fname);

%% seizure frequency and regularity, and count seizure number

%count seizures
for ipatient = 1:size(config,2)
    count.n_seizures.left(ipatient)  = size(~isnan(seizure_timing_L{ipatient}.time_start.synctime),2);
    count.n_seizures.right(ipatient) = size(~isnan(seizure_timing_R{ipatient}.time_start.synctime),2);
    count.n_seizures.all(ipatient)   = count.n_seizures.left(ipatient) + count.n_seizures.right(ipatient);
end
count.n_seizures.total = sum(count.n_seizures.all);
count.n_seizures.med = median(count.n_seizures.all);

%get seizure starts, for right and left TDW
for ipatient = 1:size(config,2)
    % sort data and remove intervals falling between 2 directories.
    % right and left TDW together :
    [data, indexes] = sort([seizure_timing_R{ipatient}.time_start.clock, seizure_timing_L{ipatient}.time_start.clock]);
    dirtemp = [seizure_timing_R{ipatient}.time_start.dir, seizure_timing_L{ipatient}.time_start.dir];
    dir_list = dirtemp(indexes);
    data = minutes(diff(data));
    data_L_R{ipatient} = data(~logical(diff(dir_list)));
    % only left TDW :
    data_L{ipatient} = seizure_timing_L{ipatient}.time_start.clock;
    dir_list = seizure_timing_L{ipatient}.time_start.dir;
    data_L{ipatient} = minutes(diff(data_L{ipatient}));
    data_L{ipatient} = data_L{ipatient}(~logical(diff(dir_list)));
    % only right TDW :
    data_R{ipatient} = seizure_timing_R{ipatient}.time_start.clock;
    dir_list = seizure_timing_R{ipatient}.time_start.dir;
    data_R{ipatient} = minutes(diff(data_R{ipatient}));
    data_R{ipatient} = data_R{ipatient}(~logical(diff(dir_list)));

end

% time between 2 seizures : boxplot of each patient
fig = figure;hold on;
y = [data_L_R{:}]';
g = [];
for i = 1:length(data_L_R)
    g = [g; repmat(i, length(data_L_R{i}), 1)];
end
g(isnan(y)) = [];
y(isnan(y)) = [];
boxplot(y, g, 'symbol', 'ok');
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
ylim([-10 45]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_R_L_sansAOU_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

% time between 2 seizures : table
h = findobj(gca,'Type','line');
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
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_R_L_sansAOU_table.csv');
writetable(tableout, fname, 'delimiter', ';');

% interval between 2 seizures, means of all patients
fig = figure; hold on;
y = [];
for ipatient = 1:size(config,2)
    y = [y, mean(data_L_R{ipatient}, 'omitnan')];
end
scatter(rand(size(y))*0.1-0.05+1, y, 50, 'ok', 'filled');
boxplot(y);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylim([-10 45]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_timebetween_2_SlowWaves_R_L_sansAOU_boxplot_average');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%cv2 of only right TDWs, only left TDWs, and pooled R&L TDWs
for ipatient = 1:size(config,2)
    data_L_R_cv2{ipatient} = cv2(data_L_R{ipatient});
    data_L_cv2{ipatient} = cv2(data_L{ipatient});
    data_R_cv2{ipatient} = cv2(data_R{ipatient});
end

% TDWs CV2 : boxplot from each patient
g = [];
y = [];
for idata = {data_L_cv2, data_R_cv2}
    y = [y, idata{1}{:}];
    for i = 1:length(idata{1})
        g = [g; repmat(i, length(idata{1}{i}), 1)];
    end
end
y = y';
g(isnan(y)) = [];
y(isnan(y)) = [];
fig = figure; hold on;
boxplot(y, g, 'symbol', 'ok');
ylim([0 2]);
x = xlim;
plot(x, [1 1], '--', 'color', [0.6 0.6 0.6]);
set(gca, 'tickdir', 'out', 'fontsize', 25);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_cv2_crises_boxplot');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

% TDWs CV2, mean of all patients
fig = figure; hold on;
avg = [];
for i = unique(g)'
    avg = [avg, mean(y(g==i), 'omitnan')];
end
scatter(rand(size(avg))*0.1-0.05+1, avg, 50, 'ok', 'filled');
boxplot(avg);
set(gca, 'tickdir', 'out', 'fontsize', 25);
x = xlim;
plot(x, [1 1], '--', 'color', [0.6 0.6 0.6]);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylim([0 2]);
fname = fullfile(config{ipatient}.imagesavedir,'..','seizures_description','allpatients_cv2_crises_boxplot_average');
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%% #DATA# read and align LFP, compute surface laplacian, visual rejection of artefacts (load precomputed data)
if false %do not recompute values if it has already been done
    for ipatient = 1:size(config,2)
        
        MuseStruct = readMuseMarkers(config{ipatient},false);
        MuseStruct = alignMuseMarkersPeaks(config{ipatient},MuseStruct,true);
        LFP = readLFP(config{ipatient},MuseStruct,false);
        
        % the following file was donwloaded from : https://github.com/fieldtrip/fieldtrip/blob/master/template/neighbours/elec1020_neighb.mat
        load('elec1020_neighb.mat','neighbours');
        clear data*
        
        for ipart = 1:size(LFP,2)
            for markername = string(config{ipatient}.LFP.name)
                if ~isfield(LFP{ipart},markername) %some patients have only TDW on one side (right or left)
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
                
                %compute surface laplacian
                cfgtemp                     = [];
                cfgtemp.method              = 'spline';
                cfgtemp.elec                = ft_read_sens('standard_1020.elc'); %standard_1020.elc was already in the Fieldtrip's template folder
                cfgtemp.neighbours          = neighbours;
                CSD_temp                    = ft_scalpcurrentdensity(cfgtemp,LFP_temp);
                
                % select only time of interest for vizualization
                cfgtemp         = [];
                cfgtemp.latency = [-1 0.5];
                cfgtemp.channel = {'all', '-EMG*'};
                data_sel_LFP    = ft_selectdata(cfgtemp, LFP_temp);
                
                %reject artefacted trials and/or channels
                waitfor(msgbox(sprintf('%s, p%d, %s : reject trials and/or channels',config{ipatient}.prefix(1:end-1),ipart,markername)));
                cfgtemp             = [];
                cfgtemp.method      = 'trial';
                cfgtemp.latency     = [-1 0.5]; %only for vizualization
                cfgtemp.keepchannel = 'no';
                cfgtemp.box         = 'yes';
                cfgtemp.axis        = 'yes'; %this cfg parameter does not work, so I modified it in the defaults of ft_plot_vector
                data_cleaned        = ft_rejectvisual(cfgtemp,data_sel_LFP); 
                
                %reject same trials/channels in LFP and CSD whole data
                cfgtemp = [];
                cfgtemp.trials = data_cleaned.cfg.trials;
                cfgtemp.channel = data_cleaned.label;
                data.LFP{ipart}.(markername) = ft_selectdata(cfgtemp,LFP_temp);
                data.CSD{ipart}.(markername) = ft_selectdata(cfgtemp,CSD_temp);
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

%% TDW topography 
for markername = string(config{ipatient}.LFP.name(1:2))
    for idata = ["LFP", "CSD"]
        fig = figure;
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
            toi = config{ipatient}.morpho.toi.(markername);
            
            iplot = iplot+1;
            subplot(4,3,iplot);hold on;
            cfgtemp = [];
            cfgtemp.layout        = 'EEG1020';
            cfgtemp.colorbar      = 'no';
            cfgtemp.zlim          = 'maxabs';
            cfgtemp.xlim          = toi;
            cfgtemp.comment       = 'xlim';
            cfgtemp.fontsize      = 15;
            cfgtemp.renderer      = 'painters';
            cfgtemp.colormap      = flip(jet(5000)); %flip so negative peak is in yellow
            cfgtemp.comment = 'no';
            ft_topoplotER(cfgtemp,to_plot);
                
            title(sprintf('pat %d (n=%d)',ipatient, size(data{ipatient}.(idata){ipart}.(markername).trial,2)), 'Interpreter','none');
        end
        
        %print to file
        fname = fullfile(config{ipatient}.imagesavedir, '..',sprintf('allpatients_topoplot_%s_%s',markername,idata));
        dtx_savefigure(fig,fname,'pdf','png','close');
        
    end
end

%save a single topoplot, to have the colorbar for the figure
fig = figure;
to_plot.avg = to_plot.avg ./ max(to_plot.avg);
cfgtemp = [];
cfgtemp.layout        = 'EEG1020';
cfgtemp.colorbar      = 'yes';
cfgtemp.zlim          = 'maxabs';
cfgtemp.xlim          = toi;
cfgtemp.comment       = 'xlim';
cfgtemp.fontsize      = 15;
cfgtemp.renderer      = 'painters';
cfgtemp.colormap      = flip(jet(5000)); %flip so negative peak is in yellow
cfgtemp.comment = 'no';
ft_topoplotER(cfgtemp,to_plot);
caxis([-1 1]);
c = colorbar;
c.TickDirection = 'out';
set(gca, 'TickDir', 'out', 'FontSize', 25);
fname = fullfile(config{ipatient}.imagesavedir, '..', 'allpatients_topoplot_colorbar');
dtx_savefigure(fig,fname,'pdf','png','close');

%% #DATA# compute TDW morphology, for each trial (load precomputed data)
if false %do not recompute values if it has already been done
    for ipatient = 1:size(config,2)
        for idata = ["LFP", "CSD"]
            for markername = string(config{ipatient}.LFP.name)
                if ~isfield(data{ipatient}.(idata){ipart},markername) %some patients have only TDW on one side (right or left)
                    continue
                end
                if isempty(data{ipatient}.(idata){ipart}.(markername))
                    continue
                end
                
                fig = figure; hold on;
                
                %compute data and plot results (to check to quality of the detection), for each trial
                for itrial = 1 : size(data{ipatient}.(idata){ipart}.(markername).trial,2)
                    cfgtemp        = [];
                    cfgtemp.trials = itrial;
                    data_trial     = ft_selectdata(cfgtemp, data{ipatient}.(idata){ipart}.(markername));
                    
                    config{ipatient}.morpho.toiac    = [config{ipatient}.morpho.toi.(markername)(1) - 0.5, config{ipatient}.morpho.toi.(markername)(2)+0.2];
                    config{ipatient}.morpho.channame = config{ipatient}.align.channel.(markername);
                    try
                        [hw, amp] = plot_morpho(config{ipatient}, data_trial);
                    catch
                        warning('cannot find slowwave in trial %d', itrial);
                        hw = nan;
                        amp = nan;
                    end
                    morpho.(markername).(idata).halfwidth{ipatient}(itrial) = hw;
                    morpho.(markername).(idata).amplitude{ipatient}(itrial) = amp;
                end
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
    %load precomputed amplitude and halfwidth
    load(fullfile(config{ipatient}.datasavedir,'allpatients_slowwave_morpho.mat'), 'morpho');
end

%% plot TDW halfwidth

iparam = "halfwidth";
idata = "LFP";
y = [];
g = [];
for markername = ["SlowWave_L", "SlowWave_R"]
    y = [y; (abs([morpho.(markername).(idata).(iparam){:}]))'];
    for i = 1:length(morpho.(markername).(idata).(iparam))
        g = [g; repmat(i, length(morpho.(markername).(idata).(iparam){i}), 1)];
    end
end

% boxplot for each patient
fig = figure;hold on;
g(isnan(y)) = [];
y(isnan(y)) = [];
boxplot(y, g, 'symbol', 'ok');
set(gca, 'tickdir', 'out', 'fontsize', 15);
set(findall(gca, 'type', 'line'), 'linewidth', 2, 'markerfacecolor', 'k', 'markersize', 1);
ylabel('s');
ylim([-0.1 1.5]);
fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_RL_%s_%s_boxplot',idata,iparam));
dtx_savefigure(fig,fname,'pdf','png','fig');

% output table
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
fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_RL_%s_%s_table.csv',idata,iparam));
writetable(tableout, fname, 'delimiter', ';');

% boxplot with the mean of each patient
avg = [];
for i = unique(g)'
    avg = [avg, mean(y(g==i), 'omitnan')];
end
fig = figure; hold on;
scatter(rand(size(avg))*0.1-0.05+1, avg, 50, 'ok', 'filled');
boxplot(avg, 'symbol', 'k');
set(gca, 'tickdir', 'out', 'fontsize', 15);
set(findall(gca, 'type', 'line'), 'linewidth', 2);
ylabel('s');
ylim([-0.1 1.5]);
fname = fullfile(config{ipatient}.imagesavedir,'..',sprintf('allpatients_morpho_RL_%s_%s_boxplot_average',idata,iparam));
dtx_savefigure(fig,fname,'pdf','png','fig','close');

%% plot TDW waveforms

%overdraw of all trials, per patient
for markername = ["SlowWave_R", "SlowWave_L"]
    for ipatient = 1:size(config, 2)
        if ~isfield(data{ipatient}.LFP{ipart},markername) %some patients have only TDW on one side (right or left)
            continue
        end
        if isempty(data{ipatient}.LFP{ipart}.(markername))
            continue
        end
        cfgtemp         = [];
        cfgtemp.channel = config{ipatient}.align.channel.(markername);
        to_plot         = ft_selectdata(cfgtemp, data{ipatient}.LFP{1}.(markername));
        y   = cat(1,to_plot.trial{:});
        x   = to_plot.time{1};
        fig = figure; hold on;
        plot(x,y,'k');
        y   = mean(y, 1, 'omitnan');
        plot(x,y, 'Color', [0.2, 0.6, 1], 'LineWidth', 2);
        axis tight;
        xlim([-2 2]);
        set(gca, 'TickDir', 'out', 'FontSize', 15, 'FontWeight', 'bold');
        xlabel('Time (s)');
        ylabel('uV');
        fname = fullfile(config{ipatient}.imagesavedir, '..', 'figure_morpho', [config{ipatient}.prefix,'SlowWave_morpho_', char(markername)]);
        dtx_savefigure(fig,fname, 'png', 'pdf', 'close');
    end
end

% overdraw of all averages (normalized)
fig = figure; hold on;
removed = 0;
kept = 0;
for ipatient = 1:size(config, 2)
    y = {};
    for markername = ["SlowWave_R", "SlowWave_L"]
        if ~isfield(data{ipatient}.LFP{ipart},markername)
            continue
        end
        if isempty(data{ipatient}.LFP{ipart}.(markername))
            continue
        end
        cfgtemp         = [];
        cfgtemp.channel = config{ipatient}.align.channel.(markername);
        to_plot         = ft_selectdata(cfgtemp, data{ipatient}.LFP{1}.(markername));
        y{end+1} = cat(1,to_plot.trial{:});
    end
    y   = mean(cat(1,y{:}), 1, 'omitnan');
    x   = to_plot.time{1};
    
    [normval, normloc] = min(y(x>-0.4 & x<0.4), [], 'omitnan');
    temp = x(x>-0.4 & x<0.4);
    x = x - temp(normloc);
    y = y ./ -normval;
    plot(x, y, 'LineWidth', 2);
end

axis tight;
xlim([-2 2]);
set(gca, 'TickDir', 'out', 'FontSize', 25, 'FontWeight', 'bold');
xlabel('Time (s)');
ylabel('uV');
fname = fullfile(config{ipatient}.imagesavedir, '..', 'figure_morpho', 'Allpatients-SlowWave_morpho_averages_colored');
dtx_savefigure(fig,fname, 'png', 'pdf', 'close');

%% #DATA# read and align EMG data (load precomputed data)
ipatient = 1;
if false
    %get EMG data
    clear data_EMG
    for ipatient = 1:size(config,2)
        MuseStruct = readMuseMarkers(config{ipatient},false);
        MuseStruct = alignMuseMarkersPeaks(config{ipatient},MuseStruct,false);
        LFP = readLFP(config{ipatient},MuseStruct,false);
        %remove EMG with BAD markers
        LFP = removeArtefactedTrials(config{ipatient},LFP);
        for markername = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align", "SlowWave_R", "SlowWave_L"]
            if isempty(LFP{1}.(markername))
                data_EMG{ipatient}.(markername) = [];
                continue
            end
            data_EMG{ipatient}.(markername) = LFP{1}.(markername);
        end
    end
    save(fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'), 'data_EMG', '-v7.3');
else
    fprintf('Reading %s\n',fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'));
    load(fullfile(config{ipatient}.datasavedir, 'allpatients_emg_data.mat'), 'data_EMG');
end

%% EMG : compute envelopes
clear env*
for markername = ["SlowWave_R_EMG_align", "SlowWave_L_EMG_align", "SlowWave_R", "SlowWave_L"]
    for ipatient = 1:size(config,2)
        if contains(markername, 'align')
            markertemp = markername;
        else
            markertemp = sprintf('%s_EMG_align', markername);
        end
        
        %some patients without EMG, or EMG only on one side : 
        if ipatient > size(data_EMG,2)
            continue
        end
        if ~isfield(data_EMG{ipatient}, markertemp)
            continue
        end
        if isempty(data_EMG{ipatient}.(markertemp))
            continue
        end
        if isempty(data_EMG{ipatient}.(markertemp).label)
            continue
        end
        
%         %find non-artefacted EMG trials
%         if contains(markername, "align")
%             triallist = 1:size(data_EMG{ipatient}.(markername).trial, 2);
%         else
%             triallist = [];
%             MuseStruct = readMuseMarkers(config{ipatient}, false);
%             for itrial = 1:size(data_EMG{ipatient}.(markername).trial, 2)
%                 data_cleaned_start  = data_EMG{ipatient}.(markername).trialinfo.begsample/data_EMG{ipatient}.(markername).fsample;
%                 idir                = data_EMG{ipatient}.(markername).trialinfo.idir(itrial);
%                 data_cleaned_start  = seconds(data_cleaned_start) + MuseStruct{1}{idir}.starttime;
%                 starttime = data_EMG{ipatient}.(markername).trialinfo.starttime(itrial);
%                 difftime = seconds(min(abs(data_cleaned_start - starttime)));
%                 if difftime < 10
%                     triallist = [triallist itrial];
%                 end
%             end
%         end
        
        cfgtemp         = [];
        cfgtemp.channel = config{ipatient}.EMG.(markername);
        EMG             = ft_selectdata(cfgtemp,data_EMG{ipatient}.(markername));
        
        if isempty(EMG.label)
            continue
        end
        
        t = EMG.time{1};
        
        %compute all envelopes, and average
        i=0;
        for itrial = 1 : size(EMG.trial,2)
            i = i+1;
            rect_emg                    = abs(EMG.trial{itrial}(1,:));
            [env{ipatient}.(markername).trial{i}, ~] = envelope(rect_emg, config{ipatient}.EMG.envparam, config{ipatient}.EMG.envmethod);
            env{ipatient}.(markername).time{i}       = t;
        end
        env{ipatient}.(markername).label     = {'dummy'};
        env{ipatient}.(markername).trialinfo = EMG.trialinfo;
        
        env_avg{ipatient}.(markername) = ft_timelockanalysis([], env{ipatient}.(markername));
        
        %value to normalize later
        bl_avg.(markername)(ipatient)        = nanmean(env_avg{ipatient}.(markername).avg(t>-2&t<-0.5));
        norm_factor.(markername)(ipatient)   = max(env_avg{ipatient}.(markername).avg - bl_avg.(markername)(ipatient));

    end
end

%% EMG : plot overdraw of each envelope, and avg, for each patient 

for ipatient = 1:size(config,2)
    for markername = ["SlowWave_R_EMG_align", "SlowWave_R", "SlowWave_L_EMG_align", "SlowWave_L"]
        
        if ipatient > size(env,2)
            continue
        end
        if ~isfield(env{ipatient}, markername)
            continue
        end
        if isempty(env{ipatient}.(markername))
            continue
        end
        
        fig = figure; hold on
        sgtitle(sprintf('%s : %d trials', config{ipatient}.prefix(1:end-1), size(env{ipatient}.(markername).trial,2)),'Interpreter','none','FontSize',18,'FontWeight','bold');
        
        %plot each trial
        for itrial = 1:size(env{ipatient}.(markername).trial,2)
            bl_trial   = mean(env{ipatient}.(markername).trial{itrial}(t>-2&t<-0.5), 'omitnan');
            trial = env{ipatient}.(markername).trial{itrial} - bl_trial;
            plot(t,trial,'color', [0.5 0.25 0]);
        end
        
        %plot avg
        EMG_avg_blcorrected = env_avg{ipatient}.(markername).avg-bl_avg.(markername)(ipatient);
        p = plot(t,EMG_avg_blcorrected, 'color', [0.2 0.6 1], 'LineWidth', 2);
        
        %set figure display
        ax = axis;
        xlim([-2 2]);
        xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
        ylabel('uV');
        set(gca, 'FontWeight','bold', 'Fontsize',15);
        set(gca,'TickDir','out');
        
        fname = fullfile(config{ipatient}.imagesavedir,'..','emg_morpho', sprintf('%semg_morphology_figure_%s',config{ipatient}.prefix, markername));
        dtx_savefigure(fig,fname, 'png','close');
    end
end

%% EMG : plot mean envelope of each patient
for markerlist = {["SlowWave_R_EMG_align", "SlowWave_L_EMG_align"], ["SlowWave_R", "SlowWave_L"]}
    fig = figure; hold on
    subplot(2,1,1); hold on;
    clear eegplot
    for ipatient = 1:size(config,2)
        for iside = 1:2 %right and left
            %plot EEG
            if ~isfield(data_EMG{ipatient}, markerlist{1}(iside))
                eegplot{iside} = [];
                continue
            end
            if isempty(data_EMG{ipatient}.(markerlist{1}(iside)))
                eegplot{iside} = [];
                continue
            end
            switch markerlist{1}(iside)
                case "SlowWave_R_EMG_align"
                    electrodetoplot = config{ipatient}.align.channel.SlowWave_R;
                    trialinfo_cleaned = data{ipatient}.LFP{1}.SlowWave_R.trialinfo;
                case "SlowWave_L_EMG_align"
                    electrodetoplot = config{ipatient}.align.channel.SlowWave_L;
                    trialinfo_cleaned = data{ipatient}.LFP{1}.SlowWave_L.trialinfo;
                otherwise
                    electrodetoplot = config{ipatient}.align.channel.(markerlist{1}(iside));
                    trialinfo_cleaned = data_EMG{ipatient}.(markerlist{1}(iside)).trialinfo;
            end
            cfgtemp            = [];
            cfgtemp.channel    = electrodetoplot;
            cfgtemp.keeptrials = 'yes';
            eegplot{iside}     = ft_timelockanalysis(cfgtemp, data_EMG{ipatient}.(markerlist{1}(iside)));
        end
        eegplot = eegplot(~cellfun( @isempty, eegplot));
        if isempty(eegplot)
            continue
        end
        if length(eegplot) > 1
            x = cat(1, eegplot{1}.time, eegplot{2}.time);
            y = cat(1, permute(eegplot{1}.trial, [1 3 2]), permute(eegplot{2}.trial, [1 3 2]));
        else
            x = eegplot{1}.time;
            y = squeeze(eegplot{1}.trial);
        end
        x = x(1,:);
        y = mean(y, 1);
        [normval, normloc] = min(y(x>-0.4 & x<0.4), [], 'omitnan');
        if ~contains(markerlist{1}(iside), 'EMG')
            temp = x(x>-0.4 & x<0.4);
            x = x - temp(normloc);
        end
        y = y ./ -normval; 
        if contains(markerlist{1}(iside), 'EMG') && ipatient == 1
            y = y - mean(y(x>-2 & x<-1.5));
        end
        plot(x, y, 'k');
    end %ipatient eeg
    
    %set figure display
    ax = axis;
    xlim([-2 2]);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    ylabel('uV');
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    set(gca,'TickDir','out');
    
    % plot EMG
    for markername = markerlist{1}
        for ipatient = 1:size(config,2)
            if ipatient > size(env,2)
                continue
            end
            if ~isfield(env{ipatient}, markername)
                continue
            end
            if isempty(env{ipatient}.(markername))
                continue
            end
            %plot normalized EMG
            subplot(2,1,2); hold on;
            env_avg_norm = (env_avg{ipatient}.(markername).avg - bl_avg.(markername)(ipatient)) ./ norm_factor.(markername)(ipatient);
            env_avg_norm = movmean(env_avg_norm, 10);
            p = plot(t,env_avg_norm, 'k');%'color', [0.2 0.6 1]);
        end
    end
    %set figure display
    ax = axis;
    xlim([-2 2]);
    xlabel('Time (s)','Interpreter','none', 'Fontsize',15);
    ylabel('uV');
    set(gca, 'FontWeight','bold', 'Fontsize',15);
    set(gca,'TickDir','out');
    
    label = erase(markername, '_L');
    fname = fullfile(config{ipatient}.imagesavedir,'..','emg_begin_alignment', sprintf('allpatients-semg_morphology_figure_%s', markername));
    dtx_savefigure(fig,fname, 'png', 'pdf', 'close'); 
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

%% count eeg duration per patient and save count data
for ipatient = 1:size(config,2)
    ipart = 1;
    for idir = 1:size(config{ipatient}.directorylist{ipart},2)
        [isNeuralynx, isMicromed, isBrainvision] = get_data_format(config{ipatient});
        % find data file
        if isNeuralynx
            temp     = dir(fullfile(config{ipatient}.rawdir, config{ipatient}.directorylist{ipart}{idir},'*.ncs'));
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
save(fullfile(config{ipatient}.datasavedir, 'allpatients_data_length.mat'), 'count'); 