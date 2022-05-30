function [config] = Figure_LFP_stages

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epishare-master'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/SPIKY_apr_2021'))
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\git\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end

ft_defaults

%% load data
config = hspike_setparams;
for ipatient = 1:8
    MuseStruct{ipatient}        = readMuseMarkers(config{ipatient}, false);
    config{ipatient}.LFP.name   = {'template1','template2','template3','template4','template5','template6'};
%     LFP{ipatient}             = readLFP(config{ipatient}, MuseStruct{ipatient}(1:3), true);
    LFP{ipatient}               = rerefLFP(config{ipatient}, MuseStruct{ipatient}(1:3), false);    
    [clusterindx{ipatient}, trialindx, LFP_cluster{ipatient}] = clusterLFP(config{ipatient});
    LFP_cluster{ipatient}       = LFP_cluster{ipatient}{1}.Hspike.kmedoids{6};    
end

%% rereference clusters

for ipatient = 1 : 8
    % rereference to bipolar
    if strcmp(config{ipatient}.cluster.reref, 'yes')
        for itemplate = 1 : 6
            if contains(LFP_cluster{ipatient}{itemplate}.label, '-')
                disp('Rereferencing already done');
                continue
            end
            disp('Rereferencing');
            labels_nonum    = regexprep(LFP_cluster{ipatient}{itemplate}.label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.channel     = LFP_cluster{ipatient}{itemplate}.label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
            end
            LFP_cluster{ipatient}{itemplate} = ft_appenddata([],group{:});
        end
    end
end

%% demean/baseline correct clusters
for ipatient = 1 : 8
    for itemplate = 1 : 6
        cfgtemp                 = [];
        cfgtemp.demean          = 'yes';
        cfgtemp.baselinewindow  = [-0.15, -0.05];
        LFP_cluster{ipatient}{itemplate} = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
    end
end

%% demean/baseline correct LFP
for ipatient = 1 : 8
    for ipart = 1 : size(LFP{ipatient}, 2)
        for template = ["template1","template2","template3","template4","template5","template6"]
            if isempty(LFP{ipatient}{ipart}.(template))
                continue
            end
            cfgtemp                 = [];
            cfgtemp.demean          = 'yes';
            cfgtemp.baselinewindow  = [-0.15, -0.05];
            LFP{ipatient}{ipart}.(template) = ft_preprocessing(cfgtemp, LFP{ipatient}{ipart}.(template));
        end
    end
end

%% select latency LFPs and clusters
latency = [-0.2 0.5];
cfgtemp = [];
cfgtemp.latency = latency;
for ipatient = 1 : 8
    for ipart = 1 : size(LFP{ipatient}, 2)
        for template = ["template1","template2","template3","template4","template5","template6"]
            if isempty(LFP{ipatient}{ipart}.(template))
                 continue
            end
            LFP{ipatient}{ipart}.(template) = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(template));
        end
    end
    for itemplate = 1:6
        LFP_cluster{ipatient}{itemplate} = ft_selectdata(cfgtemp, LFP_cluster{ipatient}{itemplate});
    end
end

%% concatinate LFP: parts and over templates
for ipatient = 1 : 8
    temp = {};
    for template = ["template1","template2","template3","template4","template5","template6"]
        for ipart = 1 : 3
            if ~isempty(LFP{ipatient}{ipart}.(template))
                LFP{ipatient}{ipart}.(template).trialinfo.part = ones(height(LFP{ipatient}{ipart}.(template).trialinfo), 1) * ipart;
                temp{size(temp, 2)+1} = LFP{ipatient}{ipart}.(template);
            end
        end
    end
    if isempty(temp)
        continue
    end
    LFPstage_overtemplate{ipatient} = ft_appenddata([], temp{:});
end

%% average & sem LFP over templates 
clear temp
for ipatient = 1 : 8
    for stage = unique(LFPstage_overtemplate{ipatient}.trialinfo.hyplabel)'
        cfgtemp = [];
        cfgtemp.trials = strcmp(LFPstage_overtemplate{ipatient}.trialinfo.hyplabel, stage);
        cfgtemp.avgoverrpt = 'yes';
        cfgtemp.channel = 'all';
        LFPstage_overtemplate_avg{ipatient}.(stage) = ft_selectdata(cfgtemp, LFPstage_overtemplate{ipatient});
        
        for ichan = 1:size(LFPstage_overtemplate_avg{ipatient}.(stage).trial{1}, 1)
            cfgtemp = [];
            cfgtemp.trials = strcmp(LFPstage_overtemplate{ipatient}.trialinfo.hyplabel, stage);          
            cfgtemp.channel = ichan;
            cfgtemp.avgoverrpt = 'no';
            temp = ft_selectdata(cfgtemp, LFPstage_overtemplate{ipatient});
            LFPstage_overtemplate_avg{ipatient}.(stage).sem(ichan, :) = std(vertcat(temp.trial{:})) ./ sqrt(size(temp.trial, 2));
        end
    end
end

%% rename and clean up LFP sleepstage names over templates
for ipatient = 1 : 8
    if isempty(LFPstage_overtemplate_avg{ipatient})
        continue
    end
    if isfield(LFPstage_overtemplate_avg{ipatient}, "PHASE_1")
        LFPstage_overtemplate_avg{ipatient}.S1 = LFPstage_overtemplate_avg{ipatient}.PHASE_1;
        LFPstage_overtemplate_avg{ipatient}    = rmfield(LFPstage_overtemplate_avg{ipatient}, 'PHASE_1');
    end
    
    if isfield(LFPstage_overtemplate_avg{ipatient}, "PHASE_2")
        LFPstage_overtemplate_avg{ipatient}.S2      = LFPstage_overtemplate_avg{ipatient}.PHASE_2;
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'PHASE_2');
    end
    if isfield(LFPstage_overtemplate_avg{ipatient}, "PHASE_3")
        LFPstage_overtemplate_avg{ipatient}.S3      = LFPstage_overtemplate_avg{ipatient}.PHASE_3;
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'PHASE_3');
    end
    
    if isfield(LFPstage_overtemplate_avg{ipatient}, "AWAKE")
        LFPstage_overtemplate_avg{ipatient}.Wake    = LFPstage_overtemplate_avg{ipatient}.AWAKE;
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'AWAKE');
    end
    
    if isfield(LFPstage_overtemplate_avg{ipatient}, "POSTSLEEP")
        LFPstage_overtemplate_avg{ipatient}.Post    = LFPstage_overtemplate_avg{ipatient}.POSTSLEEP;
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'POSTSLEEP');
    end
    
    if isfield(LFPstage_overtemplate_avg{ipatient}, "PRESLEEP")
        LFPstage_overtemplate_avg{ipatient}.Pre    = LFPstage_overtemplate_avg{ipatient}.PRESLEEP;
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'PRESLEEP');
    end
    
    if isfield(LFPstage_overtemplate_avg{ipatient}, "NO_SCORE")
        LFPstage_overtemplate_avg{ipatient} = rmfield(LFPstage_overtemplate_avg{ipatient}, 'NO_SCORE');
    end
end

% flip polarity to align spike and wave
for ipatient = 6
    for hyplabel = ["S3", "S2", "S1", "REM", "Pre", "Post", "Wake"]
        if isfield(LFPstage_overtemplate_avg{ipatient}, hyplabel)
            LFPstage_overtemplate_avg{ipatient}.(hyplabel).trial{1} = -LFPstage_overtemplate_avg{ipatient}.(hyplabel).trial{1};
        end
    end
end

%% concatinate LFP: parts per stage
for ipatient = 1 : 8
    for template = ["template1","template2","template3","template4","template5","template6"]
        temp = {};
        for ipart = 1 : 3
            if ~isempty(LFP{ipatient}{ipart}.(template))
                LFP{ipatient}{ipart}.(template).trialinfo.part = ones(height(LFP{ipatient}{ipart}.(template).trialinfo), 1) * ipart;
                temp{size(temp, 2)+1} = LFP{ipatient}{ipart}.(template);
            end
        end
        if isempty(temp)
            continue
        end
        LFPstage{ipatient}.(template) = ft_appenddata([], temp{:});
    end
end

%% average & sem LFP per stage

for ipatient = 1 : 8
    for template = ["template1","template2","template3","template4","template5","template6"]
        cfgtemp = [];
        
        if ~isfield(LFPstage{ipatient}, template)
            continue
        end
        if isempty(LFPstage{ipatient}.(template))
            continue
        end
        for stage = unique(LFPstage{ipatient}.(template).trialinfo.hyplabel)'
            cfgtemp.trials = strcmp(LFPstage{ipatient}.(template).trialinfo.hyplabel, stage);
            cfgtemp.avgoverrpt = 'yes';
            cfgtemp.channel = 'all';
            LFPstageavg{ipatient}.(template).(stage) = ft_selectdata(cfgtemp, LFPstage{ipatient}.(template));
            
            for ichan = 1:size(LFPstageavg{ipatient}.(template).(stage).trial{1}, 1)
                cfgtemp.channel = ichan;            
                cfgtemp.avgoverrpt = 'no';            
                temp = ft_selectdata(cfgtemp, LFPstage{ipatient}.(template)); 
                LFPstageavg{ipatient}.(template).(stage).sem(ichan, :) = std(vertcat(temp.trial{:})) ./ sqrt(size(temp.trial, 2));
            end
        end
    end
end

fname = fullfile(config{1}.datasavedir, 'LFPstageavg.mat');
save(fname, 'LFPstageavg', '-v7.3');

%% rename and clean up sleepstages
for ipatient = 1 : 8
    for template = string(fields(LFPstageavg{ipatient}))'
        if isempty(LFPstageavg{ipatient}.(template))
            continue
        end
        if isfield(LFPstageavg{ipatient}.(template), "PHASE_1")
            LFPstageavg{ipatient}.(template).S1 = LFPstageavg{ipatient}.(template).PHASE_1;
            LFPstageavg{ipatient}.(template)    = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_1');
        end
        
        if isfield(LFPstageavg{ipatient}.(template), "PHASE_2")
            LFPstageavg{ipatient}.(template).S2      = LFPstageavg{ipatient}.(template).PHASE_2;
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_2');
        end
        if isfield(LFPstageavg{ipatient}.(template), "PHASE_3")
            LFPstageavg{ipatient}.(template).S3      = LFPstageavg{ipatient}.(template).PHASE_3;
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PHASE_3');
        end
        
        if isfield(LFPstageavg{ipatient}.(template), "AWAKE")
            LFPstageavg{ipatient}.(template).Wake    = LFPstageavg{ipatient}.(template).AWAKE;
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'AWAKE');
        end
        
        if isfield(LFPstageavg{ipatient}.(template), "POSTSLEEP")
            LFPstageavg{ipatient}.(template).Post    = LFPstageavg{ipatient}.(template).POSTSLEEP;
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'POSTSLEEP');
        end

        if isfield(LFPstageavg{ipatient}.(template), "PRESLEEP")
            LFPstageavg{ipatient}.(template).Pre    = LFPstageavg{ipatient}.(template).PRESLEEP;
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'PRESLEEP');
        end
        
        if isfield(LFPstageavg{ipatient}.(template), "NO_SCORE")
            LFPstageavg{ipatient}.(template) = rmfield(LFPstageavg{ipatient}.(template), 'NO_SCORE');
        end
    end
end

%% flip polarity to align spike and waves between patients
for ipatient = 6
    for itemplate = 1 : 6
        LFP_cluster{ipatient}{itemplate}.avg = -LFP_cluster{ipatient}{itemplate}.avg;
    end
    for template = ["template1","template2","template3","template4","template5","template6"]
        for hyplabel = ["S3", "S2", "S1", "REM", "Pre", "Post", "Wake"]
            if isfield(LFPstageavg{ipatient}.(template), hyplabel)
                
                LFPstageavg{ipatient}.(template).(hyplabel).trial{1} = -LFPstageavg{ipatient}.(template).(hyplabel).trial{1};
                for itrial = 1 : size(LFPstage{ipatient}.(template).trial, 2)
                    LFPstage{ipatient}.(template).trial{itrial} = -LFPstage{ipatient}.(template).trial{itrial};
                end
            end
        end
        
    end
    for hyplabel = ["S3", "S2", "S1", "REM", "Pre", "Post", "Wake"]
        LFPstage_overtemplate_avg{ipatient}.(hyplabel).trial{1} = -LFPstage_overtemplate_avg{ipatient}.(hyplabel).trial{1};
    end
end

%% determine time periods based on direction of deflection
clear tpos tneg
for ipatient = 1 : 8
    for itemplate = 1 : 6
        maxchan = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');
        
        m1 = max(LFP_cluster{ipatient}{itemplate}.avg(maxchan, :)) / 2;
        m2 = min(LFP_cluster{ipatient}{itemplate}.avg(maxchan, :)) / 2;
        
        % add edges in case timecourse starts or ends above half-height
        t1 = [0, LFP_cluster{ipatient}{itemplate}.avg(maxchan, :), 0] > m1;
        t2 = [0, LFP_cluster{ipatient}{itemplate}.avg(maxchan, :), 0] < m2;
        
        % adjust for start and end edge
        d1 = find(diff(t1)) - 1;
        d1(d1 < 1) = 1;
        d1(d1 > size(LFP_cluster{ipatient}{itemplate}.time, 2)) = size(LFP_cluster{ipatient}{itemplate}.time, 2);
        d2 = find(diff(t2)) - 1;
        d2(d2 < 1) = 1;
        d2(d2 > size(LFP_cluster{ipatient}{itemplate}.time, 2)) = size(LFP_cluster{ipatient}{itemplate}.time, 2);
        
        posindx{ipatient}{itemplate} = zeros(1, size(LFP_cluster{ipatient}{itemplate}.time, 2));
        negindx{ipatient}{itemplate} = zeros(1, size(LFP_cluster{ipatient}{itemplate}.time, 2));
        
        % add timing
        ii = 1;
        for i = 1 : 2 : size(d1 == 1, 2)
            if LFP_cluster{ipatient}{itemplate}.time(d1(i)) <= 0.15 && LFP_cluster{ipatient}{itemplate}.time(d1(i)) > -0.15
                tpos{ipatient}{itemplate}{ii} = LFP_cluster{ipatient}{itemplate}.time([d1(i), d1(i+1)]);
                posindx{ipatient}{itemplate}(d1(i):d1(i+1)) = 1;
                ii = ii + 1;                
            end
        end
        ii = 1;
        for i = 1 : 2 : size(d2 == 1, 2)
            if LFP_cluster{ipatient}{itemplate}.time(d2(i)) <= 0.05 && LFP_cluster{ipatient}{itemplate}.time(d2(i)) > -0.15
                tneg{ipatient}{itemplate}{ii} = LFP_cluster{ipatient}{itemplate}.time([d2(i), d2(i+1)]);
                negindx{ipatient}{itemplate}(d2(i):d2(i+1)) = 1;
                ii = ii + 1;
            end
        end
    end
end

save(strcat(config{ipatient}.datasavedir, "LFP_posneg"),'posindx','negindx');

%% plot figure template timings
nrow = 8;
ncol = 6;
ylims = [2000, 1000, 800, 2000, 4000, 800, 200, 1200] / 1000; % go to mV

for showrejected = [false]
    
    fig = figure('visible', true);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'position', get(0,'ScreenSize'));
    set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
    set(fig, 'position', [10 10 1000/sqrt(2) 1000]);
    set(fig, 'PaperOrientation', 'portrait');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer', 'Painters');
    
    for ipatient = 1:8
        
        col = 1;
        for itemplate = 1:6

            if ~showrejected && any(itemplate == config{ipatient}.template.rejected)
                continue
            end
            
            subplot(nrow, ncol, col + (ipatient-1) * ncol);
%             set(gca, 'clipping', 'off');
            hold on
            
            for i = 1 : size(tpos{ipatient}{itemplate}, 2)
                patch([tpos{ipatient}{itemplate}{i}, ...
                    tpos{ipatient}{itemplate}{i}(2:-1:1)], ...
                    [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                    [0, 0, 1], 'facealpha', 0.1, 'edgecolor', 'none');
            end
            for i = 1 : size(tneg{ipatient}{itemplate}, 2)
                patch([tneg{ipatient}{itemplate}{i}, ...
                    tneg{ipatient}{itemplate}{i}(2:-1:1)], ...
                    [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                    [1, 0, 0], 'facealpha', 0.1, 'edgecolor', 'none');
            end
            
            maxchan = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');

%             plot(LFP_cluster{ipatient}{itemplate}.time, LFP_cluster{ipatient}{itemplate}.avg/1000, 'color', [0.8, 0.8, 0.8]);
            plot(LFP_cluster{ipatient}{itemplate}.time, LFP_cluster{ipatient}{itemplate}.avg(maxchan, :)'/1000, 'color', [0, 0, 0], 'linewidth', 1);
            plot(latency, [0, 0], ':k');
            
            ylim([-ylims(ipatient), ylims(ipatient)]);
            yticks([-ylims(ipatient), ylims(ipatient)]);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0.03]);
                
            % for clinical montage flip Y-axis
            set(gca, 'YDir','reverse')
    
            if col == 1 && ipatient == 8
                xlabel('Time (s)');
                lh = ylabel('mV');
                set(lh, 'units', 'normalized')
                get(lh, 'position');
                lh.Position(1) = -0.1; % change horizontal position of ylabel
                lh.Position(2) =  0.5; % change vertical position of ylabel
            end
            if col ~= 1
                set(gca,'yticklabel', []);
            end
            
            xlim(latency);
            xticks([latency(1), 0, latency(2)]);
            
            if ipatient ~= 8
                set(gca, 'xticklabel', [], 'XColor', 'none');
            end
            
            col = col + 1;
        end
        
        % plot for legend
        if ipatient == 1
            ic = 1;
            clear h
            h(1) = patch([0, 0], [0,0], [0, 0, 1], 'facealpha', 0.1, 'edgecolor', 'none');
            h(2) = patch([0, 0], [0,0], [1, 0, 0], 'facealpha', 0.1, 'edgecolor', 'none');
            spos = get(gca, 'Position');
            l = legend(h, {'Positive', 'Negative'}, 'box', 'off', 'location', 'eastoutside');
            set(gca, 'Position', spos);
        end
        
        set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 10);
        
    end
    
    % write to figure for article
    if ~ispc
        if showrejected
            fname = fullfile(config{1}.imagesavedir, 'template_timing_rejected');
        else
            fname = fullfile(config{1}.imagesavedir, 'template_timing');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
    else
        if showrejected
            fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'template_timing_rejected');
        else
            fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'template_timing');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
    end
    
end

%% plot figure sleep stages
hyplabels = ["S3", "S2", "S1", "REM", "Wake"];
cm = cbrewer('qual', 'Set2', 8);
ylims = [2000, 1000, 500, 900, 2800, 500, 200, 1000] / 1000; % go to mV
nrow = 8;
ncol = 6;
    
for showrejected = [false, true]
    
    fig = figure('visible', true);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'position', get(0,'ScreenSize'));
    set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
    set(fig, 'position', [10 10 1000/sqrt(2) 1000]);
    set(fig, 'PaperOrientation', 'portrait');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer', 'Painters');

    for ipatient = 1:8
        
        col = 1;
        for itemplate = 1:6
                        
            if ~isfield(LFPstageavg{ipatient}, sprintf('template%d', itemplate))
                col = col + 1;
                continue
            end
            if ~showrejected && any(itemplate == config{ipatient}.template.rejected)
                col = col + 1;
                continue
            end

            subplot(nrow, ncol, col + (ipatient-1) * ncol);
            set(gca, 'clipping', 'off');
            hold on
 
            for i = 1 : size(tpos{ipatient}{itemplate}, 2)
                patch([tpos{ipatient}{itemplate}{i}, ...
                       tpos{ipatient}{itemplate}{i}(2:-1:1)], ...
                       [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                       [0, 0, 1], 'facealpha', 0.1, 'edgecolor', 'none');
            end
            for i = 1 : size(tneg{ipatient}{itemplate}, 2)
                patch([tneg{ipatient}{itemplate}{i}, ...
                       tneg{ipatient}{itemplate}{i}(2:-1:1)], ...
                       [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                       [1, 0, 0], 'facealpha', 0.1, 'edgecolor', 'none');
            end     
            
            maxchan     = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');
            plot(LFP_cluster{ipatient}{itemplate}.time, LFP_cluster{ipatient}{itemplate}.avg(maxchan, :)/1000, 'color', [0.7, 0.7, 0.7]);
            plot(latency, [0, 0], ':k');

            for hyplabel = hyplabels(end:-1:1)
                if ~isfield(LFPstageavg{ipatient}.(sprintf('template%d', itemplate)), hyplabel)
                    continue
                end
                maxchan     = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');   
                sem         = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).sem(maxchan, :) / 1000;
                x           = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).time{1};
                y           = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).trial{1}(maxchan, :) / 1000; % go to mV
                patch([x, x(end:-1:1)], [y-sem, y(end:-1:1)+sem(end:-1:1)], cm(9-find(hyplabels == hyplabel), :), 'edgecolor', 'none', 'facealpha', 0.3);
                lh          = plot(x, y, 'color', cm(9-find(hyplabels == hyplabel), :), 'linewidth', 1); 
            end
            
            ylim([-ylims(ipatient), ylims(ipatient)]);
            yticks([-ylims(ipatient), ylims(ipatient)]);
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0.03]);
            
            if col == 1 && ipatient == 8
                xlabel('Time (s)');
                lh = ylabel('mV');
                set(lh, 'units', 'normalized')
                get(lh, 'position');
                lh.Position(1) = -0.1; % change horizontal position of ylabel
                lh.Position(2) = 0.5; % change vertical position of ylabel
            end
            if col ~= 1
                set(gca,'yticklabel', []);
            end
                        
            xlim(latency);
            xticks([latency(1), 0, latency(2)]);
            if ipatient == 8
            else
                set(gca, 'xticklabel', [], 'XColor', 'none');
            end
            
            if any(itemplate == config{ipatient}.template.rejected)
                ax = axis;
                plot([ax(1), ax(2)], [ax(3), ax(4)], 'r', 'linewidth', 2);
            end
 
            col = col + 1;
        end
        
        % plot for legend
        if ipatient == 1
            ic = 1;
            clear h
            for hyplabel = hyplabels
                h(ic) = patch([0, 0], [0,0], cm(9-ic, :), 'facealpha', 1, 'edgecolor', 'none');
                ic = ic + 1;
            end
            
            spos = get(gca, 'Position');
            l = legend(h, hyplabels, 'box', 'off', 'location', 'eastoutside');
            title(l, "Stage");
            l.Title.FontWeight = 'normal';
            set(gca, 'Position', spos);
        end
        
        set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 10);
    end
    
    % write to figure for article %% REPLACE WITH AVERAGE PLOTS
%     if ~ispc
%         if showrejected
%             fname = fullfile(config{1}.imagesavedir, 'LFP_stages_rejected');
%         else
%             fname = fullfile(config{1}.imagesavedir, 'LFP_stages');
%         end
%         isdir_or_mkdir(fileparts(fname));
%         exportgraphics(fig, strcat(fname, '.pdf'));
%     else
%         if showrejected
%             fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'LFP_stages_rejected');
%         else
%             fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'LFP_stages');
%         end
%         isdir_or_mkdir(fileparts(fname));
%         exportgraphics(fig, strcat(fname, '.pdf'));
%     end
end

%% plot average LFPs per sleep stage

fig = figure('visible', true);
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'Renderer', 'Painters');
set(0,'DefaultAxesTitleFontWeight','normal');

hyplabels = ["S3", "S2", "S1", "REM", "Pre", "Post", "Wake"];
hyplabels = ["S3", "S2", "S1", "REM", "Wake"];
cm = cbrewer('qual', 'Set2', 8);
ylims = [600, 1000, 300, 500, 2000, 400, 250, 900] / 1000; % go to mV

for ipatient = 1 : 8
    
    s1 = subplot(2, 4, ipatient);
    hold on
    
    for hyplabel = hyplabels(end:-1:1)
        maxchan     = find(~cellfun(@isempty, strfind(LFPstage_overtemplate_avg{ipatient}.(hyplabel).label, config{ipatient}.align.zerochannel)), 1, 'first');
        sem         = LFPstage_overtemplate_avg{ipatient}.(hyplabel).sem(maxchan, :) / 1000;
        x           = LFPstage_overtemplate_avg{ipatient}.(hyplabel).time{1};
        y           = LFPstage_overtemplate_avg{ipatient}.(hyplabel).trial{1}(maxchan, :) / 1000; % go to mV
        patch([x, x(end:-1:1)], [y-sem, y(end:-1:1)+sem(end:-1:1)], cm(9-find(hyplabels == hyplabel), :), 'edgecolor', 'none', 'facealpha', 0.3);
        lh          = plot(x, y, 'color', cm(9-find(hyplabels == hyplabel), :), 'linewidth', 1);
    end
    
    ylim([-ylims(ipatient), ylims(ipatient)]);
    yticks([-ylims(ipatient), 0, ylims(ipatient)]);
    yticklabels([-ylims(ipatient), "", ylims(ipatient)]);
    set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0.03]);
    
    % for clinical montage flip Y-axis
    set(gca, 'YDir','reverse')
    
    if ipatient == 1 || ipatient == 5
        xlabel('Time (s)');
        lh = ylabel('\muV');
        lh = ylabel('mV');
        set(lh, 'units', 'normalized')
        get(lh, 'position');
        lh.Position(1) = -0.1; % change horizontal position of ylabel
        lh.Position(2) = 0.5; % change vertical position of ylabel
    end
    
    xlim(latency);
    xticks([latency(1), 0, latency(2)]);
    set(gca, 'xticklabel', [], 'XColor', 'none');
    
    % baseline line
    if ipatient == 4
        y = ylim;
        y = y(2);
        x = xlim;
        line([-0.15, -0.05], [y, y], 'linewidth', 2, 'color', 'k');
        text(-0.15, y, 'BL', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        
        % time line
        duration = 0.300;
        t = sprintf('300 ms');
        line([x(2) - duration, config{ipatient}.spike.toi.template1(end)], [y, y], 'linewidth', 2, 'color', 'k');
        text(x(2) * 0.99, y, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
    
    % plot for legend
    if ipatient == 4
        ic = 1;
        clear h
        for hyplabel = hyplabels
            h(ic) = patch([0, 0], [0,0], cm(9-ic, :), 'facealpha', 1, 'edgecolor', 'none');
            ic = ic + 1;
        end
        
        spos = get(gca, 'Position');
        l = legend(h, hyplabels, 'box', 'off', 'location', 'eastoutside');
        title(l, "Stage");
        l.Title.FontWeight = 'normal';
        set(gca, 'Position', spos);
    end
    
    % plot legend (empty to keep horizontal stretch between rows)
    if ipatient == 8
        ic = 1;
        clear h
        for hyplabel = hyplabels
            h(ic) = patch([0, 0], [0,0], cm(9-ic, :), 'facealpha', 0, 'edgecolor', 'none');
            ic = ic + 1;
        end
        
        spos = get(gca, 'Position');
        l = legend(h, ["","","","",""], 'box', 'off', 'location', 'eastoutside');
        title(l, "");
        l.Title.FontWeight = 'normal';
        set(gca, 'Position', spos);
    end
    
    title(sprintf('Patient %d', ipatient));
    col = col + 1;
end
set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 24);

% write to figure for article
if ~ispc
    fname = fullfile(config{1}.imagesavedir, 'LFP_stages');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
else
    fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'LFP_stages');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
end

%% export data to file for R

t = table;
for ipatient = 1 : 8
    for ipart = 1 : 3
        for itemplate = 1 : 6
            fprintf('Extracting data patient %d, part %d, template %d\n', ipatient, ipart, itemplate);
            if ~isfield(LFP{ipatient}{ipart}, sprintf('template%d', itemplate))
                continue
            end
%             if ~showrejected && any(itemplate == config{ipatient}.template.rejected)
%                 continue
%             end
            if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
                continue
            end
            t_temp          = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trialinfo;
            t_temp.patient  = ones(height(t_temp), 1) * ipatient;
            t_temp.part     = ones(height(t_temp), 1) * ipart;
            t_temp.template = ones(height(t_temp), 1) * itemplate;
            maxchan         = find(~cellfun(@isempty, strfind(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, config{ipatient}.align.zerochannel)), 1, 'first');
            
            cfgtemp         = [];
            cfgtemp.channel = maxchan;
            temp            = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
            d               = vertcat(temp.trial{:});
            t_temp.posamp   = mean(d(:, posindx{ipatient}{itemplate} == 1), 2);
            t_temp.negamp   = mean(d(:, negindx{ipatient}{itemplate} == 1), 2);
            
            % add trialinfo
            
            % flip polarity for patient 6
            if ipatient == 6
                t_temp.posamp = -t_temp.posamp;
                t_temp.negamp = -t_temp.negamp;
            end
            
            t = [t; t_temp];
            clear d temp t_temp
        end
    end
end

% save data to table for R
fname   = fullfile(config{ipatient}.datasavedir, 'amplitude_table');
writetable(t, fname);

disp('done');