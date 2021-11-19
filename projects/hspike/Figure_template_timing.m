function [config] = Figure_template_timing

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
MuseStruct{8}           = [];
LFP_cluster_detected{8} = [];

%% load data
config = hspike_setparams;
for ipatient = 1:8
    
    % extracted data
    config{ipatient}.LFP.name   = {'template1','template2','template3','template4','template5','template6'};
    LFP{ipatient}               = readLFP(config{ipatient});
    
    % clusters
    [clusterindx{ipatient}, LFP_cluster{ipatient}] = clusterLFP(config{ipatient});
    LFP_cluster{ipatient} = LFP_cluster{ipatient}{1}.Hspike.kmedoids{6};    
end

%% select latency
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
end
for ipatient = 1 : 8
    for itemplate = 1:6
        LFP_cluster{ipatient}{itemplate} = ft_selectdata(cfgtemp, LFP_cluster{ipatient}{itemplate});
    end
end

%% rereference and baseline data
for ipatient = 1:8
    if strcmp(config{ipatient}.template.reref, 'yes')
        for ipart = 1 : size(LFP{ipatient}, 2)
            for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
                
                if isempty(LFP{ipatient}{ipart}.(markername))
                    continue
                end
                if ~contains(LFP{ipatient}{ipart}.(markername).label{1}, '-')
                    
                    labels_nonum    = regexprep(LFP{ipatient}{ipart}.(markername).label, '[0-9_]', '');
                    [~, ~, indx]    = unique(labels_nonum);
                    
                    % average per part (day) then reref
                    clear group
                    for i = 1 : max(indx)
                        cfgtemp             = [];
                        cfgtemp.reref       = 'yes';
                        cfgtemp.refmethod   = 'bipolar';
                        cfgtemp.demean      = 'yes';
                        %                     cfgtemp.baselinewindow = config{ipatient}.LFP.baselinewindow.(markername);
                        cfgtemp.baselinewindow = [-0.3, -0.1];
                        cfgtemp.channel     = LFP{ipatient}{ipart}.(markername).label(indx==i);
                        group{i}            = ft_preprocessing(cfgtemp, LFP{ipatient}{ipart}.(markername));
                    end
                    LFP{ipatient}{ipart}.(markername) = ft_appenddata([], group{:});
                end
            end
        end
    end
end
    
fname = fullfile(config{1}.datasavedir, 'LFP_reref.mat');
save(fname, 'LFP', '-v7.3');
    
%% rereference and baseline cluster
for ipatient = 1 : 8
    % rereference to bipolar
    if strcmp(config{ipatient}.cluster.reref, 'yes')
        for itemplate = 1 : 6
            cfgtemp             = [];
            cfgtemp.demean      = 'yes';
            cfgtemp.baselinewindow = [-0.3, -0.1];
            LFP_cluster{ipatient}{itemplate} = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
            
            if contains(LFP_cluster{ipatient}{itemplate}.label, '-')
                disp('akready done');
                continue
            end
            disp('rereferencing');
            labels_nonum    = regexprep(LFP_cluster{ipatient}{itemplate}.label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.reref       = 'yes';
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.demean      = 'yes';
                cfgtemp.baselinewindow = [-0.3, -0.1];
                cfgtemp.channel     = LFP_cluster{ipatient}{itemplate}.label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
            end
            LFP_cluster{ipatient}{itemplate} = ft_appenddata([],group{:});
        end
    end
end

%% determine time periods based on direction of deflection
clear tpos tneg
for ipatient = 1 : 8
    for itemplate = 1 : 6
        maxchan     = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)));

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

        ii = 1;
        for i = 1 : 2 : size(d1 == 1, 2)
            tpos{ipatient}{itemplate}{ii} = [d1(i), d1(i+1)];
            ii = ii + 1;
        end
        ii = 1;
        for i = 1 : 2 : size(d2 == 1, 2)
            tneg{ipatient}{itemplate}{ii} = [d2(i), d2(i+1)];
            ii = ii + 1;
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

hyplabels = ["S3", "S2", "S1", "REM", "Pre", "Post", "Wake"];
hyplabels = ["S3", "S2", "S1", "REM", "Wake"];
cm = cbrewer('qual', 'Set2', 7);

%% plot figure
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
        
    if showrejected
        ncol = 6;
    else
        ncol = 5; % depends on maximum amount of selected templates
    end
    nrow = 8;
    ylims = [1000, 1000, 1000, 1500, 3000, 300, 300, 500] / 1000; % go to mV
    
    for ipatient = 1:8
        
        col = 1;
        for itemplate = 6:-1:1
            
            if ~isfield(LFPstageavg{ipatient}, sprintf('template%d', itemplate))
                continue
            end
            if ~showrejected && any(itemplate == config{ipatient}.template.rejected)
                continue
            end
            
            subplot(nrow, ncol, col + (ipatient-1) * ncol);
            set(gca, 'clipping', 'off');
            hold on
            
            for i = 1 : size(tpos{ipatient}{itemplate}, 2)
                patch([LFP_cluster{ipatient}{itemplate}.time(tpos{ipatient}{itemplate}{i}), ...
                       LFP_cluster{ipatient}{itemplate}.time(tpos{ipatient}{itemplate}{i}(2:-1:1))], ...
                       [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                       [0, 1, 0], 'facealpha', 0.1, 'edgecolor', 'none');
            end
            for i = 1 : size(tneg{ipatient}{itemplate}, 2)
                patch([LFP_cluster{ipatient}{itemplate}.time(tneg{ipatient}{itemplate}{i}), ...
                       LFP_cluster{ipatient}{itemplate}.time(tneg{ipatient}{itemplate}{i}(2:-1:1))], ...
                       [-ylims(ipatient), -ylims(ipatient), ylims(ipatient), ylims(ipatient)], ...
                       [1, 0, 0], 'facealpha', 0.1, 'edgecolor', 'none');
            end    
            
            maxchan     = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)));
            plot(LFP_cluster{ipatient}{itemplate}.time, LFP_cluster{ipatient}{itemplate}.avg(maxchan, :)/1000, 'color', [0.5, 0.5, 0.5]);
            plot(latency, [0, 0], ':k');
            

            for hyplabel = hyplabels
                if ~isfield(LFPstageavg{ipatient}.(sprintf('template%d', itemplate)), hyplabel)
                    continue
                end
                                
                maxchan     = ~cellfun(@isempty, strfind(LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).label, config{ipatient}.align.zerochannel));               
                sem         = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).sem(maxchan, :) / 1000;
                x           = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).time{1};
                y           = LFPstageavg{ipatient}.(sprintf('template%d', itemplate)).(hyplabel).trial{1}(maxchan, :) / 1000; % go to mV
                patch([x, x(end:-1:1)], [y-sem, y(end:-1:1)+sem(end:-1:1)], cm(hyplabels == hyplabel, :), 'edgecolor', 'none', 'facealpha', 0.3);
                lh          = plot(x, y, 'color', cm(hyplabels == hyplabel, :), 'linewidth', 1); 
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
            %             if (col == 1 && ipatient == 8)
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
                h(ic) = patch([0, 0], [0,0], cm(ic, :), 'facealpha', 1, 'edgecolor', 'none');
                ic = ic + 1;
            end
            
            spos = get(gca, 'Position');
            l = legend(h, hyplabels, 'box', 'off', 'location', 'eastoutside');
            title(l, "Stage");
            set(gca, 'Position', spos);
        end
        
        set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 10);
    end
    
    % write to figure for article (if not on Desktop PC)
    if ~ispc
        if showrejected
            fname = fullfile(config{1}.imagesavedir, 'LFP_stages_rejected');
        else
            fname = fullfile(config{1}.imagesavedir, 'LFP_stages');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
    else
        % write to figure for Article
        if showrejected
            fname = fullfile('D:\Dropbox\Apps\Overleaf\images\Hspike', 'LFP_stages_rejected');
        else
            fname = fullfile('D:\Dropbox\Apps\Overleaf\images\Hspike', 'LFP_stages');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
    end
    
end




    
    
