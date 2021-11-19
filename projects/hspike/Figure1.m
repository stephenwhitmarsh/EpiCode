function Figure1

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
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
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end

ft_defaults
MuseStruct{8}           = [];
LFP_cluster{8}          = [];
LFP_cluster_detected{8} = [];
clusterindx{8}          = [];
ipart = 1;

%% Why are templates missing in part 1

for ipatient = 1:8
    config                                                      = hspike_setparams;
    [clusterindx{ipatient}, LFP_cluster{ipatient}]              = clusterLFP(config{ipatient});
%     [config{ipatient}, LFP_cluster{ipatient}]                   = alignClusters(config{ipatient},  LFP_cluster{ipatient}{1}.Hspike.kmedoids{6});
    LFP_cluster{ipatient}                   = LFP_cluster{ipatient}{1}.Hspike.kmedoids{6};
    
    [MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}]   = detectTemplate(config{ipatient});
    
    
    % now get some single event data (without overwriting)
    config{ipatient}.LFP.name                                   = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
    LFP{ipatient}                                               = readLFP(config{ipatient});
    
    % rereference to bipolar
    if strcmp(config{ipatient}.template.reref, 'yes')
        for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
            if isempty(LFP{ipatient}{ipart}.(markername))
                continue
            end
            if contains(LFP{ipatient}{ipart}.(markername).label{1}, '-')
                continue
            end
            labels_nonum    = regexprep(LFP{ipatient}{ipart}.(markername).label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.reref       = 'yes';
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.demean      = 'yes';
                cfgtemp.channel     = LFP{ipatient}{ipart}.(markername).label(indx==i);
                cfgtemp.trials      = randperm(size(LFP{ipatient}{ipart}.(markername).trial, 2), min(100, size(LFP{ipatient}{ipart}.(markername).trial, 2)));
                group{i}            = ft_preprocessing(cfgtemp,LFP{ipatient}{ipart}.(markername));
            end
            LFP{ipatient}{ipart}.(markername) = ft_appenddata([], group{:});
        end
    end
    
    % rereference to bipolar
    if strcmp(config{ipatient}.cluster.reref, 'yes')
        for itemplate = 1 : 6
            if contains(LFP_cluster{ipatient}{itemplate}.label, '-')
                continue
            end
            labels_nonum    = regexprep(LFP_cluster{ipatient}{itemplate}.label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.reref       = 'yes';
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.demean      = 'yes';
                cfgtemp.channel     = LFP_cluster{ipatient}{itemplate}.label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
            end
            LFP_cluster{ipatient}{itemp 297,39late} = ft_appenddata([],group{:});
        end
    end
    
end

% plot figure

fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
set(fig, 'position', [10 10 1000/sqrt(2) 1000]);
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

for ipatient = 1:7
    
%     labelled = false;
%     maxrange = 0;
%     for itemplate = 1 : 6
%         if isfield(LFP_cluster{ipatient}{itemplate}, 'avg')
%             [~, maxchan] = max(max(LFP_cluster{ipatient}{itemplate}.avg'));
%         else
%             [~, maxchan] = max(max(LFP_cluster{ipatient}{itemplate}.trial{1}'));
%         end
%         if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%             continue
%         end
%         
%         cfgtemp = [];
%         cfgtemp.channel = maxchan; 297,39
%         temp = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
%         temp = vertcat(temp.trial{:});
%         maxrange = max(max(max(abs(temp))), maxrange);
%     end
    
    for itemplate = 1 : 6
        
        if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
            continue
        end
%         if isfield(LFP_cluster{ipatient}{itemplate}, 'avg')
%             [~, maxchan] = max(max(LFP_cluster{ipatient}{itemplate}.avg'));
%         else
%             [~, maxchan] = max(max(LFP_cluster{ipatient}{itemplate}.trial{1}'));
%         end
        maxchan = find(ismember(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel));
        
%         subplot(8, 6, itemplate + (ipatient-1) * 6);
        figure
        hold;
        
        for itrial = 1 : min(100, size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial, 2))
            
            x           = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{itrial};
            y           = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{itrial}(maxchan, :);
            lh          = plot(x, y, 'k');
            lh.Color    = [lh.Color 0.1];
            label       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label{maxchan};
            n           = n + 1;
        end
        
        cfgtemp             = [];
        cfgtemp.channel     = maxchan;
        cfgtemp.avgoverrpt  = 'yes';
        temp                = ft_selectdata(cfgtemp, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
        plot(x, temp.trial{1}, 'b', 'linewidth', 1);
        
        n = 1;
        if isfield(LFP_cluster{ipatient}{itemplate}, 'avg')
            y = LFP_cluster{ipatient}{itemplate}.avg(maxchan, :);
            x = LFP_cluster{ipatient}{itemplate}.time;
            
        else
            y = LFP_cluster{ipatient}{itemplate}.trial{1}(maxchan, :);
            x = LFP_cluster{ipatient}{itemplate}.time{1};
        end
        
        plot(x, y, 'r', 'linewidth', 1);
        
        xlim(config{ipatient}.cluster.latency);
%         ylim([-maxrange, maxrange]);
        
%         if labelled == false
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', {'', label, ''}, 'TickDir', 'out')
%             labelled = true;
%         else
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', [], 'TickDir', 'out')
%         end
        
    end
end






% 
% 
% 
% % plot figure
% 
% fig = figure('visible', true);
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
% set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
% set(fig, 'PaperOrientation', 'portrait');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0 0 1 1]);
% set(fig, 'Renderer', 'Painters');
% 
% for ipatient = 1:7
%         
%     labelled = false;
%     maxrange = 0;
%     for itemplate = 1 : 6
%         
%         if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%             continue
%         end
%         
%         temp = vertcat(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{:});
%         maxrange = max(max(max(abs(temp))), maxrange);
%     end
%     maxrange = maxrange / 4;
%     
%     for itemplate = 1 : 6
%         
%         if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%             continue
%         end
%         
%         subplot(8, 6, itemplate + (ipatient-1) * 6);
%         if labelled == false
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', label, 'TickDir', 'out')
%             labelled = true;
%         else
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', [], 'TickDir', 'out')
%         end
%         
%         hold;        
%         
%         for itrial = 1 : min(10, size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial, 2))
%             n = 1; ytick = []; label = [];
%             for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
%                 ytick = [ytick, n*maxrange];
%                 x       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{itrial};
%                 y       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{itrial}(ichan, :);
%                 lh = plot(x, y + n*maxrange, 'k');
%                 lh.% 
% 
% 
% % plot figure
% 
% fig = figure('visible', true);
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
% set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
% set(fig, 'PaperOrientation', 'portrait');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0 0 1 1]);
% set(fig, 'Renderer', 'Painters');
% 
% for ipatient = 1:7
%         
%     labelled = false;
%     maxrange = 0;
%     for itemplate = 1 : 6
%         
%         if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%             continue
%         end
%         
%         temp = vertcat(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{:});
%         maxrange = max(max(max(abs(temp))), maxrange);
%     end
%     maxrange = maxrange / 4;
%     
%     for itemplate = 1 : 6
%         
%         if isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%             continue
%         end
%         
%         subplot(8, 6, itemplate + (ipatient-1) * 6);
%         if labelled == false
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', label, 'TickDir', 'out')
%             labelled = true;
%         else
%             set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', [], 'TickDir', 'out')
%         end
%         
%         hold;        
%         
%         for itrial = 1 : min(10, size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial, 2))
%             n = 1; ytick = []; label = [];
%             for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
%                 ytick = [ytick, n*maxrange];
%                 x       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{itrial};
%                 y       = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{itrial}(ichan, :);
%                 lh = plot(x, y + n*maxrange, 'k');
%                 lh.Color = [lh.Color 0.1];
%                 label{ichan} = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label{ichan};
%                 n = n + 1;
%             end
%         end
%         n = 1;
%         for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
%             x       = LFP_cluster{ipatient}{itemplate}.time;
%             y       = LFP_cluster{ipatient}{itemplate}.avg(ichan, :);
%             plot(x, mean(LFP_cluster{ipatient}{itemplate}.avg(ichan, :)) + n*maxrange, 'b', 'linewidth', 1);
%             plot(x, y + n*maxrange, 'r', 'linewidth', 1);
%             n = n + 1;
%         end
%         xlim([LFP_cluster{ipatient}{itemplate}.time(1), LFP_cluster{ipatient}{itemplate}.time(end)]);
%         ylim([n/2, maxrange * n + n/2]);
%     end
% endColor = [lh.Color 0.1];
%                 label{ichan} = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label{ichan};
%                 n = n + 1;
%             end
%         end
%         n = 1;
%         for ichan = 1 : size(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, 1)
%             x       = LFP_cluster{ipatient}{itemplate}.time;
%             y       = LFP_cluster{ipatient}{itemplate}.avg(ichan, :);
%             plot(x, mean(LFP_cluster{ipatient}{itemplate}.avg(ichan, :)) + n*maxrange, 'b', 'linewidth', 1);
%             plot(x, y + n*maxrange, 'r', 'linewidth', 1);
%             n = n + 1;
%         end
%         xlim([LFP_cluster{ipatient}{itemplate}.time(1), LFP_cluster{ipatient}{itemplate}.time(end)]);
%         ylim([n/2, maxrange * n + n/2]);
%     end
% end