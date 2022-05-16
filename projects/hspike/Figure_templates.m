function [config] = Figure_templates

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
nr_examples = 100;

%% load data
for ipatient = 1:8
    config = hspike_setparams;
    config{ipatient}                            = addparts(config{ipatient});
    
    [clusterindx{ipatient}, trialindx{ipatient}, LFP_cluster{ipatient}] = clusterLFP(config{ipatient});
    LFP_cluster{ipatient} = LFP_cluster{ipatient}{1}.Hspike.kmedoids{6};

    [MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}] = detectTemplate(config{ipatient});
    
    % get the averages
    config{ipatient}.LFP.name = {'template1', 'template2', 'template3', 'template4', 'template5', 'template6'};
    LFPavg{ipatient} = readLFPavg(config{ipatient});
    
    % get LFP of the original  annotations
    config{ipatient}.LFP.name = {'Hspike'};
    LFP_original{ipatient} = readLFP(config{ipatient});

    % read timestamps of aligned Hspike events
    MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);

end

%% select random examples
for ipatient = 1:8
    cfgtemp                 = [];
    cfgtemp.trials          = randperm(size(LFP_original{ipatient}{1}.Hspike.trial, 2), min(nr_examples, size(LFP_original{ipatient}{1}.Hspike.trial, 2)));
    LFP_original_examples{ipatient}{1}.Hspike = ft_selectdata(cfgtemp, LFP_original{ipatient}{1}.Hspike);
end

%% rereference random examples
for ipatient = 1:8

    if ~contains(LFP_original_examples{ipatient}{1}.Hspike.label{1}, '-')
        labels_nonum    = regexprep(LFP_original_examples{ipatient}{1}.Hspike.label, '[0-9_]', '');
        [~, ~, indx]    = unique(labels_nonum);
        clear group
        for i = 1 : max(indx)
            cfgtemp             = [];
            cfgtemp.reref       = 'yes';
            cfgtemp.refmethod   = 'bipolar';
            cfgtemp.channel     = LFP_original_examples{ipatient}{1}.Hspike.label(indx==i);
            cfgtemp.trials      = randperm(size(LFP_original_examples{ipatient}{1}.Hspike.trial, 2), min(nr_examples, size(LFP_original_examples{ipatient}{1}.Hspike.trial, 2)));
            group{i}            = ft_preprocessing(cfgtemp, LFP_original_examples{ipatient}{1}.Hspike);
        end
        LFP_original_examples{ipatient}{1}.Hspike = ft_appenddata([], group{:});
    end
    
end


%% average original
for ipatient = 1 : 8
    cfgtemp             = [];
    cfgtemp.avgoverrpt  = 'yes';
    LFP_original_avg{ipatient}{1}.Hspike = ft_selectdata(cfgtemp, LFP_original{ipatient}{1}.Hspike);
end


%% rereference original average
for ipatient = 1 : 8
    if ~contains(LFP_original_avg{ipatient}{1}.Hspike.label{1}, '-')
        labels_nonum    = regexprep(LFP_original_examples{ipatient}{1}.Hspike.label, '[0-9_]', '');
        [~, ~, indx]    = unique(labels_nonum);
        clear group
        for i = 1 : max(indx)
            cfgtemp             = [];
            cfgtemp.reref       = 'yes';
            cfgtemp.refmethod   = 'bipolar';
            cfgtemp.channel     = LFP_original_avg{ipatient}{1}.Hspike.label(indx==i);
            group{i}            = ft_preprocessing(cfgtemp, LFP_original_avg{ipatient}{1}.Hspike);
        end
        LFP_original_avg{ipatient}{1}.Hspike = ft_appenddata([], group{:});
    end
end

%% rereference template average
for ipatient = 1:8
    for ipart = 1 : size(LFP_cluster_detected{ipatient}, 2)
        for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
            
            if ~isfield(LFPavg{ipatient}{ipart}, markername)
                continue
            end            
            if isempty(LFPavg{ipatient}{ipart}.(markername))
                continue
            end
            LFPavg_reref{ipatient}{ipart}.(markername) = LFPavg{ipatient}{ipart}.(markername);

            if strcmp(config{ipatient}.template.reref, 'yes') & ~contains(LFPavg{ipatient}{ipart}.(markername).label{1}, '-')
                
                labels_nonum    = regexprep(LFPavg{ipatient}{ipart}.(markername).label, '[0-9_]', '');
                [~, ~, indx]    = unique(labels_nonum);
                
                % average per part (day) then reref
                clear group
                for i = 1 : max(indx)
                    cfgtemp             = [];
                    cfgtemp.reref       = 'yes';
                    cfgtemp.refmethod   = 'bipolar';
                    cfgtemp.channel     = LFPavg{ipatient}{ipart}.(markername).label(indx==i);
                    group{i}            = ft_preprocessing(cfgtemp, LFPavg{ipatient}{ipart}.(markername));
                end
                LFPavg_reref{ipatient}{ipart}.(markername) = ft_appenddata([], group{:});
            end
        end
    end
end

%% rereference cluster
for ipatient = 1 : 8
    % rereference to bipolar
    for itemplate = 1 : 6
        if strcmp(config{ipatient}.cluster.reref, 'yes')
            if isempty(LFP_cluster{ipatient}{itemplate})
                continue
            end
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
                cfgtemp.channel     = LFP_cluster{ipatient}{itemplate}.label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
            end
            LFP_cluster{ipatient}{itemplate} = ft_appenddata([],group{:});
        end
    end
end

%% baseline correct
cfgtemp                 = [];
cfgtemp.demean          = 'yes';
cfgtemp.baselinewindow  = [-0.3, -0.1];

for ipatient = 1 : 8
    LFP_original_examples{ipatient}{1}.Hspike = ft_preprocessing(cfgtemp, LFP_original_examples{ipatient}{1}.Hspike);
    LFP_original_avg{ipatient}{1}.Hspike      = ft_preprocessing(cfgtemp, LFP_original_avg{ipatient}{1}.Hspike);

    for ipart = 1 : size(LFP_cluster_detected{ipatient}, 2)
        for markername = ["template1", "template2", "template3", "template4", "template5", "template6"]
           
            if ~isfield(LFPavg{ipatient}{ipart}, markername)
                continue
            end
            if isempty(LFPavg{ipatient}{ipart}.(markername))
                continue
            end
            
            LFPavg_reref{ipatient}{ipart}.(markername) = ft_preprocessing(cfgtemp, LFPavg_reref{ipatient}{ipart}.(markername));
        end

        for itemplate = 1 : 6
            if isempty(LFP_cluster{ipatient}{itemplate})
                continue
            end
            LFP_cluster{ipatient}{itemplate} = ft_preprocessing(cfgtemp, LFP_cluster{ipatient}{itemplate});
        end
    end
end

%% select latency
latency = [-0.2 0.5];
cfgtemp = [];
cfgtemp.latency = latency;
for ipatient = 1 : 8
    LFP_original_examples{ipatient}{1}.Hspike = ft_selectdata(cfgtemp, LFP_original_examples{ipatient}{1}.Hspike);
    for itemplate = 1 : 6
        for ipart = 1 : size(LFPavg_reref{ipatient}, 2)
            if isfield(LFPavg_reref{ipatient}{ipart}, sprintf('template%d', itemplate))
                if ~isempty(LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)))
                LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_selectdata(cfgtemp, LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)));
                end
            end
        end
        LFP_cluster{ipatient}{itemplate} = ft_selectdata(cfgtemp, LFP_cluster{ipatient}{itemplate});
    end
end

%% calculate maximum range
maxrange(8) = 0;
for ipatient = 1:8
    maxrange(ipatient) = 0;
    for itemplate = 1 : 6
        cfgtemp = [];
        cfgtemp.channel = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');
        cfgtemp.latency = [0, 0.5];
        temp = ft_selectdata(cfgtemp, LFP_cluster{ipatient}{itemplate});
        maxrange(ipatient) = max(max(abs(temp.avg)), maxrange(ipatient));
    end
end

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
    
    ncol = 8; % depends on maximum amount of selected templates
    nrow = 8;
    ylims = [1000, 1000, 1000, 1500, 3000, 300, 300, 500] / 1000; % go to mV

    for ipatient = 1:8
        
        
        % originals
        subplot(nrow, ncol, 1 + (ipatient-1) * ncol);
        set(gca, 'clipping', 'on');
        hold on;
        maxchan = find(~cellfun(@isempty, strfind(LFP_original_examples{ipatient}{1}.Hspike.label, config{ipatient}.align.zerochannel)), 1, 'first');
        for itrial = 1 : min(100, size(LFP_original_examples{ipatient}{1}.Hspike.trial, 2))
            x           = LFP_original_examples{ipatient}{1}.Hspike.time{itrial};
            y           = LFP_original_examples{ipatient}{1}.Hspike.trial{itrial}(maxchan, :) / 1000; % got to mV
            lh          = plot(x, y, 'k');
            lh.Color    = [lh.Color 0.05];
        end
        set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.03, 0.03])
        colorbar('off');
        
        xlim(latency);
        if ipatient == 8
            xticks([latency(1), 0, latency(end)]);
            if ipatient == nrow
                xlabel('Time (s)')
            end
        else
            xticks([latency(1), 0, latency(end)]);
            set(gca,'XTickLabels', [])
            %             set(gca,'XColor', 'none');
        end
        
        ylim([-ylims(ipatient), ylims(ipatient)]);
        yticks([-ylims(ipatient), ylims(ipatient)]);
        lh = ylabel('mV');
        set(lh, 'units', 'normalized')
        get(lh, 'position');
        lh.Position(1) = -0.1; % change horizontal position of ylabel
        lh.Position(2) = 0.5; % change vertical position of ylabel
        
        
        % Recovered events
        row = 1;
        for itemplate = 1 : 6
            
            if ~showrejected & any(itemplate == config{ipatient}.template.rejected)
                continue
            end
             
            subplot(nrow, ncol, 1 + row + (ipatient-1) * ncol);
            set(gca, 'clipping', 'off');
            
            hold on
            cm = cool(size(LFP_cluster_detected{ipatient}, 2));
            cm = cm(end:-1:1,:);
            
            for ipart = 1 : size(LFP_cluster_detected{ipatient}, 2)
                if ~isfield(LFPavg_reref{ipatient}{ipart},sprintf('template%d', itemplate))
                    continue
                end
                if isempty(LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)))
                    continue
                end
                maxchan     = find(~cellfun(@isempty, strfind(LFP_cluster{ipatient}{itemplate}.label, config{ipatient}.align.zerochannel)), 1, 'first');
                x           = LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{1};
                y           = LFPavg_reref{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{1}(maxchan, :) / 1000; % go to mV
                lh          = plot(x, y, 'color', cm(ipart, :));
            end
            
            if isfield(LFP_cluster{ipatient}{itemplate}, 'avg')
                y = LFP_cluster{ipatient}{itemplate}.avg(maxchan, :) / 1000; % go to mV
                x = LFP_cluster{ipatient}{itemplate}.time;
            else
                y = LFP_cluster{ipatient}{itemplate}.trial{1}(maxchan, :) / 1000; % go to mV
                x = LFP_cluster{ipatient}{itemplate}.time{1};
            end
            plot(x, y, 'k');
            
            xlim(latency);
            xticks([latency(1), 0, latency(2)]);  
%             ylim([-ylims(ipatient), ylims(ipatient)]);
            yticks([-ylims(ipatient), ylims(ipatient)]);
            
            if any(itemplate == config{ipatient}.template.rejected)
                ax = axis;
                plot([ax(1), ax(2)], [ax(3), ax(4)], 'r', 'linewidth', 2);
            end
                
            set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'yticklabel', [], 'xticklabel', [], 'TickDir', 'out', 'XColor', 'none')
            
            row = row + 1;
            
        end
        
        spos = get(gca, 'Position');
        h = colorbar;
        h.TickDirection = 'out';
        h.TickLength = 0.05;
        h.Box = 'off';
        hl = ylabel(h, 'Day');
        hl.Position = [2, 0.5, 0];
        set(gca, 'Position', spos);
        colormap(cm);      
        h.Ticks = [0, 1];
        h.TickLabels = ["1", num2str(size(cm, 1))];
        h.YDir = 'reverse';
        
    end
    set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 10);
    
    % write to figure for article (if not on Desktop PC)
    if ~ispc
        if showrejected
            fname = fullfile(config{1}.imagesavedir, 'templates_overview_rejected');
        else
            fname = fullfile(config{1}.imagesavedir, 'templates_overview');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
    else
        % write to figure for Article
        if showrejected
            fname = fullfile('D:\Dropbox\Apps\Overleaf\images\Hspike', 'templates_overview_rejected');
        else
            fname = fullfile('D:\Dropbox\Apps\Overleaf\images\Hspike', 'templates_overview');
        end
        isdir_or_mkdir(fileparts(fname));
        exportgraphics(fig, strcat(fname, '.pdf'));
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