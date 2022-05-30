function Figure_raster

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

config = hspike_setparams;

for ipatient = 1 : 8
    SpikeRaw{ipatient}           = readSpikeRaw_Phy(config{ipatient});
    config{ipatient}.spike.name  = ["template1", "template2", "template3", "template4", "template5", "template6"];
    SpikeTrials{ipatient}        = readSpikeTrials(config{ipatient});
    SpikeStats{ipatient}         = spikeTrialStats(config{ipatient});
    config{ipatient}.LFP.name    = ["template1", "template2", "template3", "template4", "template5", "template6"];
    LFPavg{ipatient}             = readLFPavg(config{ipatient});
    LFP{ipatient}                = rerefLFP(config{ipatient});    
end

% 
% for ipatient = 1 : 8
%     
%     for itemplate = 1 : 6
%         for ipart = 1 : 3
%             i = 1;
%             clear temp
%             if ~isempty(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)))
%                 temp{i} = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate));
%                 i = i + 1;
%             end
%             
%             if exist('temp','var')
%                 cfg = [];
%                 LFPavg_temp{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_appenddata(cfg, temp{:});
%                 
%                 cfg = [];
%                 cfg.avgoverrpt = 'yes';
%                 LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_selectdata(cfg, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
%                 if strcmp(config{ipatient}.cluster.reref, 'yes')
%                     disp('Rereferencing');
%                     labels_nonum    = regexprep(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, '[0-9_]', '');
%                     [~,~,indx]      = unique(labels_nonum);
%                     clear group
%                     for i = 1 : max(indx)
%                         cfgtemp             = [];
%                         cfgtemp.reref       = 'yes';
%                         cfgtemp.refmethod   = 'bipolar';
%                         cfgtemp.demean      = 'yes';
%                         cfgtemp.baselinewindow = [-0.3, -0.1];
%                         cfgtemp.channel     = LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).label(indx==i);
%                         group{i}            = ft_preprocessing(cfgtemp, LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)));
%                     end
%                     LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_appenddata([], group{:});
%                 end
%             end
%         end
%     end
% end




% rereference average LFP to biolar

for ipatient = 1 : 8
    
    if strcmp(config{ipatient}.cluster.reref, 'yes')
        for ipart = 1 : 3
            for itemplate = 1 : 6
                
                if ~isfield(LFPavg{ipatient}{ipart}, sprintf('template%d', itemplate))
                    continue
                end
                if isempty(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)))
                    continue
                end
                if contains(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, '-')
                    disp('Rereferencing already done');
                    continue
                end
                disp('Rereferencing');
                labels_nonum    = regexprep(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, '[0-9_]', '');
                [~,~,indx]      = unique(labels_nonum);
                clear group
                for i = 1 : max(indx)
                    cfgtemp             = [];
                    cfgtemp.reref       = 'yes';
                    cfgtemp.refmethod   = 'bipolar';
                    cfgtemp.demean      = 'yes';
                    cfgtemp.baselinewindow = [-0.3, -0.1];
                    cfgtemp.channel     = LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).label(indx==i);
                    group{i}            = ft_preprocessing(cfgtemp, LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)));
                end
                LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_appenddata([], group{:});
            end
        end
    end
end

%% append spikes

for ipatient = 1 : 8
    
    if ~isempty(SpikeTrials{ipatient}{ipart}.template1)
        SpikeTrials{ipatient}{ipart}.template0 = SpikeTrials{ipatient}{ipart}.template1;
        startindx = 2;
    else
        SpikeTrials{ipatient}{ipart}.template0 = SpikeTrials{ipatient}{ipart}.template2;
        startindx = 3;
    end

    for itemplate = startindx : 6
        if isempty(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)))
            continue
        end
        if isempty(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)))
            continue
        end
        
        for ispike = 1 : size(SpikeTrials{ipatient}{ipart}.template0.label, 2)
            for fn = ["amplitude", "sample", "time", "timestamp", "trial"]
%                 SpikeTrials{ipatient}{ipart}.template0.(fn){1} = {};
                SpikeTrials{ipatient}{ipart}.template0.(fn){1} = [SpikeTrials{ipatient}{ipart}.template0.(fn){1}, SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).(fn){ispike}];
            end
        end
%         for fn = ["trial","sample"]
%             SpikeTrials{}{ipart}.template0.(fn) = [SpikeTrials{ipatient}{ipart}.template0.(fn), SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).(fn)];
%             SpikeTrials{ipatient}{ipart}.template0.(fn) = [SpikeTrials{ipatient}{ipart}.template0.(fn), SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).(fn)];
%         end
    end
end
%%



% configure figure
fig = figure('visible', true);
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
% set(fig, 'PaperOrientation', 'portrait');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0 0 1 1]);
% set(fig, 'Renderer', 'Painters');

% parameters for subplots
ncols       = 6;
nrows       = 8 * 2;
w           = 0.90/ncols;
hratio      = 0.98/ncols;
vratio      = 1/nrows;
htop        = 1/nrows * 0.8;
hbottom     = 1/nrows * 0.8;
spacehor    = 0.3;
rshift      = 0;
upshift     = -0.01;
imarker     = 1;
xlimits     = [-0.2, 0.8];
rownr       = 0;
ylims = [2000, 1000, 1500, 1500, 3000, 1000, 300, 500]; % go to mV

for ipatient = 1 : 8
    rownr = rownr + 1;
    
    itemplate = 1;
    
    % LFP
    row     = (rownr-1)*2 + 1;
    
%     if ~isfield(LFPavg{ipatient}{ipart}, 'template0')
%         continue
%     end
%     if isempty(LFPavg{ipatient}{ipart}.template0)
%         continue
%     end
    %         maxchan = find(~cellfun(@isempty, strfind(LFPavg{ipatient}{ipart}.template0.label, config{ipatient}.align.zerochannel)), 1, 'first');
    %         s1      = axes('Position', [hratio*(itemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
    %
    %         plot(LFPavg{ipatient}{ipart}.template0.time{1}, LFPavg{ipatient}{ipart}.template0.trial{1}(maxchan, :), 'k');
    %         set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
    %         xlim(xlimits);
    %         xticks([xlimits(1), 0, xlimits(2)])
    %         ylim([-ylims(ipatient), ylims(ipatient)]);
    %         yticks([-ceil(max(abs(y))), 0, ceil(max(abs(y)))]);
    %
    %         if itemplate == 1
    %             l1 = ylabel('Amplitude (μV)');
    %             set(l1, 'Units', 'normalized');
    %             lpos = l1.Position;
    %             set(l1, 'Position', [-0.1, 0.5, 0]);
    %         end
    
%     next row
    row     = (rownr-1)*2 + 2;
    s2      = axes('Position', [hratio*(ipatient-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
    set(s2, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
    
    if isempty(SpikeTrials{ipatient}{ipart}.template0)
        continue
    end

    % raster
    cfg_raster              = [];
    cfg_raster.spikechannel = 1;
    cfg_raster.trialborders = 'no';
    %     cmap = repmat([linspace(0, 0.8, size(SpikeTrials{ipatient}{ipart}.SEIZURE.label, 2))], 3, 1)';
    %     cfg_raster.cmapneurons  = cmap;
    %         cfg_raster.cmapneurons  = gray(size(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).label, 2));
    ft_spike_plot_raster(cfg_raster, SpikeTrials{ipatient}{ipart}.template0);
    
    xlim(xlimits);
    xticks([xlimits(1), 0, xlimits(2)])
    y = ylim;
    ylim([1, ceil(max(abs(y)))]);
    yticks([1, ceil(max(abs(y)))]);
    set(gca,'ticklength', [0.005 0.100]);
    
    if itemplate == 1
        l1 = ylabel('Trial');
        set(l1, 'Units', 'normalized');
        lpos = l1.Position;
        set(l1, 'Position', [-0.1, 0.5, 0]);
    else
        ylabel([]);
    end
    
end

% 
% 
% % write to figure
% fname = fullfile(config{1}.imagesavedir, 'article', 'rasterplots_templates');
% exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
% exportgraphics(fig, strcat(fname, '.pdf'));


%%

% configure figure
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

% parameters for subplots
ncols       = 6;
nrows       = 8 * 2;
w           = 0.90/ncols;
hratio      = 0.98/ncols;
vratio      = 1/nrows;
htop        = 1/nrows * 0.8;
hbottom     = 1/nrows * 0.8;
spacehor    = 0.3;
rshift      = 0;
upshift     = -0.01;
imarker     = 1;
xlimits     = [-0.2, 0.8];
rownr       = 0;
ylims = [2000, 1000, 1500, 1500, 3000, 1000, 300, 500]; % go to mV

for ipatient = 1 : 8
    rownr = rownr + 1;

    for itemplate = 1 : 6
        
        
        % LFP
        row     = (rownr-1)*2 + 1;
        
        if ~isfield(LFPavg{ipatient}{ipart}, sprintf('template%d', itemplate))
            continue
        end
        if isempty(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)))
            continue
        end      
        maxchan = find(~cellfun(@isempty, strfind(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, config{ipatient}.align.zerochannel)), 1, 'first');
        s1      = axes('Position', [hratio*(itemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
        
        plot(LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).time{1}, LFPavg{ipatient}{ipart}.(sprintf('template%d', itemplate)).trial{1}(maxchan, :), 'k');
        set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
        xlim(xlimits);
        xticks([xlimits(1), 0, xlimits(2)])
        ylim([-ylims(ipatient), ylims(ipatient)]);
        yticks([-ceil(max(abs(y))), 0, ceil(max(abs(y)))]);
        
        if itemplate == 1
            l1 = ylabel('Amplitude (μV)');
            set(l1, 'Units', 'normalized');
            lpos = l1.Position;
            set(l1, 'Position', [-0.1, 0.5, 0]);
        end
 
        % next row
        row     = (rownr-1)*2 + 2;
        s2      = axes('Position', [hratio*(itemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
        set(s2, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);

        if isempty(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)))
            continue
        end
        % raster
        cfg_raster              = [];
        cfg_raster.spikechannel = 1;
        cfg_raster.trialborders = 'no';
        %     cmap = repmat([linspace(0, 0.8, size(SpikeTrials{ipatient}{ipart}.SEIZURE.label, 2))], 3, 1)';
        %     cfg_raster.cmapneurons  = cmap;
%         cfg_raster.cmapneurons  = gray(size(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).label, 2));
        ft_spike_plot_raster(cfg_raster, SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)));
        
        xlim(xlimits);
        xticks([xlimits(1), 0, xlimits(2)])
        y = ylim;
        ylim([1, ceil(max(abs(y)))]);
        yticks([1, ceil(max(abs(y)))]);
        set(gca,'ticklength', [0.005 0.100]);

        if itemplate == 1
            l1 = ylabel('Trial');
            set(l1, 'Units', 'normalized');
            lpos = l1.Position;
            set(l1, 'Position', [-0.1, 0.5, 0]);
        else
            ylabel([]);
        end

    end
    
end


% write to figure
fname = fullfile(config{1}.imagesavedir, 'article', 'rasterplots_templates');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));
