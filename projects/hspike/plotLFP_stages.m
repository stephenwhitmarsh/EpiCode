function LFP_avg = plotLFP_stages(cfg, LFP, marker, hypnogram, force)

cfg.visible = ft_getopt(cfg, 'visible', 'on');

fname = fullfile(cfg.datasavedir, sprintf('%sLFP_stages.mat', cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('**********************************\n');
    fprintf('** loading %s  stats *****\n', fname);
    fprintf('**********************************\n\n');
    load(fname, 'LFP_avg');
    return
end

for markername = string(cfg.LFP.name)
    
    for istage = 0 : 4
        
        hasdata = [];
        for ipart = 1 : 3
            if ipart > numel(LFP)
                continue
            end
            if ~isfield(LFP{ipart}, markername)
                continue
            end
            if isempty(LFP{ipart}.(markername))
                continue
            end
            
            cfgtemp = [];
            cfgtemp.trials = find(LFP{ipart}.(markername).trialinfo(:,4) == istage);
            if ~isempty(cfgtemp.trials)
                LFP_sel{ipart} = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                hasdata = [hasdata, ipart];
            end
        end
        if isempty(hasdata)
            continue
        end
        LFP_app = ft_appenddata([], LFP_sel{hasdata});
        clear LFP_sel
        
        cfgtemp = [];
        cfgtemp.baseline = 'yes';
        cfgtemp.baseline = cfg.LFP.baselinewindow.(markername);
        LFP_avg.(markername){istage+1} = ft_timelockbaseline(cfgtemp, LFP_app);
        LFP_avg.(markername){istage+1} = ft_timelockanalysis([], LFP_avg.(markername){istage+1});
        LFP_avg.(markername){istage+1}.nr = size(LFP_app.trial, 2);
        
    end
end

clear LFP


% fig = figure('visible', cfg.visible);
fig = figure;
fig.Renderer = 'Painters';
cm = cool(5);

% max height for scaling
h = 0;

for markername = string(fieldnames(LFP_avg))'
    if isempty(LFP_avg.(markername){1})
        continue
    end
    for ichan = 1 : size(LFP_avg.(markername){1}.label,1)
        for istage = [4, 0, 1, 2, 3]
            try
                y = LFP_avg.(markername){istage+1}.avg(ichan,:);
                ystd = sqrt(LFP_avg.(markername){istage+1}.var(ichan,:));
                h = max(abs([h, y+ystd, y-ystd]));
            catch
            end
        end
    end
end

% plot separate LFPs for templates and sleepstage
nrtemplates = size(fieldnames(LFP_avg), 1);
iplot = 1;

for markername = string(fieldnames(LFP_avg))'
    
    clear totaldur totalsum
    for i = 0 : 4
        totaldur(i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i)));
    end
    for i = 0 : 4
        totalsum(i+1) = sum(marker.stage == i & strcmp(string(marker.name), markername));
    end
    IEDrate = totalsum ./ totaldur * 60;
    
    cm  = cool(5);
    x   = [2 3 4 5 1];
    
    subplot(2, nrtemplates, iplot); hold;
    title(sprintf('%s', markername));

    % IED rate normalized to wake
    IEDrateNorm = IEDrate ./ IEDrate(1);
    for ib = 1:5
        hb = bar(x(ib), IEDrateNorm(ib), 1);
        xData = hb.XData + hb.XOffset;
        text(xData,IEDrateNorm(ib), num2str(round(IEDrateNorm(ib),1)),'vert','bottom','horiz','center');
        set(hb, 'FaceColor', cm(x(ib),:));
    end
    xlim([0,6]);
    ylabel('IED rate vs. wake');
    box off
    if max(IEDrateNorm) > 0
        ylim([0, max(IEDrateNorm) * 1.3]);
    end
    set(gca,'TickLength',[0 0])
    set(gca,'Xticklabels', []);    
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';   
    
   % legend  
    l = {'W n=0','S1 n=0','S2 n=0','S3 n=0',' n=0','W n=0'};
    for istage = 1 : size(LFP_avg.(markername),2)
        switch istage
            case 1
                if ~isempty(LFP_avg.(markername){istage})
                    l{2} = sprintf('W n=%d', totalsum(istage));
                else
                    l{2} = 'W n=0';
                end            
            case 2
                if ~isempty(LFP_avg.(markername){istage})
                    l{3} = sprintf('S1 n=%d',totalsum(istage));
                else
                    l{3} = 'S1 n=0';
                end
            case 3
                if ~isempty(LFP_avg.(markername){istage})
                    l{4} = sprintf('S2 n=%d', totalsum(istage));
                else
                    l{4} = 'S2 n=0';
                end
            case 4
                if ~isempty(LFP_avg.(markername){istage})
                    l{5} = sprintf('S3 n=%d', totalsum(istage));
                else
                    l{5} = 'S3 n=0';
                end
            case 5
                if ~isempty(LFP_avg.(markername){istage})
                    l{1} = sprintf('REM n=%d', totalsum(istage));
                else
                    l{1} = 'REM n=0';
                end
        end
    end
    
    clear p
    for ii = 1:size(cm,1)
        p(ii) = patch(NaN, NaN, cm(ii,:));
    end
    legend(p, l, 'location','northoutside');
    
    
    % plot LFP
    subplot(2, nrtemplates, nrtemplates + iplot); hold;
    title(sprintf('%s', markername));
    
    n = 1;
    ytick = [];
    
    for ichan = 1 : size(LFP_avg.(markername){1}.label,1)
        ytick = [ytick, n*h];
        icolor = 1;
        for istage = [4, 0, 1, 2, 3]
            try
                x = LFP_avg.(markername){istage+1}.time;
                y = LFP_avg.(markername){istage+1}.avg(ichan,:);
                ystd = sqrt(LFP_avg.(markername){istage+1}.var(ichan,:));
                patch([x x(end:-1:1)], [y-ystd y(end:-1:1)+ystd] + n*h, cm(icolor,:), 'edgecolor', 'none', 'facealpha', 0.2);
            catch
            end
            icolor = icolor + 1;
        end
        icolor = 1;
        
        for istage = [4, 0, 1, 2, 3]
            try
                x = LFP_avg.(markername){istage+1}.time;
                y = LFP_avg.(markername){istage+1}.avg(ichan,:);
                ystd = sqrt(LFP_avg.(markername){istage+1}.var(ichan,:));
                plot(x,y + n*h, 'color', cm(icolor,:),'linewidth', 1);
            catch
            end
            icolor = icolor + 1;
        end
        temp = strsplit(LFP_avg.(markername){1}.label{ichan}, '-');
        label{ichan} = temp{1};
        n = n + 1;
    end
    yticks(ytick);
    set(gca,'TickLabelInterpreter', 'none')
%     set(gca,'fontsize', 8)
    yticklabels(label);
    xlabel('Time (s)');
    axis tight;
    xlim([LFP_avg.(markername){1}.time(1), LFP_avg.(markername){1}.time(end)]);
    
    % legend  
    l = {'W n=0','S1 n=0','S2 n=0','S3 n=0',' n=0','W n=0'};
    for istage = 1 : size(LFP_avg.(markername),2)
        switch istage
            case 1
                if ~isempty(LFP_avg.(markername){istage})
                    l{2} = sprintf('W n=%d', LFP_avg.(markername){istage}.nr);
                else
                    l{2} = 'W n=0';
                end            
            case 2
                if ~isempty(LFP_avg.(markername){istage})
                    l{3} = sprintf('S1 n=%d', LFP_avg.(markername){istage}.nr);
                else
                    l{3} = 'S1 n=0';
                end
            case 3
                if ~isempty(LFP_avg.(markername){istage})
                    l{4} = sprintf('S2 n=%d', LFP_avg.(markername){istage}.nr);
                else
                    l{4} = 'S2 n=0';
                end
            case 4
                if ~isempty(LFP_avg.(markername){istage})
                    l{5} = sprintf('S3 n=%d', LFP_avg.(markername){istage}.nr);
                else
                    l{5} = 'S3 n=0';
                end
            case 5
                if ~isempty(LFP_avg.(markername){istage})
                    l{1} = sprintf('REM n=%d', LFP_avg.(markername){istage}.nr);
                else
                    l{1} = 'REM n=0';
                end

        end
    end
    
    clear p
    for ii = 1:size(cm,1)
        p(ii) = patch(NaN, NaN, cm(ii,:));
    end
    legend(p, l, 'location','northoutside');
    
    iplot = iplot + 1;
end

set(fig,'PaperOrientation', 'landscape');
set(fig,'PaperUnits', 'normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'LFPxSleep.pdf')), '-r600');
print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'LFPxSleep.png')), '-r600');
close all

% save results
save(fname, 'LFP_avg');