function plotGrandAverageTimelocked(cfg, LFP, dat, dat_hyp, trialinfo)

hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

% baseline correction
clear dat_bl
baseline_all = [-0.5 -0.2];
baseline_hyp = [dat{1}.time(1) dat{1}.time(end)]; % extract from dat - max 


for iunit = 1 : size(dat, 2)
    
    fprintf('Processing %d of %d', iunit, size(dat, 2));
    
    cfgtemp                 = [];
    cfgtemp.latency         = baseline_all;
    sel                     = ft_selectdata(cfgtemp, dat{iunit});
    bl                      = nanmean(sel.avg);
    dat_bl{iunit}           = dat{iunit};
    dat_bl{iunit}.avg       = (dat{iunit}.avg ./ bl) * 100;
    dat_bl_norm{iunit}      = dat{iunit};
    dat_bl_norm{iunit}.avg  = ((dat{iunit}.avg - bl) ./ (dat{iunit}.avg + bl)) * 100;
    
    for hyplabel = hyplabels
        cfgtemp                                 = [];
        cfgtemp.latency                         = baseline_hyp;
        sel                                     = ft_selectdata(cfgtemp, dat_hyp.(hyplabel){iunit});
        bl                                      = nanmean(sel.avg);
        dat_hyp_bl.(hyplabel){iunit}            = dat_hyp.(hyplabel){iunit};
        dat_hyp_bl.(hyplabel){iunit}.avg        = (dat_hyp_bl.(hyplabel){iunit}.avg ./ bl) * 100;
        dat_hyp_bl_norm.(hyplabel){iunit}       = dat_hyp.(hyplabel){iunit};
        dat_hyp_bl_norm.(hyplabel){iunit}.avg   = ((dat_hyp.(hyplabel){iunit}.avg + bl) ./ (dat_hyp.(hyplabel){iunit}.avg - bl)) * 100;
    end
    
end

% combine all with units as repetitions
cfgtemp                     = [];
cfgtemp.keepindividual      = 'yes';
GA                          = ft_timelockgrandaverage(cfgtemp, dat{:}); 
GA_bl                       = ft_timelockgrandaverage(cfgtemp, dat_bl{:}); 
GA_bl_norm                  = ft_timelockgrandaverage(cfgtemp, dat_bl_norm{:}); 

clear dat dat_bl dat_nl_norm

for hyplabel = hyplabels
    GA_hyp.(hyplabel)           = ft_timelockgrandaverage(cfgtemp, dat_hyp.(hyplabel){:});
    GA_hyp_bl.(hyplabel)        = ft_timelockgrandaverage(cfgtemp, dat_hyp_bl.(hyplabel){:});
    GA_hyp_bl_norm.(hyplabel)   = ft_timelockgrandaverage(cfgtemp, dat_hyp_bl_norm.(hyplabel){:});
end

cfgtemp                                         = [];
cfgtemp.statistic                               = 'ft_statfun_correlationT';
cfgtemp.alpha                                   = 0.05;
cfgtemp.clusteralpha                            = 0.05;
cfgtemp.method                                  = 'montecarlo';
cfgtemp.computestat                             = 'yes';
cfgtemp.correctm                                = 'cluster';
cfgtemp.parameter                               = 'individual';
% cfgtemp.latency                                 = [cfg.stats.bl.(markername)(2) sdf_bar{ipart}.time(end)]; % active perio starts after baseline
cfgtemp.ivar                                    = 1;
% cfgtemp.uvar                                    = 2;
cfgtemp.design(1, :)                            = [ones(1, size(GA_bl.individual, 1)) ones(1, size(GA_bl.individual, 1))*2 ones(1, size(GA_bl.individual, 1))*3 ones(1, size(GA_bl.individual, 1))*4 ones(1, size(GA_bl.individual, 1))*5];
% cfgtemp.design(2, :)                            = [1 :  size(GA_bl.individual, 1), 1 :  size(GA_bl.individual, 1), 1 :  size(GA_bl.individual, 1), 1 :  size(GA_bl.individual, 1), 1 :  size(GA_bl.individual, 1) ];
cfgtemp.numrandomization                        = 1000;
stats                                           = ft_timelockstatistics(cfgtemp, GA_hyp.AWAKE, GA_hyp.REM, GA_hyp.PHASE_1, GA_hyp.PHASE_2, GA_hyp.PHASE_3);

% create average LFPs
for ipatient = 1 : size(LFP, 2)
    for ipart = 1 : 3
        for markername = string(unique(cfg{ipatient}.editmarkerfile.torename(:, 2))')
            LFP_avg{ipatient}{ipart}.(markername) = ft_timelockanalysis([], LFP{ipatient}{ipart}.(markername));
        end
    end
end

%% PLOTTING

xl = [-0.2 0.6];
nrows = 4;
fig = figure;
set(gcf, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);

channel = 1;
cm  = cool(5);

% FIRST ROW: LFP 
for icol = 1 : 6
    subplot(nrows, 6, icol);
    title('LFP');
    hold;
    itemp = 1;
    clear y
    for ipatient = 1 : size(LFP, 2)
        for ipart = 1 : 3
            for markername = string(unique(cfg{ipatient}.editmarkerfile.torename(:, 2))')
                x = LFP_avg{ipatient}{ipart}.(markername).time;
                y(itemp, :) = LFP_avg{ipatient}{ipart}.(markername).avg(channel, :);
                lh = plot(x, y(itemp, :), 'k');
                lh.Color        = [0,0,0,0.3];
                itemp = itemp + 1;
            end
        end
    end
    plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
    plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
    ylabel('Current');
    axis tight
    xlim(xl);
end

% SECOND ROW: FIRING RATE
subplot(nrows, 6, 7);
title('Firing rate');
hold on
clear y
iunit = 1;
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "good ")'
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl.individual(indx, 1, :));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1;
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate vs. baseline');
axis tight
xlim(xl);

subplot(nrows, 6, 8);
hold on
clear y
iunit = 1;
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "good ")'
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl_norm.individual(indx, 1, :));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1;    
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate change (norm. %)');
ylim([-100, 100]);
axis tight
xlim(xl);

subplot(nrows, 6, 9);
hold on
clear y
iunit = 1;
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "good ")'
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl.individual(indx, 1, :));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1;       
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate (log(Hz))');
axis tight
set(gca, 'YScale', 'log');

subplot(nrows, 6, 10);
hold on
clear y
iunit = 1;
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")'
    x               = GA.time;
    y(iunit,:)      = (squeeze(GA_bl.individual(indx, 1, :)));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1;
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate (Hz)');
axis tight

subplot(nrows, 6, 11);
hold on
clear y
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")'
    x               = GA.time;
    y(iunit,:)      = (squeeze(GA_bl_norm.individual(indx, 1, :)));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1; 
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate change (norm. %)');
ylim([-100, 100]);
axis tight

subplot(nrows, 6, 12);
hold on
clear y
iunit = 1;
for indx = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")'
    x               = GA.time;
    y(iunit,:)      = (squeeze(GA_bl.individual(indx, 1, :)));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
    iunit           = iunit + 1;    
end
plot(x, nanmean(y, 1), 'k', 'linewidth', 3);
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate (log(Hz))');
axis tight
set(gca, 'YScale', 'log');

% THIRD ROW: FIRING RATE PER HYPNOGRAM
subplot(nrows, 6, 13);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "good ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
ylabel('Firing rate');
xlim(xl);


subplot(nrows, 6, 14);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "good ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp_bl.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
axis tight
ylabel('Firing rate');
xlim(xl);


subplot(nrows, 6, 15);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "good ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp_bl.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
axis tight
ylabel('Firing rate');
xlim(xl);
set(gca, 'YScale', 'log');


subplot(nrows, 6, 16);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
ylabel('Firing rate');
xlim(xl);


subplot(nrows, 6, 17);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp_bl.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
axis tight
ylabel('Firing rate');
xlim(xl);


subplot(nrows, 6, 18);
hold on
clear y
ihyp = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    iunit = 1;
    indxs = find(trialinfo.responsive & trialinfo.cluster_group == "mua  ")';
    for indx = indxs
        lh = plot(x, squeeze(GA_hyp_bl.(hyplabel).individual(indx, 1, :)), 'k');
        lh.Color = [cm(ihyp, :), 0.1];
        iunit = iunit + 1;
    end
    ihyp = ihyp + 1;
end

i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual(indxs, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
axis tight
ylabel('Firing rate');
xlim(xl);
set(gca, 'YScale', 'log');

for iplot = 0 : 5
    
    % Forth row: statistics
    subplot(nrows, 6, 19 + iplot); hold;
    plot(stats.time, stats.rho, 'k');
    ylabel('Correlation {\it r}');
    yyaxis right
    set(gca,'YColor','k');
    lh = plot(stats.time, stats.stat);
    lh.Color        = [0,0,0,0]; % turn fully transparant
    
    ylabel('t');
    axis tight
    ax = axis;
    for ineg = 1 : size(stats.negclusters)
        if stats.negclusters(ineg).prob < 0.025
            x(1) = stats.time(find(stats.negclusterslabelmat == ineg, 1, 'first'));
            x(2) = stats.time(find(stats.negclusterslabelmat == ineg, 1, 'last'));
            height = ax(4) - (ax(4) - ax(3))*0.05;
            text(x(1)+0.1, height, sprintf('p = %.3f', stats.negclusters(ineg).prob), 'HorizontalAlignment', 'left', 'fontsize', 8);
            patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [1 0 0], 'facealpha', 0.1, 'edgecolor', 'none');
        end
    end
    for ipos = 1 : size(stats.posclusters)
        if stats.posclusters(ipos).prob < 0.025
            x(1) = stats.time(find(stats.posclusterslabelmat == ipos, 1, 'first'));
            x(2) = stats.time(find(stats.posclusterslabelmat == ipos, 1, 'last'));
            height = ax(4) - (ax(4) - ax(3))*0.05;
            text(x(1)+0.1, height, sprintf('p = %.3f', stats.posclusters(ipos).prob), 'HorizontalAlignment', 'left', 'fontsize', 6);
            patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [0 1 0], 'facealpha', 0.1, 'edgecolor', 'none');
        end
    end
    plot([ax(1), ax(2)],[0 0],'k:');
    set(gca,'children',flipud(get(gca,'children')));
    xlim(xl);
    
end

% print ISI to file
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
fname = fullfile(cfg{1}.imagesavedir, 'GA');
print(fig, '-dpdf', fname);
print(fig, '-djpeg', fname, '-r1200');

disp(['Export fig ', fname])
export_fig(fname, '-png'); % need to install https://www.ghostscript.com/download/gsdnld.html
close all



