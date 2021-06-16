function [dat, stats_all, stats_responsive] = stats_unit_density(cfg, force)

set(groot,'defaultAxesTickLabelInterpreter','none');
hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];
fname_out   = fullfile(cfg{1}.datasavedir, 'stats_unit_density.mat');

%% load data
if force == false && exist(fname_out, 'file')
    load(fname_out, 'dat', 'stats_all', 'stats_responsive');
else
    dat = [];
    dat.trialinfo = table;
    dat.time = [];
    dat.label{1} = 'all';
    dat.dimord = 'rpt_time';
    iunit = 1; 
    for ipatient = 1 : 7
        %
        cfg{ipatient}.postfix   = [];
        SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg{ipatient});
        %         cfg{ipatient}.postfix   = '-windowed';
        %         SpikeStats_windowed     = spikeTrialStats(cfg{ipatient});
        SpikeDensity_timelocked = spikeTrialDensity(cfg{ipatient});
        
        for ipart = 1 : size(SpikeDensity_timelocked, 2)
            
            % if no good units were found
            if isempty(SpikeDensity_timelocked{ipart})
                continue
            end
            
            for markername = string(fields(SpikeDensity_timelocked{ipart}.stat))'
                
                % loop over units for spikedensity per unit/trial
                for ilabel = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername), 2)
                    
                    for hyplabel = ["all", hyplabels]
                        
                        fprintf('Adding: patient %d, part %d, unit %d of %d, stage: %s\n', ipatient, ipart, ilabel, size(SpikeDensity_timelocked{ipart}.stat.(markername), 2), hyplabel);
                        
                        if strcmp(hyplabel, "all")
                            trials = 1 : size(SpikeTrials_timelocked{ipart}.(markername).trialinfo, 1);
                        else
                            trials = SpikeTrials_timelocked{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                        end
                        temp                = table;
                        
                        if any(trials)
                            dat.trial{iunit} = squeeze(nanmean(SpikeDensity_timelocked{ipart}.sdf_lin.(markername).trial(trials, ilabel, :), 1))';
                            dat.time{iunit}  = SpikeDensity_timelocked{ipart}.sdf_lin.(markername).time;                
                        else
                            sprintf('No trials found!\n');
                        end
                        temp.cluster_group  = string(deblank(SpikeTrials_timelocked{ipart}.(markername).cluster_group{ilabel}));
                        temp.ipatient       = repmat(ipatient, size(temp, 1), 1);
                        temp.ipart          = repmat(ipart, size(temp, 1), 1);
                        temp.ilabel         = repmat(ilabel, size(temp, 1), 1);
                        temp.markername     = repmat(markername, size(temp, 1), 1);
                        temp.hyplabel       = repmat(hyplabel, size(temp, 1), 1);
                        temp.label          = repmat(string(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.label), size(temp, 1), 1);
                        
                        % check if unit responds statistically
                        responsive_pos  = false;
                        responsive_neg  = false;
                        if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}, 'posclusters')
                            for ipos = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.posclusters, 2)
                                if SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.posclusters(ipos).prob < 0.01
                                    responsive_pos  = true;
                                end
                            end
                        end
                        if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}, 'negclusters')
                            for ineg = 1 : size(SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.negclusters, 2)
                                if SpikeDensity_timelocked{ipart}.stat.(markername){ilabel}.negclusters(ineg).prob < 0.01
                                    responsive_neg  = true;
                                end
                            end
                        end
                        temp.responsive     = repmat(responsive_pos | responsive_neg, size(temp, 1), 1);
                        temp.responsive_pos = repmat(responsive_pos, size(temp, 1), 1);
                        temp.responsive_neg = repmat(responsive_neg, size(temp, 1), 1);
                        
                        % concatinate
                        dat.trialinfo = vertcat(dat.trialinfo, temp);
                        iunit = iunit + 1;
                    end
                end
            end
        end
    end
    

%% prepare data for stats and plotting

% stats on timecourses
cfgtemp                                         = [];
cfgtemp.trials                                  = dat.trialinfo.hyplabel ~= "all";
dat_all                                         = ft_selectdata(cfgtemp, dat);

cfgtemp.trials                                  = dat.trialinfo.hyplabel ~= "all" & dat.trialinfo.responsive;
dat_responsive                                  = ft_selectdata(cfgtemp, dat);

cfgtemp                                         = [];
cfgtemp.statistic                               = 'ft_statfun_correlationT';
cfgtemp.type                                    = 'Spearman';
cfgtemp.alpha                                   = 0.05;
cfgtemp.clusteralpha                            = 0.05;
cfgtemp.method                                  = 'montecarlo';
cfgtemp.computestat                             = 'yes';
cfgtemp.correctm                                = 'cluster';
cfgtemp.tail                                    = 0;
cfgtemp.correcttail                             = 'prob';
cfgtemp.numrandomization                        = 1000;
% cfgtemp.latency                                 = [cfg.stats.bl.(markername)(2) sdf_bar{ipart}.time(end)]; % active perio starts after baseline
cfgtemp.ivar                                    = 1;

cfgtemp.design                                  = nan(1, size(dat_all.trial, 2));
for i = 1 : size(hyplabels, 2)
    cfgtemp.design(dat_all.trialinfo.hyplabel == hyplabels{i}) = i;
end
stats_all                                       = ft_timelockstatistics(cfgtemp, dat_all);

cfgtemp.design                                  = nan(1, size(dat_responsive.trial, 2));
for i = 1 : size(hyplabels, 2)
    cfgtemp.design(dat_responsive.trialinfo.hyplabel == hyplabels{i}) = i;
end
stats_responsive                                = ft_timelockstatistics(cfgtemp, dat_responsive);

% save data
save(fname_out, 'dat', 'stats_all', 'stats_responsive', '-v7.3');

end

%% plotting

% baseline correction
dat_bl = dat;
for i = 1 : size(dat.trial, 2)
   dat_bl.trial{i} = dat_bl.trial{i} - nanmedian(dat.trial{i});
end


fig = figure;
set(gcf, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0.2 0.2 0.8 0.8]);
set(fig, 'Renderer', 'Painters');

selindx = dat.trialinfo.responsive & dat.trialinfo.cluster_group == "good";
% selindx = ones(size(dat.trialinfo, 1), 1);

cm  = cool(5);
xl = cfg{1}.epoch.toi.Hspike;
    
% firing rate all
subplot(2, 3, 1); hold;
indx = find(dat.trialinfo.hyplabel == "all" & selindx);
for i = indx'
    lh = plot(dat.time{i}, dat.trial{i}); lh.Color  = [0, 0, 0, 0.3];
end
x = dat.time{1}; y = nanmean(cat(1, dat.trial{indx}), 1); s = nanstd(cat(1, dat.trial{indx}), 1);
i = ~isnan(y); x = x(i); y = y(i); s = s(i);
fill([x, x(end:-1:1)], [y-s, y(end:-1:1)+s(end:-1:1)], [1, 1, 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
plot(dat.time{1}, nanmean(cat(1, dat.trial{indx}), 1), 'k', 'linewidth', 3);
plot(dat.time{1}, nanmean(cat(1, dat.trial{indx}), 1), 'w', 'linewidth', 2);
ylabel('Firing Rate'); axis tight; xlim(xl); ylim([0, 30]);

% firing rate per sleepstage
subplot(2, 3, 4); hold;
ic = 1;
for hyplabel = hyplabels
    indx = find(dat.trialinfo.hyplabel == hyplabel & selindx);
    for i = indx'
        lh = plot(dat.time{i}, dat.trial{i}); lh.Color  = [cm(ic, :), 0.3];
    end
    ic = ic + 1;
end
ic = 1;
for hyplabel = hyplabels
    indx = find(dat.trialinfo.hyplabel == hyplabel & selindx);    
    x = dat.time{1}; y = nanmean(cat(1, dat.trial{indx}), 1); s = nanstd(cat(1, dat.trial{indx}), 1);
    i = ~isnan(y); x = x(i); y = y(i); s = s(i);
    fill([x, x(end:-1:1)], [y-s, y(end:-1:1)+s(end:-1:1)], cm(ic, :), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);    
    ic = ic + 1;
end
ic = 1;
for hyplabel = hyplabels
    indx = find(dat.trialinfo.hyplabel == hyplabel & selindx);
    plot(dat.time{1}, nanmean(cat(1, dat.trial{indx}), 1), 'k', 'linewidth', 3);
    plot(dat.time{1}, nanmean(cat(1, dat.trial{indx}), 1), 'color', cm(ic, :), 'linewidth', 2);
    ic = ic + 1;    
end
ylabel('Relative Firing Rate'); axis tight; xlim(xl); ylim([0, 20]);

% FR Baseline corrected
subplot(2, 3, 2); hold;
indx = find(dat_bl.trialinfo.hyplabel == "all" & selindx);
for i = indx'
    lh = plot(dat_bl.time{i}, dat_bl.trial{i}); lh.Color  = [0, 0, 0, 0.3];
end
x = dat_bl.time{1}; y = nanmean(cat(1, dat_bl.trial{indx}), 1); s = nanstd(cat(1, dat_bl.trial{indx}), 1);
i = ~isnan(y); x = x(i); y = y(i); s = s(i);
fill([x, x(end:-1:1)], [y-s, y(end:-1:1)+s(end:-1:1)], [1, 1, 1], 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
plot(dat_bl.time{1}, nanmean(cat(1, dat_bl.trial{indx}), 1), 'k', 'linewidth', 3);
plot(dat_bl.time{1}, nanmean(cat(1, dat_bl.trial{indx}), 1), 'w', 'linewidth', 2);
ylabel('Firing Rate'); axis tight; xlim(xl); ylim([-20, 20]);

% firing rate per sleepstage
subplot(2, 3, 5); hold;
ic = 1;
for hyplabel = hyplabels
    indx = find(dat_bl.trialinfo.hyplabel == hyplabel & selindx);
    for i = indx'
        lh = plot(dat_bl.time{i}, dat_bl.trial{i}); lh.Color  = [cm(ic, :), 0.3];
    end
    ic = ic + 1;
end
ic = 1;
for hyplabel = hyplabels
    indx = find(dat_bl.trialinfo.hyplabel == hyplabel & selindx);    
    x = dat_bl.time{1}; y = nanmean(cat(1, dat_bl.trial{indx}), 1); s = nanstd(cat(1, dat_bl.trial{indx}), 1);
    i = ~isnan(y); x = x(i); y = y(i); s = s(i);
    fill([x, x(end:-1:1)], [y-s, y(end:-1:1)+s(end:-1:1)], cm(ic, :), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);    
    ic = ic + 1;
end
ic = 1;
for hyplabel = hyplabels
    indx = find(dat_bl.trialinfo.hyplabel == hyplabel & selindx);
    plot(dat_bl.time{1}, nanmean(cat(1, dat_bl.trial{indx}), 1), 'k', 'linewidth', 3);
    plot(dat_bl.time{1}, nanmean(cat(1, dat_bl.trial{indx}), 1), 'color', cm(ic, :), 'linewidth', 2);
    ic = ic + 1;    
end
ylabel('Relative Firing Rate'); axis tight; xlim(xl); ylim([-10, 10]);

% stats all
subplot(2, 3, 3); hold;
plot(stats_all.time, stats_all.rho, 'k');
ylabel('Correlation {\it r}');
yyaxis right
set(gca,'YColor','k');
lh = plot(stats_all.time, stats_all.stat);
lh.Color = [0,0,0,0]; % turn fully transparant

ylabel('t');
axis tight
ax = axis;
for ineg = 1 : size(stats_all.negclusters)
    if stats_all.negclusters(ineg).prob < 0.05
        x(1) = stats_all.time(find(stats_all.negclusterslabelmat == ineg, 1, 'first'));
        x(2) = stats_all.time(find(stats_all.negclusterslabelmat == ineg, 1, 'last'));
        height = ax(4) - (ax(4) - ax(3))*0.05;
        text(x(1)+0.1, height, sprintf('p=%.3f', stats_all.negclusters(ineg).prob), 'HorizontalAlignment', 'left', 'fontsize', 8);
        patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [1 0 0], 'facealpha', 0.3, 'edgecolor', 'none');
    end
end
for ipos = 1 : size(stats_all.posclusters)
    if stats_all.posclusters(ipos).prob < 0.05
        x(1) = stats_all.time(find(stats_all.posclusterslabelmat == ipos, 1, 'first'));
        x(2) = stats_all.time(find(stats_all.posclusterslabelmat == ipos, 1, 'last'));
        height = ax(4) - (ax(4) - ax(3))*0.05;
        text(x(1)+0.1, height, sprintf('p=%.3f', stats_all.posclusters(ipos).prob), 'HorizontalAlignment', 'left', 'fontsize', 8);
        patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [0 1 0], 'facealpha', 0.3, 'edgecolor', 'none');
    end
end
plot([ax(1), ax(2)],[0 0],'k:');
set(gca,'children',flipud(get(gca,'children')));
xlim(xl);


% stats responsive
subplot(2, 3, 6); hold;
plot(stats_responsive.time, stats_responsive.rho, 'k');
ylabel('Correlation {\it r}');
yyaxis right
set(gca,'YColor','k');
lh = plot(stats_responsive.time, stats_responsive.stat);
lh.Color = [0,0,0,0]; % turn fully transparant

ylabel('t');
axis tight
ax = axis;
for ineg = 1 : size(stats_responsive.negclusters)
    if stats_responsive.negclusters(ineg).prob < 0.05
        x(1) = stats_responsive.time(find(stats_responsive.negclusterslabelmat == ineg, 1, 'first'));
        x(2) = stats_responsive.time(find(stats_responsive.negclusterslabelmat == ineg, 1, 'last'));
        height = ax(4) - (ax(4) - ax(3))*0.05;
        text(x(1)+0.1, height, sprintf('p=%.3f', stats_responsive.negclusters(ineg).prob), 'HorizontalAlignment', 'left', 'fontsize', 8);
        patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [1 0 0], 'facealpha', 0.3, 'edgecolor', 'none');
    end
end
for ipos = 1 : size(stats_responsive.posclusters)
    if stats_responsive.posclusters(ipos).prob < 0.05
        x(1) = stats_responsive.time(find(stats_responsive.posclusterslabelmat == ipos, 1, 'first'));
        x(2) = stats_responsive.time(find(stats_responsive.posclusterslabelmat == ipos, 1, 'last'));
        height = ax(4) - (ax(4) - ax(3))*0.05;
        text(x(1)+0.1, height, sprintf('p=%.3f', stats_responsive.posclusters(ipos).prob), 'HorizontalAlignment', 'left', 'fontsize', 8);
        patch([x(1), x(2), x(2), x(1)], [ax(3), ax(3), ax(4), ax(4)], [0 1 0], 'facealpha', 0.3, 'edgecolor', 'none');
    end
end
plot([ax(1), ax(2)],[0 0],'k:');
set(gca,'children',flipud(get(gca,'children')));
xlim(xl);

% print to file
fname = fullfile(cfg{1}.imagesavedir, 'stats', 'stats_unit_density_SHARP');
disp(['Exporting figure ', fname])
exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff'], 'resolution', 600);
close all

disp('Done');

