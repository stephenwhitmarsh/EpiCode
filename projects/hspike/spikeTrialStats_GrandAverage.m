function [GA, stats] = spikeTrialStats_GrandAverage(cfg, SpikeTrials_timelocked, SpikeDensity_timelocked)

hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];

% organize all data in single struct array

iunit = 1;
clear dat
trialinfo = table;

for ipatient = 1 : size(SpikeDensity_timelocked, 2)
    
    if isempty(SpikeDensity_timelocked{ipatient})
        continue
    end
    
    for ipart = 1 : size(SpikeDensity_timelocked{ipatient}, 2)
        if isempty(SpikeDensity_timelocked{ipatient}{ipart})
            continue
        end
        
        for markername = string(fields(SpikeDensity_timelocked{ipatient}{ipart}.sdf_bar))'
            
            for ilabel = 1 : size(SpikeDensity_timelocked{ipatient}{ipart}.sdf_bar.combined1.label, 2)
                
                fprintf('Adding: patient %d, part %d, %s, %s\n', ipatient, ipart, markername, SpikeDensity_timelocked{ipatient}{ipart}.sdf_bar.combined1.label{ilabel});
                dat{iunit}                          = [];
                dat{iunit}.avg                      = SpikeDensity_timelocked{ipatient}{ipart}.sdf_lin.(markername).avg(ilabel, :);
                dat{iunit}.time                     = SpikeDensity_timelocked{ipatient}{ipart}.sdf_lin.(markername).time;
                dat{iunit}.label{1}                 = 'Cluster';
                dat{iunit}.dimord                   = 'chan_time';
                
                trialinfo.ipatient(iunit)           = ipatient;
                trialinfo.ipart(iunit)              = ipart;
                trialinfo.ilabel(iunit)             = ilabel;
                trialinfo.markername(iunit)         = markername;
                trialinfo.cluster_group{iunit}      = SpikeTrials_timelocked{ipatient}{ipart}.(markername).cluster_group{ilabel};
                if strfind(SpikeTrials_timelocked{ipatient}{ipart}.(markername).cluster_group{ilabel}, "good")
                    trialinfo.good(iunit)           = true;
                else
                    trialinfo.good(iunit)           = false;
                end

                trialinfo.responsive(iunit)         = false;    

                if isfield(SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}, 'posclusters')
                    for ipos = 1 : size(SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}.posclusters, 2)
                        if SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}.posclusters(ipos).prob < 0.01
                            trialinfo.responsive(iunit) = true;
                        end
                    end
                end
                if isfield(SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}, 'negclusters')
                    for ineg = 1 : size(SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}.negclusters, 2)                     
                        if SpikeDensity_timelocked{ipatient}{ipart}.stat.(markername){ilabel}.negclusters(ineg).prob < 0.01
                            trialinfo.responsive(iunit) = true;
                        end
                    end
                end
                
                SpikeTrials_timelocked{ipatient}{ipart}.(markername).trialinfo.hyplabel(SpikeTrials_timelocked{ipatient}{ipart}.(markername).trialinfo.hyplabel == "NO_SCORE") = "AWAKE";
                
                for hyplabel = hyplabels
                    
                    trials = SpikeTrials_timelocked{ipatient}{ipart}.(markername).trialinfo.hyplabel == hyplabel;
                    dat_hyp.(hyplabel){iunit}.avg                      = squeeze(nanmean(SpikeDensity_timelocked{ipatient}{ipart}.sdf_lin.(markername).trial(trials, ilabel, :), 1))';
                    dat_hyp.(hyplabel){iunit}.time                     = SpikeDensity_timelocked{ipatient}{ipart}.sdf_lin.(markername).time;
                    dat_hyp.(hyplabel){iunit}.label{1}                 = 'Cluster';
                    dat_hyp.(hyplabel){iunit}.dimord                   = 'chan_time';
                    dat_hyp.(hyplabel){iunit}.trialinfo.ipatient       = ipatient;
                    dat_hyp.(hyplabel){iunit}.trialinfo.ipart          = ipart;
                    dat_hyp.(hyplabel){iunit}.trialinfo.markername     = markername;
                    dat_hyp.(hyplabel){iunit}.trialinfo.cluster_group  = SpikeTrials_timelocked{ipatient}{ipart}.(markername).cluster_group{ilabel};
                    
                end
                
                iunit = iunit + 1;
            end
        end
    end
end


clear trialinfo_sum
trialinfo_sum = table;
i = 1;
for ipatient = unique(trialinfo.ipatient)'
    for ipart = 1 : 3
        trialinfo_sum.ipatient(i) = ipatient;
        trialinfo_sum.ipart(i)              = ipart;
        trialinfo_sum.units(i)              = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart);

        trialinfo_sum.responsive_p(i)       = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart)) * 100, 0);
 
        trialinfo_sum.sua(i)                = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true);
        trialinfo_sum.mua(i)                = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false);        

        %         trialinfo_sum.sua_responsive(i)     = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true & trialinfo.responsive == true);
%         trialinfo_sum.sua_unresponsive(i)   = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true & trialinfo.responsive == false);
        trialinfo_sum.sua_responsive_p(i)   = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == true)) * 100, 0);

%         trialinfo_sum.mua_responsive(i)     = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false & trialinfo.responsive == true);        
%         trialinfo_sum.mua_unresponsive(i)   = sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false & trialinfo.responsive == false);        
        trialinfo_sum.mua_responsive_p(i)   = round((sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false & trialinfo.responsive == true) / sum(trialinfo.ipatient == ipatient & trialinfo.ipart == ipart & trialinfo.good == false)) * 100, 0);

        i = i + 1;
    end
end

fname = fullfile(cfg{1}.datasavedir, 'overview_GA_all');
writetable(trialinfo, fname);
fname = fullfile(cfg{1}.datasavedir, 'overview_GA_sum');
writetable(trialinfo_sum, fname);


% baseline correction
clear dat_bl
for iunit = 1 : size(dat, 2)
    cfgtemp             = [];
%     cfgtemp.latency     = cfg{1}.stats.bl.(markername);
    sel                 = ft_selectdata(cfgtemp, dat{iunit});
    bl                  = nanmean(sel.avg);
    dat_bl{iunit}        = dat{iunit};
    dat_bl{iunit}.avg    = (dat_bl{iunit}.avg ./ bl);
    fprintf('Processing %d of %d', iunit, size(dat, 2));
    for hyplabel = hyplabels
        cfgtemp                          = [];
%         cfgtemp.latency                  = cfg{1}.stats.bl.(markername);
        sel                              = ft_selectdata(cfgtemp, dat_hyp.(hyplabel){iunit});
        bl                               = nanmean(sel.avg);
        dat_hyp_bl.(hyplabel){iunit}     = dat_hyp.(hyplabel){iunit};
%         dat_hyp_bl.(hyplabel){iunit}.avg = (dat_hyp_bl.(hyplabel){iunit}.avg ./ (dat_hyp_bl.(hyplabel){iunit}.avg + bl));
        dat_hyp_bl.(hyplabel){iunit}.avg = (dat_hyp_bl.(hyplabel){iunit}.avg ./ bl);
        
    end
end

% combine all with units as repetitions
cfgtemp                     = [];
cfgtemp.keepindividual      = 'yes';
GA                          = ft_timelockgrandaverage(cfgtemp, dat{trialinfo.responsive});
GA_bl                       = ft_timelockgrandaverage(cfgtemp, dat_bl{trialinfo.responsive});
for hyplabel = hyplabels
    GA_hyp.(hyplabel)       = ft_timelockgrandaverage(cfgtemp, dat_hyp.(hyplabel){trialinfo.responsive});    
    GA_hyp_bl.(hyplabel)    = ft_timelockgrandaverage(cfgtemp, dat_hyp_bl.(hyplabel){trialinfo.responsive});
end

% select only responsive units
trialinfo = trialinfo(trialinfo.responsive, :);

% plot individual units

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])
set(0,'defaultfigurecolor',[1 1 1]);

fig = figure;
set(gcf,'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);

%% FIRST COLUMN
subplot(3, 4, 1);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.1];
end
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);


subplot(3, 4, 5);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "mua")
        lh              = plot(x, y(iunit, :), 'k');
        lh.Color = [1,0,0, 0.1];
    end
end
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "good")
        lh              = plot(x, y(iunit, :), 'k');
        lh.Color = [0,1,0, 0.1];
    end
end
plot(x, nanmean(GA.individual(trialinfo.good == false, :), 1), 'k', 'linewidth', 3);
plot(x, nanmean(GA.individual(trialinfo.good == true, :), 1), 'k', 'linewidth', 3);
plot(x, nanmean(GA.individual(trialinfo.good == false, :), 1), 'r', 'linewidth', 2);
plot(x, nanmean(GA.individual(trialinfo.good == true, :), 1), 'g', 'linewidth', 2);
ylabel('Firing rate');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);


subplot(3, 4, 9);
cm  = cool(5);
hold on
clear y
i = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    for iunit = 1 : size(GA_bl.individual, 1)
        lh = plot(x, squeeze(GA_hyp.(hyplabel).individual(iunit, 1, :)), 'k');
        lh.Color = [cm(i, :), 0.05];
    end
    i = i + 1;    
end
i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(:, 1, :)), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp.(hyplabel).individual(:, 1, :)), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
ylabel('Firing rate');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);

%% SECOND COLUMN (LOG)

subplot(3, 4, 2);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    lh              = plot(x, log(y(iunit, :)), 'k');
    lh.Color        = [0,0,0,0.2];
end
plot(x, log(nanmean(y, 1)), 'w', 'linewidth', 2);
ylabel('Firing rate (log)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);

subplot(3, 4, 6);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "mua")
        lh              = plot(x, log(y(iunit, :)), 'k');
        lh.Color = [1,0,0, 0.2];
    end
end
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "good")
        lh              = plot(x, log(y(iunit, :)), 'k');
        lh.Color = [0,1,0, 0.2];
    end
end
plot(x, log(nanmean(GA.individual(trialinfo.good == false, :), 1)), 'k', 'linewidth', 3);
plot(x, log(nanmean(GA.individual(trialinfo.good == true, :), 1)), 'k', 'linewidth', 3);
plot(x, log(nanmean(GA.individual(trialinfo.good == false, :), 1)), 'r', 'linewidth', 2);
plot(x, log(nanmean(GA.individual(trialinfo.good == true, :), 1)), 'g', 'linewidth', 2);

ylabel('Firing rate (log)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);

subplot(3, 4, 10);
cm  = cool(5);
hold on
clear y
i = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    for iunit = 1 : size(GA_bl.individual, 1)
        lh = plot(x, log(squeeze(GA_hyp.(hyplabel).individual(iunit, 1, :))), 'k');
        lh.Color = [cm(i, :), 0.05];
    end
    i = i + 1;
end
i = 1;
for hyplabel = hyplabels
    plot(x, log(nanmean(squeeze(GA_hyp.(hyplabel).individual(:, 1, :)), 1)), 'k', 'linewidth', 3);
    plot(x, log(nanmean(squeeze(GA_hyp.(hyplabel).individual(:, 1, :)), 1)), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
ylabel('Firing rate (log)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);


%% THIRD COLUMN
subplot(3, 4, 3);
hold on
clear y
for iunit = 1 : size(GA_bl.individual, 1)
    x               = GA_bl.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    lh              = plot(x, y(iunit, :), 'k');
    lh.Color        = [0,0,0,0.2];
end
plot(x, nanmean(y, 1), 'w', 'linewidth', 2);
ylabel('Firing rate (relative)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);
ylim([0, 5]);

subplot(3, 4, 7);
hold on
clear y
for iunit = 1 : size(GA_bl.individual, 1)
    x               = GA_bl.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "mua")
        lh              = plot(x, y(iunit, :), 'k');
        lh.Color = [1,0,0, 0.2];
    end
end
clear y
for iunit = 1 : size(GA_bl.individual, 1)
    x               = GA_bl.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "good")
        lh              = plot(x, y(iunit, :), 'k');
        lh.Color = [0,1,0, 0.2];
    end
end
plot(x, nanmean(GA_bl.individual(trialinfo.good == false, :), 1), 'k', 'linewidth', 3);
plot(x, nanmean(GA_bl.individual(trialinfo.good == true, :), 1), 'k', 'linewidth', 3);
plot(x, nanmean(GA_bl.individual(trialinfo.good == false, :), 1), 'r', 'linewidth', 2);
plot(x, nanmean(GA_bl.individual(trialinfo.good == true, :), 1), 'g', 'linewidth', 2);
ylabel('Firing rate (relative)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);
ylim([0, 5]);

subplot(3, 4, 11);
cm  = cool(5);
hold on
clear y
i = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    for iunit = 1 : size(GA_bl.individual, 1)
        lh = plot(x, squeeze(GA_hyp_bl.(hyplabel).individual(iunit, 1, :)), 'k');
        lh.Color = [cm(i, :), 0.05];
    end
    i = i + 1;
end
i = 1;
for hyplabel = hyplabels
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual), 1), 'k', 'linewidth', 3);
    plot(x, nanmean(squeeze(GA_hyp_bl.(hyplabel).individual), 1), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);
ylabel('Firing rate (relative)');
ylim([0, 5]);

%% FOURTH COLUMN (LOG)

subplot(3, 4, 4);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    lh              = plot(x, log(y(iunit, :)), 'k');
    lh.Color        = [0,0,0,0.2];
end
plot(x, log(nanmean(y, 1)), 'w', 'linewidth', 2);
ylabel('Firing rate (relative log)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);

subplot(3, 4, 8);
hold on
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "mua")
        lh              = plot(x, log(y(iunit, :)), 'k');
        lh.Color = [1,0,0, 0.2];
    end
end
clear y
for iunit = 1 : size(GA.individual, 1)
    x               = GA.time;
    y(iunit,:)      = squeeze(GA_bl.individual(iunit, 1, :));
    if strfind(trialinfo.cluster_group{iunit}, "good")
        lh              = plot(x, log(y(iunit, :)), 'k');
        lh.Color = [0,1,0, 0.2];
    end
end
plot(x, log(nanmean(GA_bl.individual(trialinfo.good == false, :), 1)), 'k', 'linewidth', 3);
plot(x, log(nanmean(GA_bl.individual(trialinfo.good == true, :), 1)), 'k', 'linewidth', 3);
plot(x, log(nanmean(GA_bl.individual(trialinfo.good == false, :), 1)), 'r', 'linewidth', 2);
plot(x, log(nanmean(GA_bl.individual(trialinfo.good == true, :), 1)), 'g', 'linewidth', 2);

ylabel('Firing rate (relative log)');
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);

subplot(3, 4, 12);
cm  = cool(5);
hold on
clear y
i = 1;
for hyplabel = hyplabels
    x = GA_bl.time;
    for iunit = 1 : size(GA_bl.individual, 1)
        lh = plot(x, log(squeeze(GA_hyp_bl.(hyplabel).individual(iunit, 1, :))), 'k');
        lh.Color = [cm(i, :), 0.2];
    end
    lh = plot(x, log(nanmean((squeeze(GA_hyp_bl.(hyplabel).individual(:, 1, :))))), 'k', 'linewidth', 2);
    lh.Color = [cm(i, :), 1];
    i = i + 1;
end
i = 1;
for hyplabel = hyplabels
    plot(x, log(nanmean(squeeze(GA_hyp_bl.(hyplabel).individual), 1)), 'k', 'linewidth', 3);
    plot(x, log(nanmean(squeeze(GA_hyp_bl.(hyplabel).individual), 1)), 'color', cm(i, :), 'linewidth', 2);
    i = i + 1;
end
axis tight
xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);
ylabel('Firing rate (relative log)');

% set(findall(fig, '-property', 'fontsize'), 'fontsize', 8);
fname = fullfile(cfg{1}.imagesavedir, 'overview_GA');
export_fig(fname, '-png'); % need to install https://www.ghostscript.com/download/gsdnld.html
close all

%             
% 
%  
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % create dummy for baseline
% dummy = dat;
% for itemp = 1 : size(dummy, 2)
%     slim(1) = find(dat{itemp}.time > cfg{1}.stats.bl.(markername)(1), 1, 'first');
%     slim(2) = find(dat{itemp}.time < cfg{1}.stats.bl.(markername)(2), 1, 'last');
%     bl      = nanmean(dat{itemp}.avg(slim(1):slim(2)));
%     dummy{itemp}.avg = ones(size(dat{itemp}.avg)) .* bl;
% end
% 
% 
% 
%     % stats
%     cfgtemp                                         = [];
%     cfgtemp.parameter                               = 'avg';
%     cfgtemp.statistic                               = 'ft_statfun_depsamplesT';
%     cfgtemp.alpha                                   = cfg{1}.stats.alpha;
%     cfgtemp.clusteralpha                            = 0.01;
%     cfgtemp.method                                  = 'montecarlo';
%     cfgtemp.computestat                             = 'yes';
%     cfgtemp.correctm                                = 'cluster';
%     cfgtemp.latency                                 = [cfg{1}.stats.bl.combined1(2) dat{1}.time(end)];
%     cfgtemp.ivar                                    = 1;
%     cfgtemp.uvar                                    = 2;
%     cfgtemp.design(1, :)                            = [ones(1, size(dat, 2)) ones(1, size(dat, 2)) *2];
%     cfgtemp.design(2, :)                            = [1 : size(dat, 2) 1 : size(dat, 2)];
%     cfgtemp.numrandomization                        = 1000;
%     stat                                            = ft_timelockstatistics(cfgtemp, dat{:}, dummy{:});
%     
%     
%     cfgtemp                 = [];
%     cfgtemp.keepindividual  = 'yes';
%     GA_bl.(markername)      = ft_timelockgrandaverage(cfgtemp, dat_bl.(markername){:});
%     for hyplabel = hyplabels
%         GA_hyp_bl.(markername).(hyplabel)   = ft_timelockgrandaverage(cfgtemp, dat_hyp_bl.(markername).(hyplabel){:});
%     end
%     
%     
%     %% plotting
%     figure;
%     subplot(2, 6, 1); hold on
%     bar(GA_bl.(markername).time, squeeze(nanmean(GA_bl.(markername).individual(:, 1, :), 1))', 1, 'facecolor', [0 0 0], 'Edgecolor', 'none');
%     set(gca, 'box', 'off', 'XGrid', 'on', 'TickDir', 'out'); %, 'TickLength',[0 0]); %, 'Xticklabels', []);
%     ylabel('Relative firing rate');
%     
%     % plot positive clusters
%     lag = size(GA.(markername).time, 2) - size(stat.(markername).time, 2);
%     if isfield(stat.(markername), 'posclusters')
%         for ipos = 1 : size(stat.(markername).posclusters, 2)
%             if stat.(markername).posclusters(ipos).prob < cfg{1}.stats.alpha
%                 sel = find(stat.(markername).posclusterslabelmat == ipos);
%                 bar(stat.(markername).time(sel), squeeze(nanmean(GA_bl.(markername).individual(:, 1, sel+lag))), 1, 'facecolor', 'r', 'edgecolor', 'none');
%             end
%         end
%     end
%     
%     % plot negative clusters
%     if isfield(stat.(markername), 'negclusters')
%         for ineg = 1 : size(stat.(markername).negclusters, 2)
%             if stat.(markername).negclusters(ineg).prob < cfg{1}.stats.alpha
%                 sel = find(stat.(markername).negclusterslabelmat == ineg);
%                 bar(stat.(markername).time(sel), squeeze(nanmean(GA_bl.(markername).individual(:, 1, sel+lag))), 1, 'facecolor', 'b', 'edgecolor', 'none');
%             end
%         end
%     end
%     
%     % plot baseline patch
%     x = [cfg{1}.stats.bl.(markername)(1) cfg{1}.stats.bl.(markername)(2) cfg{1}.stats.bl.(markername)(2) cfg{1}.stats.bl.(markername)(1)];
%     axis tight
%     ax = axis;
%     y = [ax(3) ax(3) ax(4) ax(4)];
%     patch('XData', x, 'YData', y, 'facecolor', [0 0 0], 'edgecolor', 'none', 'facealpha', 0.1);
%     set(gca,'XGrid', 'on', 'box', 'off', 'TickDir', 'out');
%     
%     % plot baseline indicator line
%     plot([ax(1), ax(2)], [1, 1], ':k');
%     xlim([cfg{1}.epoch.toi.(markername)(1), cfg{1}.epoch.toi.(markername)(2)]);
%     
%     % plot individual units
%     subplot(2, 6, 7);
%     hold on
%     clear y
%     for iunit = 1 : size(GA_bl.(markername).individual, 1)
%         x               = GA_bl.(markername).time;
%         y(iunit,:)      = squeeze(GA_bl.(markername).individual(iunit, 1, :));
%         lh              = plot(x, y(iunit, :), 'k');
%         lh.Color        = [0,0,0,0.2];
%     end
%     plot(x, mean(y, 1),  'k', 'linewidth', 2);
%     title('All');
%     ylabel('Relative firing rate');
%     
%     ymax = 3;
%     ymin = 0;
%     %     for hyplabel = hyplabels
%     %
% %         for iunit = 1 : size(GA_hyp_bl.(markername).(hyplabel).individual, 1)
% %             ymax = max(max(GA_hyp_bl.(markername).(hyplabel).individual(iunit, 1, :)), ymax);
% %             ymin = min(min(GA_hyp_bl.(markername).(hyplabel).individual(iunit, 1, :)), ymin);
% %             fprintf('unit %d, %s, max: %d\n', iunit, hyplabel, ymax);
% %         end
% %     end
%     
%     iplot = 1;
%     for hyplabel = hyplabels
%         subplot(2, 6, 7 + iplot);
%         hold on
%         clear y
%         for iunit = 1 : size(GA_hyp_bl.(markername).(hyplabel).individual, 1)
%             x               = GA_hyp_bl.(markername).(hyplabel).time;
%             y(iunit,:)      = squeeze(GA_hyp_bl.(markername).(hyplabel).individual(iunit, 1, :));
%             lh              = plot(x, y(iunit, :), 'k');
%             lh.Color        = [0,0,0,0.2];
%         end
%         plot(x, nanmean(y, 1),  'k', 'linewidth', 2);
%         iplot = iplot + 1;
%         title(hyplabel, 'interpreter','none');
% %         ylim([ymin, ymax]);
%         ylabel('Relative firing rate');
%     end
%     
%     
%     
%     
%     
%     
    
    
    
    
    
end

