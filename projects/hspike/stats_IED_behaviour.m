function [dat_unit, dat_lfp] = stats_IED_behaviour(cfg, force)

set(groot,'defaultAxesTickLabelInterpreter','none');
hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];
fname_out   = fullfile(cfg{1}.datasavedir, 'stats_IED_behaviour.mat');

%% load data
if force == false && exist(fname_out, 'file')
    load(fname_out, 'dat_unit', 'dat_lfp');
    
else
    
    % per patient/night/IED x unit (so not always present, but often
    % multiple)
    dat_unit = table;
    % patient/night/IED (so all)
    dat_lfp = table;
    
    for ipatient = 1 : 7
        
        cfg{ipatient}.postfix   = [];
        SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg{ipatient});
        SpikeDensity_timelocked = spikeTrialDensity(cfg{ipatient});
        cfg{ipatient}.LFP.name  = string(fields(SpikeTrials_timelocked{1})');
        LFP                     = readLFP(cfg{ipatient});
        
        for ipart = 1 : size(SpikeTrials_timelocked, 2)
            
            for markername = string(fields(LFP{ipart}))'
                
                % calculate max/min/RMS per trial
                cfgtemp         = [];
                cfgtemp.latency = [-0.2 0.5]; %cfg{ipatient}.align.latency
                dat_sel         = ft_selectdata(cfgtemp, LFP{ipart}.(markername));
                t               = reshape(cell2mat(dat_sel.trial), size(dat_sel.trial{1}, 1), size(dat_sel.trial{1}, 2), size(dat_sel.trial, 2));
                maxpeak         = max(t, [], 2);
                maxpeak         = max(maxpeak, [], 1);
                maxpeak         = squeeze(maxpeak);
                minpeak         = min(t, [], 2);
                minpeak         = min(minpeak, [], 1);
                minpeak         = squeeze(minpeak);
                
                clear dat_sel t
                LFP{ipart}.(markername).trialinfo.maxpeak    = maxpeak;
                LFP{ipart}.(markername).trialinfo.minpeak    = minpeak;
                LFP{ipart}.(markername).trialinfo.markername = repmat(markername, size(LFP{ipart}.(markername).trialinfo, 1), 1);
                LFP{ipart}.(markername).trialinfo.ipart      = repmat(ipart, size(LFP{ipart}.(markername).trialinfo, 1), 1);
                LFP{ipart}.(markername).trialinfo.ipatient   = repmat(ipatient, size(LFP{ipart}.(markername).trialinfo, 1), 1);
                
                dat_lfp         = vertcat(dat_lfp, LFP{ipart}.(markername).trialinfo);
                
                % if no good units were found
                if isempty(SpikeTrials_timelocked{ipart})
                    continue
                end
                
                % loop over units for spikedensity per unit/trial
                for ilabel = 1 : size(SpikeTrials_timelocked{ipart}.(markername).label, 2)
                    
                    fprintf('Adding: patient %d, part %d, %s, unit %d of %d\n', ipatient, ipart, markername, ilabel, size(SpikeTrials_timelocked{ipart}.(markername).label, 2));
                    
                    temp                        = SpikeTrials_timelocked{ipart}.(markername).trialinfo;
                    
                    temp.ipatient               = repmat(ipatient, size(temp, 1), 1);
                    temp.ipart                  = repmat(ipart, size(temp, 1), 1);
                    temp.ilabel                 = repmat(ilabel, size(temp, 1), 1);
                    temp.label                  = repmat(string(SpikeTrials_timelocked{ipart}.(markername).label{ilabel}), size(temp, 1), 1);
                    temp.markername             = repmat(markername, size(temp, 1), 1);
                    temp.cluster_group          = repmat(deblank(string(SpikeTrials_timelocked{ipart}.combined1.cluster_group{ilabel})), size(temp, 1), 1);
                    
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
                    
                    % map trials according to directory & trialnr.
                    temp.Properties.VariableNames{'begsample'}      = 'begsample_unit';
                    temp.Properties.VariableNames{'endsample'}      = 'endsample_unit';
                    temp.Properties.VariableNames{'starttime'}      = 'starttime_unit';
                    temp.Properties.VariableNames{'endtime'}        = 'endtime_unit';
                    temp.Properties.VariableNames{'offset'}         = 'offset_unit';
                    temp.Properties.VariableNames{'trialnr_start'}  = 'trialnr_unit';
                    temp.Properties.VariableNames{'idir'}           = 'idir_unit';
                    temp.Properties.VariableNames{'fileoffset'}     = 'fileoffset_unit';
                    temp.Properties.VariableNames{'hyplabel'}       = 'hyplabel_unit';
                    
                    LFP{ipart}.(markername).trialinfo.directory     = string(LFP{ipart}.(markername).trialinfo.directory);
                    temp.directory                                  = string(string(temp.directory));
                    
                    for i = 1 : size(temp, 1)
                        indx_dir = strcmp(temp.directory(i, :), LFP{ipart}.(markername).trialinfo.directory);
                        indx_trl = temp.trialnr_unit(i) == LFP{ipart}.(markername).trialinfo.trialnr;
                        indx = indx_trl & indx_dir;
                        
                        if sum(indx) ~= 1
                            fprintf('Cannot find 1:1 correspondance between spike and LFP data for trial indx %d\n', i);
                            temp.consistent(i)     = false;
                            temp.maxpeak_lfp(i)    = nan;
                            temp.minpeak_lfp(i)    = nan;
                            temp.hyplabel_lfp(i)   = "";
                            temp.begsample_lfp(i)  = NaN;
                            temp.starttime_lfp(i)  = NaT;
                            temp.endtime_lfp(i)    = NaT;
                            temp.directory_lfp(i)  = "";
                        else
                            temp.maxpeak_lfp(i)    = maxpeak(indx);
                            temp.minpeak_lfp(i)    = minpeak(indx);
                            temp.hyplabel_lfp(i)   = LFP{ipart}.(markername).trialinfo.hyplabel(indx);
                            temp.begsample_lfp(i)  = LFP{ipart}.(markername).trialinfo.begsample(indx);
                            temp.starttime_lfp(i)  = LFP{ipart}.(markername).trialinfo.starttime(indx);
                            temp.endtime_lfp(i)    = LFP{ipart}.(markername).trialinfo.endtime(indx);
                            temp.directory_lfp(i)  = LFP{ipart}.(markername).trialinfo.directory(indx);
                            if temp.hyplabel_unit(i) == temp.hyplabel_lfp(i)
                                temp.consistent(i)                        = true;
                            else
                                fprintf('Inconsitent hypnogram label spike and LFP data for trial indx %d, timedifference = %1.0f milliseconds \n', i,  (seconds(temp.starttime_lfp(i) - temp.starttime_unit(i)))*1000);
                                temp.consistent(i)                        = false;
                            end
                        end
                    end
                    
                    % concatinate
                    try
                        dat_unit = vertcat(dat_unit, temp);
                    catch
                        disp('meh....');
                    end
                end
            end
        end
    end
    
    % save data
    save(fname_out, 'dat_unit', 'dat_lfp', '-v7.3');
end

%% prepare data for stats and plotting
dat_unit.SUA                                    = dat_unit.cluster_group == "good";
dat_unit.sleepstage                             = dat_unit.hyplabel_unit;
dat_unit(dat_unit.sleepstage == "NO_SCORE", :)  = [];
dat_unit.sleepstage                             = categorical(dat_unit.sleepstage, hyplabels, 'Ordinal', true);

dat_lfp.sleepstage                              = dat_lfp.hyplabel;
dat_lfp(dat_lfp.sleepstage == "NO_SCORE", :)    = [];
dat_lfp.sleepstage                              = categorical(dat_lfp.sleepstage, hyplabels, 'Ordinal', true);

dat_lfp                                         = sortrows(dat_lfp, 'starttime');
dat_lfp.ISI                                     = NaN(size(dat_lfp, 1), 1);
dat_lfp.ISI(2:end, :)                           = seconds(diff(dat_lfp.starttime));

% load hynogram
dat_hyp = table;
for ipatient = 1 : 7
    [~, hyp, ~]     = hypnogramMuseStats(cfg{ipatient});
    hyp.ipatient    = repmat(ipatient, size(hyp, 1), 1);
    hyp(hyp.hyplabel == "NO_SCORE", :) = [];
    hyp.sleepstage  = categorical(hyp.hyplabel, hyplabels);
    dat_hyp         = vertcat(dat_hyp, hyp);
end

% create summary of IED rate and hypnogram
summary             = table;
[G, summary.ipatient, summary.ipart, summary.sleepstage] = findgroups(dat_hyp.ipatient, dat_hyp.part, dat_hyp.sleepstage);
summary.duration    = splitapply(@sum, dat_hyp.duration, G);

[G, ~, ~ ,~]        = findgroups(dat_lfp.ipatient, dat_lfp.ipart, dat_lfp.sleepstage);
summary.maxpeak     = splitapply(@nanmean, dat_lfp.maxpeak(~dat_lfp.artefact), G(~dat_lfp.artefact));
summary.minpeak     = splitapply(@nanmean, dat_lfp.minpeak(~dat_lfp.artefact), G(~dat_lfp.artefact));
summary.count       = groupcounts({dat_lfp.ipatient, dat_lfp.ipart, dat_lfp.sleepstage});
summary.IEDrate     = summary.count ./ minutes(summary.duration);

% normalize to average over sleepstages
for i = 1 : size(summary, 1)
    indx =  summary.ipatient == summary.ipatient(i) & ...
        summary.ipart == summary.ipart(i);
    summary.duration_rel(i) = summary.duration(i) / mean(minutes(summary.duration(indx)));
    summary.count_rel(i)    = summary.count(i)    / mean(summary.count(indx));
    summary.IEDrate_rel(i)  = summary.IEDrate(i)  / mean(summary.IEDrate(indx));
    summary.maxpeak_rel(i)  = summary.maxpeak(i)  / mean(summary.maxpeak(indx));
    summary.minpeak_rel(i)  = summary.minpeak(i)  / mean(summary.minpeak(indx));
end

%% plot results
cm = cool(5);

fig = figure;
set(gcf, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

settings = {'Notch', 'on', 'GroupByColor', summary.sleepstage, 'MarkerStyle', '.', 'JitterOutliers', 'on'};

subplot(2,5,1);
b = boxchart(summary.IEDrate, settings{:}); ylabel('IED per minute'); title('IED rate'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,2);
b = boxchart(minutes(summary.duration), settings{:}); ylabel('Minutes'); title('Total time'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,3);
b = boxchart(summary.count, settings{:}); title('Total nr. of IEDs'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,4);
b = boxchart(summary.maxpeak, settings{:}); ylabel('uV'); title('Peak max LFP'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,5);
b = boxchart(summary.minpeak, settings{:}); ylabel('uV'); title('Peak min LFP '); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,6);
b = boxchart(summary.IEDrate_rel, settings{:}); ylabel('Per minute (relative)'); title('IED rate'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,7);
b = boxchart(minutes(summary.duration_rel), settings{:}); ylabel('Minutes (relative to mean)'); title('Total time'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,8);
b = boxchart(summary.count_rel, settings{:}); ylabel('relative to mean'); title('Total nr. of IEDs'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,9);
b = boxchart(summary.maxpeak_rel, settings{:}); ylabel('Relative'); title('Peak max LFP'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,5,10);
b = boxchart(summary.minpeak_rel, settings{:}); ylabel('Relative'); title('Peak min LFP'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

% print file
fname = fullfile(cfg{1}.imagesavedir, 'stats', 'stats_timelocked');
disp(['Exporting figure ', fname])
exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff']);
close all

%% statistics

% Anova (patient x night = 21 repetitions) shows clearly that IED rate
% increases with increasing sleepstage
[p, anovatab, stats] = anova1(summary.IEDrate, summary.sleepstage);
[p, anovatab, stats] = anova1(summary.maxpeak, summary.sleepstage);
[p, anovatab, stats] = anova1(summary.minpeak, summary.sleepstage);
[p, anovatab, stats] = anova1(summary.IEDrate_rel, summary.sleepstage);
[p, anovatab, stats] = anova1(summary.maxpeak_rel, summary.sleepstage);
[p, anovatab, stats] = anova1(summary.minpeak_rel, summary.sleepstage);

% Non-parametric Kruskal-Wallis (patient x night = 21 repetitions) shows clearly that IED rate
% increases with increasing sleepstage
[p, anovatab, stats] = kruskalwallis(summary.IEDrate, summary.sleepstage)
[p, anovatab, stats] = kruskalwallis(summary.maxpeak, summary.sleepstage)
[p, anovatab, stats] = kruskalwallis(summary.minpeak, summary.sleepstage)
[p, anovatab, stats] = kruskalwallis(summary.IEDrate_rel, summary.sleepstage)
[p, anovatab, stats] = kruskalwallis(summary.maxpeak_rel, summary.sleepstage)
[p, anovatab, stats] = kruskalwallis(summary.minpeak_rel, summary.sleepstage)

% Multinomial logistic regression (ordinal response variable; 87364 rows) shows that
% IED rate (determined by ISI) explains sleepstage, as well as max/min peak
dat_lfp2 = dat_lfp(2:end,:);
[B, dev, stats] = mnrfit([dat_lfp2.ISI, dat_lfp2.maxpeak, dat_lfp2.minpeak], dat_lfp2.sleepstage, 'model', 'ordinal')

% generalized linear mixed model (87364 rows) shows that the number of IEDs is explained
% by sleep stage.
glme = fitglme(dat_lfp,...
    'ISI ~ 1 + hyplabel + (1|markername) + (1|ipatient) + (1|ipart)');
anova(glme)