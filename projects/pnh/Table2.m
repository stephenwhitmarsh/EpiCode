function T = Table2(cfg, force)

set(groot,'defaultAxesTickLabelInterpreter','none');
hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];
fname_out   = fullfile(cfg{1}.datasavedir, 'stats_unit_behaviour.mat');

%% load data
if force == false && exist(fname_out, 'file')
    load(fname_out, 'dat_avg');
    return
end

corr.channel{1}.PSW     = 'm1pNs_6';
corr.channel{1}.FA      = 'm1pNs_4';
corr.channel{1}.ES      = 'm1pNs_4';
corr.channel{2}.PSW     = 'mCasd_2';
corr.channel{2}.FA      = 'mCasd_2';
corr.channel{2}.ES      = 'mCasd_2';
corr.channel{3}.PSW     = 'mTNmi_3';
corr.channel{3}.FA      = 'mTNmi_3';
corr.channel{3}.ES      = 'mTNmi_3';
corr.channel{4}.PSW     = 'mLMI1_7';
corr.channel{4}.FA      = 'mLMI1_7';
corr.channel{4}.ES      = 'mLMI1_7';


dat = table;
for ipatient = 1 : size(cfg, 2)

    cfg{ipatient}.spike.postfix     = [];
    SpikeTrials_timelocked          = readSpikeTrials_MuseMarkers(cfg{ipatient});
    cfg{ipatient}.spike.postfix     = '-windowed';
    SpikeStats_windowed             = spikeTrialStats(cfg{ipatient});
    cfg{ipatient}.spike.postfix     = '';
    SpikeDensity_timelocked         = spikeTrialDensity(cfg{ipatient});
    
    for ipart = 1 : size(SpikeTrials_timelocked, 2)
        
        % if no good units were found
        if isempty(SpikeStats_windowed{ipart})
            continue
        end
        
        % loop over units for spikedensity per unit/trial
        for itemp = 1 : size(SpikeStats_windowed{ipart}.window, 2)
            
            fprintf('Adding: patient %d, part %d, unit %d of %d\n', ipatient, ipart, itemp, size(SpikeStats_windowed{ipart}.window, 2));
            
            temp                        = SpikeStats_windowed{ipart}.window{itemp}.trialinfo;
            temp.ipatient               = repmat(ipatient, size(temp, 1), 1);
            temp.ipart                  = repmat(ipart, size(temp, 1), 1);
            temp.itemp                  = repmat(itemp, size(temp, 1), 1);
            temp.label                  = repmat(string(SpikeStats_windowed{ipart}.window{itemp}.label), size(temp, 1), 1);
            temp.RPV                    = repmat(SpikeStats_windowed{ipart}.window{itemp}.RPV, size(temp, 1), 1);
            temp.trialfreq_corrected    = SpikeStats_windowed{ipart}.window{itemp}.trialfreq_corrected';
            temp.burst_trialsum         = SpikeStats_windowed{ipart}.window{itemp}.burst_trialsum';
            temp.CV2_trial              = SpikeStats_windowed{ipart}.window{itemp}.CV2_trial';
            temp.short                  = SpikeStats_windowed{ipart}.window{itemp}.short';
            temp.long                   = SpikeStats_windowed{ipart}.window{itemp}.long';
            temp.amplitude              = SpikeStats_windowed{ipart}.window{itemp}.amplitude';
            temp.spikecount             = SpikeStats_windowed{ipart}.window{itemp}.spikecount';
            
            % add to windowed (take first marker, data is identical)
            fn = fields(SpikeTrials_timelocked{ipart});
            temp.purity                 = repmat(SpikeTrials_timelocked{ipart}.(fn{1}).purity(itemp), size(temp, 1), 1);
            temp.cluster_group          = repmat(string(SpikeTrials_timelocked{ipart}.(fn{1}).cluster_group{itemp}), size(temp, 1), 1);
            
            % check if unit responds statistically
            responsive_pos  = false;
            responsive_neg  = false;
            for markername = string(fields(SpikeDensity_timelocked{ipart}.stat))'
                
                % add correlation with LFP
                if isfield(corr.channel{ipatient}, markername)
                    indx = ismember(cfg{ipatient}.LFP.channel, corr.channel{ipatient}.(markername));
                    if any(indx)
                        eval(sprintf('temp.%s_LFPcorr_rho = repmat(SpikeDensity_timelocked{ipart}.psth.(markername).corr_rho(itemp,  indx), size(temp, 1), 1);', markername));
                        eval(sprintf('temp.%s_LFPcorr_p   = repmat(SpikeDensity_timelocked{ipart}.psth.(markername).corr_pval(itemp, indx), size(temp, 1), 1);', markername));
                    else
                        eval(sprintf('temp.%s_LFPcorr_rho = repmat(nan, size(temp, 1), 1);', markername));
                        eval(sprintf('temp.%s_LFPcorr_p   = repmat(nan, size(temp, 1), 1);', markername));
                    end
                else
                    eval(sprintf('temp.%s_LFPcorr_rho = repmat(nan, size(temp, 1), 1);', markername));
                    eval(sprintf('temp.%s_LFPcorr_p   = repmat(nan, size(temp, 1), 1);', markername));
                end
                
                % add signification increase/decrease for each pattern
                if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'posclusters') &&  ~isempty(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters)
                    clusterindx = find(vertcat(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters.prob) < 0.025, 1, 'first');
                    if ~isempty(clusterindx)
                        responsive_pos  = true;
                        shift           = size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusterslabelmat, 2);
                        mask            = [false(1, shift), ismember(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusterslabelmat, clusterindx)];
                        Y               = max(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, mask));
                        baseline        = mean(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.baseline(:, itemp));
                        p               = SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.posclusters(clusterindx).prob;
                        val             = Y / baseline * 100 - 100;
                        eval(sprintf('temp.%s_increase_val = repmat(val, size(temp, 1), 1); ', markername));
                        eval(sprintf('temp.%s_increase_p   = repmat(p,   size(temp, 1), 1); ', markername));
                    else
                        eval(sprintf('temp.%s_increase_val = repmat(nan, size(temp, 1), 1); ', markername));
                        eval(sprintf('temp.%s_increase_p   = repmat(nan, size(temp, 1), 1); ', markername));
                    end
                else
                    eval(sprintf('temp.%s_increase_val = repmat(nan, size(temp, 1), 1); ', markername));
                    eval(sprintf('temp.%s_increase_p   = repmat(nan, size(temp, 1), 1); ', markername));
                end
                
                if isfield(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}, 'negclusters') && ~isempty(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters)
                    clusterindx = find(vertcat(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters.prob) < 0.025, 1, 'first');
                    if ~isempty(clusterindx)
                        responsive_neg  = true;
                        shift           = size(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg, 2) - size(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusterslabelmat, 2);
                        mask            = [false(1, shift), ismember(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusterslabelmat, clusterindx)];
                        Y               = max(SpikeDensity_timelocked{ipart}.sdf_bar.(markername).avg(itemp, mask));
                        baseline        = mean(SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.baseline(:, itemp));
                        p               = SpikeDensity_timelocked{ipart}.stat.(markername){itemp}.negclusters(clusterindx).prob;
                        val             = Y / baseline * 100 - 100;
                        eval(sprintf('temp.%s_decrease_val = repmat(val, size(temp, 1), 1); ', markername));
                        eval(sprintf('temp.%s_decrease_p   = repmat(p,   size(temp, 1), 1); ', markername));
                    else
                        eval(sprintf('temp.%s_decrease_val = repmat(nan, size(temp, 1), 1); ', markername));
                        eval(sprintf('temp.%s_decrease_p   = repmat(nan, size(temp, 1), 1); ', markername));
                    end
                else
                    eval(sprintf('temp.%s_decrease_val = repmat(nan, size(temp, 1), 1); ', markername));
                    eval(sprintf('temp.%s_decrease_p   = repmat(nan, size(temp, 1), 1); ', markername));
                end
            end
            
            temp.responsive     = repmat(responsive_pos | responsive_neg, size(temp, 1), 1);
            temp.responsive_pos = repmat(responsive_pos, size(temp, 1), 1);
            temp.responsive_neg = repmat(responsive_neg, size(temp, 1), 1);
            
            % concatinate over parts/patients
            temp_colmissing = setdiff(dat.Properties.VariableNames, temp.Properties.VariableNames);
            if ~isempty(temp_colmissing)
                temp = [temp array2table(nan(height(temp), numel(temp_colmissing)), 'VariableNames', temp_colmissing)];
            end
            dat = vertcat(dat, temp);
        end
    end
end

%% prepare data for stats and plotting
dat.MUA                                         = deblank(dat.cluster_group) == "mua"; % this way SUA end up before MUA (false-true)
dat.sleepstage                                  = dat.hyplabel;
dat.sleepstage(dat.sleepstage == "NO_SCORE")    = "AWAKE";
dat.sleepstage                                  = categorical(dat.sleepstage, hyplabels, 'Ordinal', true);
dat.IED(dat.overlapcount ~= 0)                  = "+IEDs";
dat.IED(dat.overlapcount == 0)                  = "-IEDs";
dat.IED                                         = categorical(dat.IED);

% select data
% dat_sel                                         = dat(~dat.artefact & dat.SUA & dat.responsive, :);
dat_sel                                         = dat(~dat.artefact & ~dat.overlapcount, :);

% calculate averages over time-windows
% [G, ipatient, ipart, itemp, sleepstage, IED]   = findgroups(dat_sel.ipatient, dat_sel.ipart, dat_sel.itemp, dat_sel.sleepstage, dat_sel.IED);
[G, ipatient, MUA, itemp]   = findgroups(dat_sel.ipatient, dat_sel.MUA, dat_sel.itemp);

t             = table;
t.ipatient    = ipatient;
t.itemp       = itemp;
t.MUA         = MUA;

t.firingrate  = splitapply(@nanmean, dat_sel.trialfreq_corrected, G);
t.amplitude   = splitapply(@nanmean, dat_sel.amplitude, G);
t.CV2         = splitapply(@nanmean, dat_sel.CV2_trial, G);
t.RPV         = splitapply(@nanmean, dat_sel.RPV, G);
t.bursts      = splitapply(@nanmean, dat_sel.burst_trialsum, G);

% nanmean is just a trick here to get the value (which is same in selection)
t.PSW_inc_p   = splitapply(@nanmean, dat_sel.PSW_increase_p, G);
t.PSW_inc_val = splitapply(@nanmean, dat_sel.PSW_increase_val, G);
t.PSW_dec_p   = splitapply(@nanmean, dat_sel.PSW_decrease_p, G);
t.PSW_dec_val = splitapply(@nanmean, dat_sel.PSW_decrease_val, G);

t.FA_inc_p   = splitapply(@nanmean, dat_sel.FA_increase_p, G);
t.FA_inc_val = splitapply(@nanmean, dat_sel.FA_increase_val, G);
t.FA_dec_p   = splitapply(@nanmean, dat_sel.FA_decrease_p, G);
t.FA_dec_val = splitapply(@nanmean, dat_sel.FA_decrease_val, G);

t.ES_inc_p   = splitapply(@nanmean, dat_sel.ES_increase_p, G);
t.ES_inc_val = splitapply(@nanmean, dat_sel.ES_increase_val, G);
t.ES_dec_p   = splitapply(@nanmean, dat_sel.ES_decrease_p, G);
t.ES_dec_val = splitapply(@nanmean, dat_sel.ES_decrease_val, G);

t.PSW_cor_p   = splitapply(@nanmean, dat_sel.PSW_LFPcorr_p, G);
t.PSW_cor_val = splitapply(@nanmean, dat_sel.PSW_LFPcorr_rho, G);

t.FA_cor_p   = splitapply(@nanmean, dat_sel.FA_LFPcorr_p, G);
t.FA_cor_val = splitapply(@nanmean, dat_sel.FA_LFPcorr_rho, G);

t.ES_cor_p   = splitapply(@nanmean, dat_sel.ES_LFPcorr_p, G);
t.ES_cor_val = splitapply(@nanmean, dat_sel.ES_LFPcorr_rho, G);



% Convert to a cell array of strings
T = table;
fun0 = @(x) sprintf('%0.0f', x);
fun1 = @(x) sprintf('%0.1f', x);
fun2 = @(x) sprintf('%0.2f', x);
fun3 = @(x) sprintf('%0.3f', x);

T.N                 = t.ipatient;
T.U                 = t.itemp;
T.type(t.MUA, :)    = "MUA";
T.type(~t.MUA, :)   = "SUA";

T.FR    = cellfun(fun1, num2cell(round(t.firingrate, 1)), 'UniformOutput',0);
T.Amp   = cellfun(fun1, num2cell(round(t.amplitude, 1)), 'UniformOutput',0);
T.RPV   = cellfun(fun1, num2cell(round(t.RPV*100, 1)), 'UniformOutput',0);
T.CV2   = cellfun(fun2, num2cell(round(t.CV2, 2)), 'UniformOutput',0);

T.PSWi  = cellfun(fun0, num2cell(round(t.PSW_inc_val)), 'UniformOutput',0);
T.PSWd  = cellfun(fun0, num2cell(abs(round(t.PSW_dec_val))), 'UniformOutput',0);
T.PSWr  = cellfun(fun2, num2cell(round(t.PSW_cor_val, 2)), 'UniformOutput',0);

T.FAi   = cellfun(fun0, num2cell(round(t.FA_inc_val)), 'UniformOutput',0);
T.FAd   = cellfun(fun0, num2cell(abs(round(t.FA_dec_val))), 'UniformOutput',0);
T.FAr   = cellfun(fun2, num2cell(round(t.FA_cor_val, 2)), 'UniformOutput',0);

T.ESi   = cellfun(fun0, num2cell(round(t.ES_inc_val)), 'UniformOutput',0);
T.ESd   = cellfun(fun0, num2cell(abs(round(t.ES_dec_val))), 'UniformOutput',0);
T.ESr   = cellfun(fun2, num2cell(round(t.ES_cor_val, 2)), 'UniformOutput',0);

T.PSWp  = cellfun(fun3, num2cell(round(t.PSW_cor_p, 3)), 'UniformOutput',0);
T.FAp   = cellfun(fun3, num2cell(round(t.FA_cor_p, 3)), 'UniformOutput',0);
T.ESp   = cellfun(fun3, num2cell(round(t.ES_cor_p, 3)), 'UniformOutput',0);

T.Properties.VariableNames(T.Properties.VariableNames == "PSWi") = "PSW+";
T.Properties.VariableNames(T.Properties.VariableNames == "PSWd") = "PSW-";
T.Properties.VariableNames(T.Properties.VariableNames == "FAi") = "FA+";
T.Properties.VariableNames(T.Properties.VariableNames == "FAd") = "FA-";
T.Properties.VariableNames(T.Properties.VariableNames == "ESi") = "ES+";
T.Properties.VariableNames(T.Properties.VariableNames == "ESd") = "ES-";

T.Properties.VariableNames(T.Properties.VariableNames == "PSWr") = sprintf("PSW%c", 961);
T.Properties.VariableNames(T.Properties.VariableNames == "FAr") = sprintf("FA%c", 961);
T.Properties.VariableNames(T.Properties.VariableNames == "ESr") = sprintf("ES%c", 961);

for fn = string(T.Properties.VariableNames)
    i = strcmp((T.(fn)), "NaN");
    if any(i)
        T.(fn)(i,:) = {''};
    end
end

fname = fullfile(cfg{1}.imagesavedir, 'Table2.xls');
writetable(T, fname);
disp('done');













% 
% 
% % calculate averages over sleep stages for normalization
% [G, ipatient, ipart, itemp, IED]                       = findgroups(dat_avg.ipatient, dat_avg.ipart, dat_avg.itemp, dat_avg.IED);
% temp                = table;
% temp.ipatient       = ipatient;
% temp.ipart          = ipart;
% temp.itemp         = itemp;
% temp.IED            = IED;
% temp.firingrate     = splitapply(@nanmean, dat_avg.firingrate, G);
% temp.amplitude      = splitapply(@nanmean, dat_avg.amplitude, G);
% temp.CV2            = splitapply(@nanmean, dat_avg.CV2, G);
% temp.bursts         = splitapply(@nanmean, dat_avg.firingrate, G);
% 
% for i = 1 : size(dat_avg, 1)
%     indx =  temp.ipatient == dat_avg.ipatient(i) & ...
%         temp.ipart == dat_avg.ipart(i) & ...
%         temp.itemp == dat_avg.itemp(i) & ...
%         temp.IED == dat_avg.IED(i);
%     dat_avg.firingrate_rel(i)   = dat_avg.firingrate(i) / temp.firingrate(indx);
%     dat_avg.amplitude_rel(i)    = dat_avg.amplitude(i) / temp.amplitude(indx);
%     dat_avg.CV2_rel(i)          = dat_avg.CV2(i) / temp.CV2(indx);
%     dat_avg.bursts_rel(i)       = dat_avg.bursts(i) / temp.bursts(indx);
% end
% 
% % save data
% save(fname_out, 'dat', 'dat_avg', '-v7.3');
% 
% % prepare colormap
% cm = cool(5);
% 
% %% plot results
% fig = figure;
% set(gcf, 'position', get(0,'ScreenSize'));
% set(fig, 'PaperOrientation', 'landscape');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0 0 1 1]);
% set(fig, 'Renderer', 'Painters');
% 
% settings = {'Notch', 'off', 'GroupByColor', dat_avg.sleepstage, 'MarkerStyle', '.', 'JitterOutliers', 'on'};
% 
% subplot(2,4,1);
% b = boxchart(dat_avg.IED, log(dat_avg.firingrate), settings{:}); ylabel('log(freq)'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,5);
% b = boxchart(dat_avg.IED, dat_avg.firingrate_rel, settings{:}); ylabel('freq'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,2);
% b = boxchart(dat_avg.IED, dat_avg.amplitude, settings{:}); ylabel('amplitude'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,6);
% b = boxchart(dat_avg.IED, dat_avg.amplitude_rel, settings{:}); ylabel('amplitude (normalized)'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,3);
% b = boxchart(dat_avg.IED, dat_avg.CV2, settings{:}); ylabel('CV2'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,7);
% b = boxchart(dat_avg.IED, dat_avg.CV2_rel, settings{:}); ylabel('CV2 (normalized)'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,4);
% b = boxchart(dat_avg.IED, dat_avg.bursts, settings{:}); ylabel('nr. of bursts'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% subplot(2,4,8);
% b = boxchart(dat_avg.IED, dat_avg.bursts_rel, settings{:}); ylabel('nr. of bursts (normalized)'); legend('interpreter','none');
% for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end
% 
% % print file
% fname = fullfile(cfg{1}.imagesavedir, 'stats', 'stats_windowed');
% disp(['Exporting figure ', fname])
% exportgraphics(fig, [fname, '.pdf']);
% exportgraphics(fig, [fname, '.tiff']);
% 
% %% statistics
% try
% % Multinomial logistic regression (ordinal response variable) shows that
% % CV2 and nr. of bursts (weakly), but not FR and amplitute, explains whether there are IEDs
% [B, dev, stats] = mnrfit([dat_avg.firingrate, dat_avg.CV2, dat_avg.amplitude, dat_avg.bursts], dat_avg.IED, 'model', 'nominal')
% 
% % Multinomial logistic regression (ordinal response variable) shows that
% % CV2 and nr. of bursts (but not FR or amplitude) explain sleep stage when looking
% % at time periods with IEDs, but not in those without
% [B, dev, stats] = mnrfit([dat_avg.firingrate(dat_avg.IED == "-IEDs", :), dat_avg.CV2(dat_avg.IED == "-IEDs", :), dat_avg.amplitude(dat_avg.IED == "-IEDs", :), dat_avg.bursts(dat_avg.IED == "-IEDs", :)], dat_avg.sleepstage(dat_avg.IED == "-IEDs", :), 'model', 'ordinal');
% [B, dev, stats] = mnrfit([dat_avg.firingrate(dat_avg.IED == "+IEDs", :), dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.amplitude(dat_avg.IED == "+IEDs", :), dat_avg.bursts(dat_avg.IED == "+IEDs", :)], dat_avg.sleepstage(dat_avg.IED == "+IEDs", :), 'model', 'ordinal');
% 
% % This shows that mnrfit reverses category order
% [B, dev, stats] = mnrfit(dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.sleepstage(dat_avg.IED == "+IEDs", :), 'model', 'ordinal');
% [B, dev, stats] = glmfit(dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.test(dat_avg.IED == "+IEDs", :));
% 
% % generalized linear mixed model shows that the number of IEDs is explained
% % by sleep stage, CV2, unit amplitude, nr. of bursts, and firingrate.
% glme = fitglme(dat(dat.SUA & dat.IED == "+IEDs", :),...
%     'IEDcount ~ 1 + hyplabel + trialfreq_corrected + CV2_trial + amplitude + burst_trialsum + (1|itemp)');
% anova(glme)
% 
% catch
% end