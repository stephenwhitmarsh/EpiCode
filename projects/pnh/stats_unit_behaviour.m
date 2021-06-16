function dat = stats_unit_behaviour(cfg, force)

set(groot,'defaultAxesTickLabelInterpreter','none');
hyplabels   = ["REM",  "AWAKE", "PHASE_1", "PHASE_2", "PHASE_3"];
fname_out   = fullfile(cfg{1}.datasavedir, 'stats_unit_behaviour.mat');

%% load data
if force == false && exist(fname_out, 'file')
    load(fname_out, 'dat');
else
    dat = table;
    for ipatient = 1 : 4
        
        cfg{ipatient}.spike.postfix   = [];
        SpikeTrials_timelocked  = readSpikeTrials_MuseMarkers(cfg{ipatient});
        cfg{ipatient}.spike.postfix   = '_window';
        SpikeTrials_windowed    = readSpikeTrials_MuseMarkers(cfg{ipatient});
        SpikeStats_windowed     = spikeTrialStats(cfg{ipatient});
        SpikeDensity_timelocked = spikeTrialDensity(cfg{ipatient});
        
        for ipart = 1 : size(SpikeTrials_timelocked, 2)
            
            % if no good units were found
            if isempty(SpikeStats_windowed{ipart})
                continue
            end
            
            % loop over units for spikedensity per unit/trial
            for ilabel = 1 : size(SpikeStats_windowed{ipart}.window, 2)
                
                fprintf('Adding: patient %d, part %d, unit %d of %d\n', ipatient, ipart, ilabel, size(SpikeStats_windowed{ipart}.window, 2));
                
                temp                        = SpikeStats_windowed{ipart}.window{ilabel}.trialinfo;
                temp.ipatient               = repmat(ipatient, size(temp, 1), 1);
                temp.ipart                  = repmat(ipart, size(temp, 1), 1);
                temp.ilabel                 = repmat(ilabel, size(temp, 1), 1);
                temp.label                  = repmat(string(SpikeStats_windowed{ipart}.window{ilabel}.label), size(temp, 1), 1);
                temp.RPV                    = repmat(SpikeStats_windowed{ipart}.window{ilabel}.RPV, size(temp, 1), 1);
                temp.trialfreq_corrected    = SpikeStats_windowed{ipart}.window{ilabel}.trialfreq_corrected';
                temp.burst_trialsum         = SpikeStats_windowed{ipart}.window{ilabel}.burst_trialsum';
                temp.CV2_trial              = SpikeStats_windowed{ipart}.window{ilabel}.CV2_trial';
                temp.short                  = SpikeStats_windowed{ipart}.window{ilabel}.short';
                temp.long                   = SpikeStats_windowed{ipart}.window{ilabel}.long';
                temp.amplitude              = SpikeStats_windowed{ipart}.window{ilabel}.amplitude';
                temp.spikecount             = SpikeStats_windowed{ipart}.window{ilabel}.spikecount';
                temp.IEDcount               = zeros(size(temp, 1), 1);
                
                try temp.IEDcount           = temp.IEDcount + temp.combined1; temp = removevars(temp,{'combined1'}); catch; end
                try temp.IEDcount           = temp.IEDcount + temp.combined2; temp = removevars(temp,{'combined2'}); catch; end
                try temp.IEDcount           = temp.IEDcount + temp.combined3; temp = removevars(temp,{'combined3'}); catch; end
                
                % add to windowed
                temp.purity                 = repmat(SpikeTrials_windowed{ipart}.window.purity(ilabel), size(temp, 1), 1);
                temp.cluster_group          = repmat(string(SpikeTrials_windowed{ipart}.window.cluster_group{ilabel}), size(temp, 1), 1);
                
                % check if unit responds statistically
                responsive_pos  = false;
                responsive_neg  = false;
                for markername = string(fields(SpikeDensity_timelocked{ipart}.stat))'
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
                end
                
                temp.responsive     = repmat(responsive_pos | responsive_neg, size(temp, 1), 1);
                temp.responsive_pos = repmat(responsive_pos, size(temp, 1), 1);
                temp.responsive_neg = repmat(responsive_neg, size(temp, 1), 1);
                
                % concatinate
                dat = vertcat(dat, temp);
            end
        end
    end
    
    % save data
    save(fname_out, 'dat', '-v7.3');
end

%% prepare data for stats and plotting
dat.SUA                                         = deblank(dat.cluster_group) == "good";
dat.sleepstage                                  = dat.hyplabel;
dat.sleepstage(dat.sleepstage == "NO_SCORE")    = "AWAKE";
dat.sleepstage                                  = categorical(dat.sleepstage, hyplabels, 'Ordinal', true);
dat.IED(dat.IEDcount ~= 0)                      = "+IEDs";
dat.IED(dat.IEDcount == 0)                      = "-IEDs";
dat.IED                                         = categorical(dat.IED);

% select data
dat_sel                                         = dat(~dat.artefact & dat.SUA & dat.responsive, :);

% calculate averages over time-windows
[G, ipatient, ipart, ilabel, sleepstage, IED]   = findgroups(dat_sel.ipatient, dat_sel.ipart, dat_sel.ilabel, dat_sel.sleepstage, dat_sel.IED);

dat_avg             = table;
dat_avg.ipatient    = ipatient;
dat_avg.ipart       = ipart;
dat_avg.ilabel      = ilabel;
dat_avg.sleepstage  = sleepstage;
dat_avg.IED         = IED;
dat_avg.firingrate  = splitapply(@nanmean, dat_sel.trialfreq_corrected, G);
dat_avg.amplitude   = splitapply(@nanmean, dat_sel.amplitude, G);
dat_avg.CV2         = splitapply(@nanmean, dat_sel.CV2_trial, G);
dat_avg.bursts      = splitapply(@nanmean, dat_sel.burst_trialsum, G);

% calculate averages over sleep stages for normalization
[G, ipatient, ipart, ilabel, IED]                       = findgroups(dat_avg.ipatient, dat_avg.ipart, dat_avg.ilabel, dat_avg.IED);
temp                = table;
temp.ipatient       = ipatient;
temp.ipart          = ipart;
temp.ilabel         = ilabel;
temp.IED            = IED;
temp.firingrate     = splitapply(@nanmean, dat_avg.firingrate, G);
temp.amplitude      = splitapply(@nanmean, dat_avg.amplitude, G);
temp.CV2            = splitapply(@nanmean, dat_avg.CV2, G);
temp.bursts         = splitapply(@nanmean, dat_avg.firingrate, G);

for i = 1 : size(dat_avg, 1)
    indx =  temp.ipatient == dat_avg.ipatient(i) & ...
        temp.ipart == dat_avg.ipart(i) & ...
        temp.ilabel == dat_avg.ilabel(i) & ...
        temp.IED == dat_avg.IED(i);
    dat_avg.firingrate_rel(i)   = dat_avg.firingrate(i) / temp.firingrate(indx);
    dat_avg.amplitude_rel(i)    = dat_avg.amplitude(i) / temp.amplitude(indx);
    dat_avg.CV2_rel(i)          = dat_avg.CV2(i) / temp.CV2(indx);
    dat_avg.bursts_rel(i)       = dat_avg.bursts(i) / temp.bursts(indx);
end

% prepare colormap
cm = cool(5);

%% plot results
fig = figure;
set(gcf, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

settings = {'Notch', 'off', 'GroupByColor', dat_avg.sleepstage, 'MarkerStyle', '.', 'JitterOutliers', 'on'};

subplot(2,4,1);
b = boxchart(dat_avg.IED, log(dat_avg.firingrate), settings{:}); ylabel('log(freq)'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,5);
b = boxchart(dat_avg.IED, dat_avg.firingrate_rel, settings{:}); ylabel('freq'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,2);
b = boxchart(dat_avg.IED, dat_avg.amplitude, settings{:}); ylabel('amplitude'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,6);
b = boxchart(dat_avg.IED, dat_avg.amplitude_rel, settings{:}); ylabel('amplitude (normalized)'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,3);
b = boxchart(dat_avg.IED, dat_avg.CV2, settings{:}); ylabel('CV2'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,7);
b = boxchart(dat_avg.IED, dat_avg.CV2_rel, settings{:}); ylabel('CV2 (normalized)'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,4);
b = boxchart(dat_avg.IED, dat_avg.bursts, settings{:}); ylabel('nr. of bursts'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

subplot(2,4,8);
b = boxchart(dat_avg.IED, dat_avg.bursts_rel, settings{:}); ylabel('nr. of bursts (normalized)'); legend('interpreter','none');
for k = 1 : size(b, 1); b(k).BoxFaceColor = cm(k, :); end; for k = 1 : size(b, 1); b(k).MarkerColor = cm(k, :); end

% print file
fname = fullfile(cfg{1}.imagesavedir, 'stats', 'stats_windowed');
disp(['Exporting figure ', fname])
exportgraphics(fig, [fname, '.pdf']);
exportgraphics(fig, [fname, '.tiff']);

%% statistics

% Multinomial logistic regression (ordinal response variable) shows that
% CV2 and nr. of bursts (weakly), but not FR and amplitute, explains whether there are IEDs
[B, dev, stats] = mnrfit([dat_avg.firingrate, dat_avg.CV2, dat_avg.amplitude, dat_avg.bursts], dat_avg.IED, 'model', 'nominal')

% Multinomial logistic regression (ordinal response variable) shows that
% CV2 and nr. of bursts (but not FR or amplitude) explain sleep stage when looking
% at time periods with IEDs, but not in those without
[B, dev, stats] = mnrfit([dat_avg.firingrate(dat_avg.IED == "-IEDs", :), dat_avg.CV2(dat_avg.IED == "-IEDs", :), dat_avg.amplitude(dat_avg.IED == "-IEDs", :), dat_avg.bursts(dat_avg.IED == "-IEDs", :)], dat_avg.sleepstage(dat_avg.IED == "-IEDs", :), 'model', 'ordinal');
[B, dev, stats] = mnrfit([dat_avg.firingrate(dat_avg.IED == "+IEDs", :), dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.amplitude(dat_avg.IED == "+IEDs", :), dat_avg.bursts(dat_avg.IED == "+IEDs", :)], dat_avg.sleepstage(dat_avg.IED == "+IEDs", :), 'model', 'ordinal');

% This shows that mnrfit reverses category order
[B, dev, stats] = mnrfit(dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.sleepstage(dat_avg.IED == "+IEDs", :), 'model', 'ordinal');
[B, dev, stats] = glmfit(dat_avg.CV2(dat_avg.IED == "+IEDs", :), dat_avg.test(dat_avg.IED == "+IEDs", :));

% generalized linear mixed model shows that the number of IEDs is explained
% by sleep stage, CV2, unit amplitude, nr. of bursts, and firingrate.
glme = fitglme(dat(dat.SUA & dat.IED == "+IEDs", :),...
    'IEDcount ~ 1 + hyplabel + trialfreq_corrected + CV2_trial + amplitude + burst_trialsum + (1|ilabel)');
anova(glme)