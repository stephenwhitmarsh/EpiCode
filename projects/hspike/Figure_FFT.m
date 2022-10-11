function [config] = Figure_FFT

set(0,'defaulttextinterpreter','none')
%% load data
for ipatient = 1:8
    config                          = hspike_setparams;
    config{ipatient}.FFT.name       = {'window'};
    config{ipatient}.FFT.postfix    = {'_noWelch'};    
    FFT{ipatient}                   = FFTtrials(config{ipatient}, false);
end

%% average over channels - all stages
cfgtemp = [];
cfgtemp.avgoverchan = 'yes';
cfgtemp.avgoverrpt = 'yes';
for ipatient = 1:8
    for ipart = 1 : 3
        FFTavg_part{ipart} = ft_selectdata(cfgtemp, FFT{ipatient}{ipart}.window);      
    end
    FFTavg{ipatient} = FFTavg_part{1};
    FFTavg{ipatient}.powspctrm = nanmean([FFTavg_part{1}.powspctrm; FFTavg_part{2}.powspctrm; FFTavg_part{3}.powspctrm]);   
end

figure;
hold;
for ipatient = 1 : 8
    plot(FFTavg{ipatient}.freq, FFTavg{ipatient}.powspctrm );
end

%% average over channels - per stage
cfgtemp = [];
cfgtemp.avgoverchan = 'yes';
cfgtemp.avgoverrpt = 'yes';
for ipatient = 1:8
    for ipart = 1 : 3
        for hyplabel = unique(FFT{ipatient}{ipart}.window.trialinfo.hyplabel)'
            cfgtemp.trials = FFT{ipatient}{ipart}.window.trialinfo.hyplabel == hyplabel;
            FFTavg_stage_part{ipatient}{ipart}.(hyplabel) = ft_selectdata(cfgtemp, FFT{ipatient}{ipart}.window);      
        end
    end
end

%% normalize by presleep
for ipatient = 1:8
    for hyplabel = unique(FFT{ipatient}{ipart}.window.trialinfo.hyplabel)'
        for ipart = 1 : 3
            FFTavg_stage_part_norm{ipatient}{ipart}.(hyplabel) = FFTavg_stage_part{ipatient}{ipart}.(hyplabel);
            FFTavg_stage_part_norm{ipatient}{ipart}.(hyplabel).powspctrm = FFTavg_stage_part{ipatient}{ipart}.(hyplabel).powspctrm ./ FFTavg_stage_part{ipatient}{ipart}.PRESLEEP.powspctrm;
        end
        FFTavg_stage_norm{ipatient}.(hyplabel) = FFTavg_stage_part_norm{ipatient}{1}.(hyplabel);
        FFTavg_stage_norm{ipatient}.(hyplabel).powspctrm = nanmedian([FFTavg_stage_part_norm{ipatient}{1}.(hyplabel).powspctrm; FFTavg_stage_part_norm{ipatient}{2}.(hyplabel).powspctrm; FFTavg_stage_part_norm{ipatient}{3}.(hyplabel).powspctrm]);
    end
end

figure;
cm = cool(size(fields(FFTavg_stage_norm{ipatient}), 1));
for ipatient = 1 : 8
    subplot(2,4,ipatient)
    hold;
    ic = 1;
    for hyplabel = string(fields(FFTavg_stage_norm{ipatient}))'
        plot(FFTavg_stage_norm{ipatient}.(hyplabel).freq, log(FFTavg_stage_norm{ipatient}.(hyplabel).powspctrm), 'color', cm(ic, :));
        ic = ic + 1;
    end
end

%% average over patients
for ipatient = 1 : 8
    for hyplabel = unique(FFT{ipatient}{ipart}.window.trialinfo.hyplabel)'
        FFT_all.(hyplabel)(ipatient, :) = log(FFTavg_stage_norm{ipatient}.(hyplabel).powspctrm);
    end
end

x = FFTavg_stage_norm{ipatient}.REM.freq;

FFT_all.S1      = FFT_all.PHASE_1;      FFT_all = rmfield(FFT_all, 'PHASE_1');
FFT_all.S2      = FFT_all.PHASE_2;      FFT_all = rmfield(FFT_all, 'PHASE_2');
FFT_all.S3      = FFT_all.PHASE_3;      FFT_all = rmfield(FFT_all, 'PHASE_3');
FFT_all.Wake    = FFT_all.AWAKE;        FFT_all = rmfield(FFT_all, 'AWAKE');
FFT_all.Post    = FFT_all.POSTSLEEP;    FFT_all = rmfield(FFT_all, 'POSTSLEEP');
FFT_all = rmfield(FFT_all, 'NO_SCORE');
FFT_all = rmfield(FFT_all, 'PRESLEEP');

% hyplabels = ["S3", "S2", "S1", "REM", "Post", "Wake"];
hyplabels = ["S3", "S2", "S1", "REM", "Wake"];

%% plot

fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
set(fig, 'position', [10 10 1000/sqrt(2) 1000]);
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');
hold on;
set(gca, 'clipping', 'on');    
ylim([-0.6, 2.8]); 

cm = cbrewer('qual', 'Set2', size(fields(FFT_all), 1));
cm(strcmp(hyplabels, 'Wake'), :) = [0, 0, 0];

% plot variance
ic = 1;
for hyplabel = hyplabels
    m = mean(FFT_all.(hyplabel));
    sem = std(FFT_all.(hyplabel)) / sqrt(size(FFT_all.(hyplabel), 1));
    m = smooth(m)';
    sem = smooth(sem)';    
    [~, indx] = max(m);
    maxfreq.(hyplabel) = x(indx); 
    maxindx.(hyplabel) = indx; 
    maxpow.(hyplabel)  = m(indx); 
    if ~strcmp(hyplabel, 'Wake')
        p = patch([x, x(end:-1:1)], [m - sem, m(end:-1:1) + sem(end:-1:1)], cm(ic, :), 'facealpha', 0.5, 'edgecolor', 'none'); uistack(p, 'bottom');
    end
    ic = ic + 1;
end

% plot mean
ic = 1;
for hyplabel = hyplabels
    m = mean(FFT_all.(hyplabel));
    m = smooth(m);
    %     plot(x, m, 'color', cm(ic, :), 'LineWidth', 2);
    if ~strcmp(hyplabel, 'Wake')
        l = plot(x, m, 'color', cm(ic, :), 'LineWidth', 2); uistack(l, 'bottom');
    else
        l = plot(x, m, 'color', cm(ic, :), 'LineWidth', 1, 'LineStyle', '--'); uistack(l, 'bottom');
    end
    ic = ic + 1;
end

ax = axis;

% add zero line
% l = plot([ax(1), ax(2)], [0, 0], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');

% plot max freq lines
ic = 1;
for hyplabel = ["S3", "S2", "S1", "REM"]
    l = plot([maxfreq.(hyplabel), maxfreq.(hyplabel)], [ax(3),maxpow.(hyplabel)], 'color', [0.5, 0.5, 0.5], 'LineStyle', '--'); uistack(l, 'bottom');
    p = plot(maxfreq.(hyplabel), maxpow.(hyplabel), 'o-','MarkerFaceColor', cm(ic, :),'MarkerEdgeColor', cm(ic, :));
    ic = ic + 1;
end


% add patches for frequency bands
% l = plot([2.5, 2.5],    [ax(3), ax(4)], 'color', [0.5, 0.5, 0.5], 'LineStyle', '-'); uistack(l, 'bottom');
% l = plot([4, 4],        [ax(3), ax(4)], 'color', [0.5, 0.5, 0.5]); uistack(l, 'bottom');
cshade = [0.8, 0.8, 0.8];
fill([0.1, 2.45, 2.45, 0.1], [2.6, 2.6, 3, 3], 'r', 'FaceColor', cshade, 'EdgeColor', cshade);
fill([2.55, 4, 4, 2.55],     [2.6, 2.6, 3, 3], 'b', 'FaceColor', cshade, 'EdgeColor', cshade);

text(1.30,  2.7, 'Slow Wave Activity', 'HorizontalAlignment', 'center');
text(3.25,  2.7, 'Delta', 'HorizontalAlignment', 'center');

xticks(sort([0.1, 2.5, 4, 5, maxfreq.REM, maxfreq.S3]));
xtickformat('%.1f')
xlim([0.1, 4]);
yticks([-1, 0, 1, 2]);

% plot for legend
ic = 1;
clear h
for hyplabel = hyplabels
    h(ic) = patch([0, 0], [0,0], cm(ic, :), 'facealpha', 1, 'edgecolor', 'none');
    ic = ic + 1;
end

xlabel('Frequency (Hz)');
ylabel('Power relative to pre-sleep');
legend(h, hyplabels, 'box', 'off', 'location', 'eastoutside'); 

% xlim([1, 20]);
set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.01, 0.01])
set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 12);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Tahoma');
axis square

%% write to figure for article (if not on Desktop PC)
if ~ispc
    fname = fullfile(config{1}.imagesavedir, 'FFT');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
else
    fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'FFT');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
end

