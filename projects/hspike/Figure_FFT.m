function [config] = Figure_FFT

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/cbrewer/cbrewer      
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
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\cbrewer\cbrewer
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end
ft_defaults

%% load data
for ipatient = 1:8
    config                      = hspike_setparams;
    config{ipatient}.FFT.name   = {'window'};
    FFT{ipatient}               = FFTtrials(config{ipatient}, false);
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

hyplabels = ["S3", "S2", "S1", "REM", "Post", "Wake"];

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


cm = cbrewer('qual', 'Set2', size(fields(FFT_all), 1));

% plot variance
ic = 1;
for hyplabel = hyplabels
    m = mean(FFT_all.(hyplabel));
    sem = std(FFT_all.(hyplabel)) / sqrt(size(FFT_all.(hyplabel), 1));
    patch([x, x(end:-1:1)], [m - sem, m(end:-1:1) + sem(end:-1:1)], cm(ic, :), 'facealpha', 0.5, 'edgecolor', 'none');
    ic = ic + 1;
end

% plot mean
ic = 1;
for hyplabel = hyplabels
    m = mean(FFT_all.(hyplabel));
    sem = std(FFT_all.(hyplabel)) / sqrt(size(FFT_all.(hyplabel), 1));
    plot(x, m, 'color', cm(ic, :), 'LineWidth', 2);
    ic = ic + 1;
end

% plot for legend
ic = 1;
for hyplabel = hyplabels
    h(ic) = patch([0, 0], [0,0], cm(ic, :), 'facealpha', 1, 'edgecolor', 'none');
    ic = ic + 1;
end

ax = axis;

% add zero line
l = plot([ax(1), ax(2)], [0, 0], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');

% add patches for frequency bands

l = plot([2.5, 2.5],    [ax(3), ax(4)], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');
l = plot([4, 4],        [ax(3), ax(4)], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');
text(1.25,  2.2, 'Delta1', 'HorizontalAlignment', 'center');
text(3.25,  2.2, 'Delta2', 'HorizontalAlignment', 'center');
xticks([0.1, 2.5, 4, 5]);
xlim([0.1, 5]);

% l = plot([4, 4],   [ax(3), ax(4)], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');
% l = plot([7, 7],   [ax(3), ax(4)], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');
% l = plot([14, 14], [ax(3), ax(4)], 'color', [0.9, 0.9, 0.9]); uistack(l, 'bottom');
% text(2.5,  2.2, 'Delta', 'HorizontalAlignment', 'center');
% text(5.5,  2.2, 'Theta', 'HorizontalAlignment', 'center');
% text(10.5, 2.2, 'Alpha', 'HorizontalAlignment', 'center');
% xticks([1, 4, 7, 14, 20]);

yticks([-1, 0, 1, 2]);

xlabel('Frequency (Hz)');
ylabel('Log power relative to pre-sleep');
legend(h, hyplabels, 'box', 'off'); 

% xlim([1, 20]);
set(gca,'TickLabelInterpreter', 'none', 'box', 'off', 'TickDir', 'out', 'TickLength', [0.01, 0.01])
set(findall(gcf, '-property', 'FontSize'), 'Fontsize', 18);

% write to figure for article (if not on Desktop PC)
if ~ispc
    fname = fullfile(config{1}.imagesavedir, 'FFT');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
else
    fname = fullfile('D:\Dropbox\Apps\Overleaf\Hspike\images', 'FFT');
    isdir_or_mkdir(fileparts(fname));
    exportgraphics(fig, strcat(fname, '.pdf'));
end

