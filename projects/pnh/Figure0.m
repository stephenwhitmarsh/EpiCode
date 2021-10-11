function Figure0(cfg, MuseStruct)

cfg.plot.units          = ft_getopt(cfg.plot, 'units', []);
cfg.spike.name          = cfg.plot.name;
cfg.LFP.name            = cfg.plot.name;

% ncols = max(3, size(cfg.LFP.name, 2));
ncols   = 2;
nrows   = 6;
        
% get file and samplinfo to recover trials
SpikeTrials_timelocked  = readSpikeTrials(cfg);
SpikeDensity_timelocked = spikeTrialDensity(cfg);
SpikeWaveforms          = readSpikeWaveforms(cfg);
cfg.spike.postfix       = '-windowed';
SpikeStats_windowed     = spikeTrialStats(cfg);

% highpassed micro (of wires used for spike sorting)
directory   = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.directory(cfg.plot.ievent, :);
temp        = dir(fullfile(cfg.rawdir, directory, ['*', cfg.circus.channel{1}, '.ncs']));
dataset     = fullfile(cfg.rawdir, directory, temp.name);
hdr         = ft_read_header(dataset);

spikeFs     = SpikeDensity_timelocked{cfg.plot.ipart}.psth.(cfg.plot.name).cfg.previous.hdr.Fs;
ss          = round(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.t0(cfg.plot.ievent) / spikeFs * hdr.Fs);
es          = round(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.t1(cfg.plot.ievent) / spikeFs * hdr.Fs);

Startsample = ss + cfg.epoch.toi.(cfg.plot.name)(1) * hdr.Fs - cfg.epoch.pad.(cfg.plot.name) * hdr.Fs;
Endsample   = es + cfg.epoch.toi.(cfg.plot.name)(2) * hdr.Fs + cfg.epoch.pad.(cfg.plot.name) * hdr.Fs;
Offset      = (cfg.epoch.toi.(cfg.plot.name)(1) - cfg.epoch.pad.(cfg.plot.name)) * hdr.Fs;

for ifile = 1 : size(cfg.circus.channel, 2)
    
    temp = dir(fullfile(cfg.rawdir, directory, ['*', cfg.circus.channel{ifile}, '.ncs']));
    if isempty(temp)
        fprintf('Could not find %s\n', cfg.circus.channel{ifile});
        continue
    end
    
    cfgtemp             = [];
    cfgtemp.dataset     = fullfile(cfg.rawdir, directory, temp.name);
    cfgtemp.trl         = round([Startsample; Endsample; Offset]');
    cfgtemp.hpfilter    = 'yes';
    cfgtemp.hpfreq      = 300;
    filedat{ifile}      = ft_preprocessing(cfgtemp);
    filedat{ifile}.label{1} = ['μ', filedat{ifile}.label{1}(end)];
    
end

cfgtemp                = [];
cfgtemp.keepsampleinfo = 'no';
LFP_micro_hp           = ft_appenddata(cfgtemp, filedat{:});
clear filedat*

% lowpass (LFP)
directory   = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.directory(cfg.plot.ievent, :);
temp        = dir(fullfile(cfg.rawdir, directory, ['*', cfg.plot.micro{1}, '.ncs']));
dataset     = fullfile(cfg.rawdir, directory, temp.name);
hdr         = ft_read_header(dataset);

spikeFs     = SpikeDensity_timelocked{cfg.plot.ipart}.psth.(cfg.plot.name).cfg.previous.hdr.Fs;
ss          = round(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.t0(cfg.plot.ievent) / spikeFs * hdr.Fs);
es          = round(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo.t1(cfg.plot.ievent) / spikeFs * hdr.Fs);

Startsample = ss + cfg.epoch.toi.(cfg.plot.name)(1) * hdr.Fs - cfg.epoch.pad.(cfg.plot.name) * hdr.Fs;
Endsample   = es + cfg.epoch.toi.(cfg.plot.name)(2) * hdr.Fs + cfg.epoch.pad.(cfg.plot.name) * hdr.Fs;
Offset      = (cfg.epoch.toi.(cfg.plot.name)(1) - cfg.epoch.pad.(cfg.plot.name)) * hdr.Fs;

for ifile = 1 : size(cfg.plot.micro, 2)
    
    temp = dir(fullfile(cfg.rawdir, directory, ['*', cfg.plot.micro{ifile}, '.ncs']));
    if isempty(temp)
        fprintf('Could not find %s\n', cfg.plot.micro{ifile});
        continue
    end
    
    cfgtemp             = [];
    cfgtemp.dataset     = fullfile(cfg.rawdir, directory, temp.name);
    cfgtemp.trl         = round([Startsample; Endsample; Offset]');
    cfgtemp.lpfilter    = 'yes';
    cfgtemp.lpfreq      = 50;
    cfgtemp.demean      = 'yes';
    filedat{ifile}      = ft_preprocessing(cfgtemp);
    filedat{ifile}.label{1} = ['μ', filedat{ifile}.label{1}(end)];
end

cfgtemp                = [];
cfgtemp.keepsampleinfo = 'no';
LFP_micro_lp           = ft_appenddata(cfgtemp, filedat{:});
clear filedat*

% get macro electrodes
cfg.LFP.postfix   = "_macro";
cfg.LFP.channel   = cfg.plot.macro;
temp              = readLFP(cfg, MuseStruct, false);
LFP_macro         = temp{cfg.plot.ipart}.(cfg.plot.name);

cfgtemp           = [];
cfgtemp.channel   = cfg.plot.macro;
cfgtemp.trials    = cfg.plot.ievent;
LFP_macro         = ft_selectdata(cfgtemp, LFP_macro);

for i = 1 : size(LFP_macro.label, 1)
    LFP_macro.label{i} = ['M', LFP_macro.label{i}(end:end)];
end

% bipolar rereferencing
cfgtemp             = [];
cfgtemp.reref       = 'yes';
cfgtemp.refmethod   = 'bipolar';
LFP_macro           = ft_preprocessing(cfgtemp, LFP_macro);

% plot figure
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [200 300 1000/sqrt(2) 1000]);
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

% LFP macro
subaxis(nrows, ncols, 1, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
n = 1; ytick = []; label = [];
maxrange = max(max(abs(LFP_macro.trial{1}))) * 2;
for ichan = 1 : size(LFP_macro.label, 1)
    ytick = [ytick, n*maxrange];
    x       = LFP_macro.time{1};
    y       = LFP_macro.trial{1}(ichan, :);
    plot(x, y + n*maxrange, 'k');
    label{ichan} = LFP_macro.label{ichan};
    n = n + 1;
end
yticks(ytick);
set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', label, 'TickDir', 'out')
axis tight;
xlim(cfg.epoch.toi.(cfg.plot.name));
t1 = title("A", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [-0.1, t1.Position(2), t1.Position(3)]);
yl = ylabel("Macro");
set(yl, 'Units', 'normalized');
yl.Position(1) = -0.08;

% micro lowpassed
subaxis(nrows, ncols, 3, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
n = 1; ytick = []; label = [];
maxrange = max(max(abs(LFP_micro_lp.trial{1}))) * 1.5;
for ichan = 1 : size(LFP_micro_lp.label, 1)
    ytick = [ytick, n*maxrange];
    x       = LFP_micro_lp.time{1};
    y       = LFP_micro_lp.trial{1}(ichan, :);
    plot(x, y + n*maxrange, 'k');
    label{ichan} = LFP_micro_lp.label{ichan};
    n = n + 1;
end
yticks(ytick);
set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', label, 'TickDir', 'out')
axis tight;
xlim(cfg.epoch.toi.(cfg.plot.name));
t1 = title("B", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [-0.1, t1.Position(2), t1.Position(3)]);
yl = ylabel("Micro");
set(yl, 'Units', 'normalized');
yl.Position(1) = -0.08;

% micro highpassed
subaxis(nrows, ncols, 5, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
n = 1; ytick = []; label = [];
maxrange = max(max(abs(LFP_micro_hp.trial{1}))) * 1.5;
for ichan = 1 : size(LFP_micro_hp.label, 1)
    ytick = [ytick, n*maxrange];
    x       = LFP_micro_hp.time{1};
    y       = LFP_micro_hp.trial{1}(ichan, :);
    plot(x, y + n*maxrange, 'k');
    label{ichan} = LFP_micro_hp.label{ichan};
    n = n + 1;
end
yticks(ytick);
set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'yticklabel', label, 'TickDir', 'out')
axis tight;
xlim(cfg.epoch.toi.(cfg.plot.name));
t1 = title("C", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [-0.1, t1.Position(2), t1.Position(3)]);
yl = ylabel("Micro");
set(yl, 'Units', 'normalized');
yl.Position(1) = -0.08;

% color units
if isempty(cfg.plot.units)
    cfg.plot.units = 1 : length(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).label);
end
    
cm = turbo(length(cfg.plot.units));
ci = 1;
for iunit = 1 : length(cfg.plot.units)
    spikeidx        = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trial{cfg.plot.units(iunit)} == cfg.plot.ievent;
    spiketime       = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).time{cfg.plot.units(iunit)}(spikeidx);
    [~, chanindx]   = max(rms(permute(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).template{cfg.plot.units(iunit)}(1,:,:), [2, 3, 1]), 2));
    for ispike = 1 : size(spiketime, 2)
        t1 = spiketime(ispike) - 0.001;
        t2 = spiketime(ispike) + 0.001;
        sel = LFP_micro_hp.time{1} >= t1 & LFP_micro_hp.time{1} <= t2;
        if any(sel)
            plot(LFP_micro_hp.time{1}(sel), LFP_micro_hp.trial{1}(chanindx, sel) + chanindx * maxrange, 'color', cm(ci, :));
        end
    end
    ci = ci + 1;
end

% spikes in single line
y = 1;
ci = 1;
clear label
subaxis(nrows, ncols, 7, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
for iunit = 1 : length(cfg.plot.units)
    
    spikeidx        = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trial{cfg.plot.units(iunit)} == cfg.plot.ievent;
    spiketime       = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).time{cfg.plot.units(iunit)}(spikeidx);
    if ~isempty(spiketime)
        for ispike = 1 : size(spiketime, 2)
            t1 = spiketime(ispike) - 0.001;
            t2 = spiketime(ispike) + 0.001;
            sel = LFP_micro_hp.time{1} >= t1 & LFP_micro_hp.time{1} <= t2;
            if any(sel)
                x = LFP_micro_hp.time{1}(sel);
                for i = 1 : size(x, 2)
                    plot([x(i); x(i)], [y-0.5, y+0.5] , 'color', cm(ci, :));
                end
            else
                plot([0, 0], [y-0.5, y+0.5], 'color', [1, 1, 1]);
            end
        end
    else
        plot([0, 0], [y-0.5, y+0.5], 'color', [1, 1, 1]);
    end
    y = y + 1;
    ci = ci + 1;
end
axis tight
xlim(cfg.epoch.toi.(cfg.plot.name));
yticks([1 : length(cfg.plot.units)])
set(gca,'TickLabelInterpreter', 'none', 'XColor', 'none', 'box', 'off', 'ydir', 'reverse')

t1 = title("E", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [-0.1, 1, t1.Position(3)]);
yl = ylabel("Unit");
set(yl, 'Units', 'normalized');
yl.Position(1) = -0.08;

% raster

subaxis(nrows, ncols, [9, 11], 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;

if size(SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name).trialinfo, 1) >= 200
    cfgtemp = [];
    cfgtemp.trials = 100:200;
    sel = ft_spike_select_rmfulltrials(cfgtemp, SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name));
else
    sel = SpikeTrials_timelocked{cfg.plot.ipart}.(cfg.plot.name);
end
cfg_raster              = [];
cfg_raster.trialborders = 'no';
cfg_raster.cmapneurons  = cm;
cfg_raster.spikechannel = cfg.plot.units;
ft_spike_plot_raster(cfg_raster, sel);
xlim(cfg.epoch.toi.(cfg.plot.name));
xticks(cfg.epoch.toi.(cfg.plot.name));
set(gca, 'XGrid', 'off', 'box', 'off', 'TickDir', 'out');
xlabel('time (s)');
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',8)
yl = ylabel("Trials"); set(yl, 'Units', 'normalized'); yl.Position(1) = -0.08;
yt = yticks; yticks([1, yt(end)]);
t1 = title("F", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [-0.1, t1.Position(2), t1.Position(3)]);

% waveforms
n = length(cfg.plot.units);
nx = ceil(sqrt(n));
ny = ceil(n/nx);
s = subaxis(nrows, ncols, [6, 8], 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
set(gca, 'color', 'none', 'XColor','none', 'YColor', 'none');

h = s.Position(4) / ny;
w = s.Position(3) / nx;

t1 = title("D", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [0.05, 1, t1.Position(3)]);

x = s.Position(1);
y = s.Position(2);
for itemp = 1 : n
    
    y_all   = vertcat(SpikeWaveforms{cfg.plot.ipart}{itemp}.trial{:})';
    r       = rms(y_all - mean(y_all, 2));
    [~, i]  = sort(r, 'ascend');
    y_sel   = y_all(:, i(1:100));
    
    col = mod(itemp-1, nx);
    row = ny - ceil(itemp / nx);
    if (row == 0) && (col == 0)
        s2{itemp} = axes('Position', [col*w+x+0.05, row*h+y, w, h], 'color', 'none'); hold on
        xticks([0, 100]);
        xl = xlabel('time (ms)'); 
        set(xl, 'Units', 'normalized');
        xl.Position(2) = -0.05;
    else
        s2{itemp} = axes('Position', [col*w+x+0.05, row*h+y, w, h], 'color', 'none', 'XColor', 'none', 'YColor', 'none', 'TickDir', 'out'); hold on
    end

    for i = 1 : size(y_sel, 2)
        lh = plot(y_sel(:, i), 'k'); lh.Color = [cm(itemp,:), 0.2];
    end
    
    plot(mean(y_all, 2), 'k', 'linewidth', 2);
    plot(mean(y_all, 2), 'w', 'linewidth', 1);
    if (row == 0) && (col == 0)
%         yt = yticks;
%         yticks([0, yt(end)]);
        yticks([]);
        ylabel('Amplitude (μV)');
    end
end

% ISI histograms
s = subaxis(nrows, ncols, [10, 12], 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
set(gca, 'color', 'none', 'XColor','none', 'YColor', 'none');

t1 = title("G", 'Fontsize', 14);
set(t1, 'Units', 'normalized');
set(t1, 'Position', [0.05, 1, t1.Position(3)]);

x = s.Position(1);
y = s.Position(2);
for itemp = 1 : n
    
    col = mod(itemp-1, nx);
    row = ny - ceil(itemp / nx);
    
    if (row == 0) && (col == 0)
        s2{itemp} = axes('Position', [col*w+x+0.05, row*h+y, w, h], 'color', 'none'); hold on
        xticks([0, 25]);
        xl = xlabel('time (ms)');
        set(xl, 'Units', 'normalized');
        xl.Position(2) = -0.05;
    else
        s2{itemp} = axes('Position', [col*w+x+0.05, row*h+y, w, h], 'color', 'none', 'XColor','none', 'YColor', 'none', 'TickDir', 'out'); hold on
    end
    
    bh = bar(SpikeStats_windowed{cfg.plot.ipart}.window{itemp}.isi_avg_time*1000, SpikeStats_windowed{cfg.plot.ipart}.window{itemp}.isi_avg, 1, 'FaceColor',cm(itemp,:), 'EdgeColor', 'none');
    xlim([0, 25]);
    if (row == 0) && (col == 0)
%         yt = yticks;
%         yticks([0, xt(end)]);
        ylabel('Count');
        yticks([]);
        set(gca,'Layer','top')
    end  
end
% set(findall(fig, '-property', 'fontsize'), 'fontsize', 8);

% write to figure
fname = fullfile(cfg.imagesavedir, 'article', [cfg.prefix, 'FigureS0_trial', num2str(cfg.plot.ievent)]);
% exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
% exportgraphics(fig, strcat(fname, '.pdf'));
savefig(fig, fname);
disp("Done");