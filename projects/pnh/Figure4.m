function Figure4

% settings
cfg         = pnh_setparams;
nodulenr    = [1, 4];
seizurenr   = [1, nan, nan, 2];
% seizurenr   = [1, nan, nan, 1];
channelnr   = [1, nan, nan, 1];
ipart       = 1;
xlimits     = [-1, 5];



% load data
x = load("\\lexport\iss01.charpier\analyses\vn_pnh\data\pnh\2230-SpikeTrials_Timelocked.mat");
SpikeTrials{1}  = x.SpikeTrials;

x = load("\\lexport\iss01.charpier\analyses\vn_pnh\data\pnh\2689-SpikeTrials_Timelocked.mat");
SpikeTrials{4}  = x.SpikeTrials;

for inodule = nodulenr

%     SpikeTrials{inodule}  = readSpikeTrials(cfg{inodule});
    LFP{inodule}          = readLFP(cfg{inodule});
end


for iunit = 1 : size(SpikeTrials{4}{ipart}.SEIZURE.trial, 2)
    i2 = SpikeTrials{4}{ipart}.SEIZURE.trial{iunit} == 2;
    i1 = SpikeTrials{4}{ipart}.SEIZURE.trial{iunit} == 1;
    SpikeTrials{4}{ipart}.SEIZURE.trial{iunit}(i2) = 1;
    SpikeTrials{4}{ipart}.SEIZURE.trial{iunit}(i1) = 2;
end
    
    
% remove artefacted and non-nodular seizure
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/modified_fieldtrip_functions
else
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development\modified_fieldtrip_functions
end
cfg_temp = [];
cfg_temp.trials = [1, 2, 4, 5, 6, 7];
SpikeTrials{4}{ipart}.SEIZURE = ft_spike_select_rmfulltrials(cfg_temp, SpikeTrials{4}{ipart}.SEIZURE);

% configure figure
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

% parameters for subplots
nrows       = 12;
w           = 0.92;
hratio      = 1/3;
vratio      = 1/nrows;
htop        = 1/nrows * 0.8;
hbottom     = 1/nrows * 0.8;
spacehor    = 0.1;
rshift      = 0.02;
upshift     = -0.01;
imarker     = 1;

% loop over nodules
rownr   = 0;
for ipatient = nodulenr
    
    % LFP
    rownr   = rownr + 1;
    row     = (rownr-1)*2 + 1;
    s1      = axes('Position', [hratio*(imarker-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
    plot(LFP{ipatient}{ipart}.SEIZURE.time{seizurenr(ipatient)}, LFP{ipatient}{ipart}.SEIZURE.trial{seizurenr(ipatient)}(channelnr(ipatient), :), 'k');
    set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
    xlim(xlimits);
    xticks([xlimits(1), 0, xlimits(2)]);

    y = ylim;
    yticks([ceil(y(1)), 0, floor(y(end))]);
    l1 = ylabel('Amplitude (μV)');
    set(l1, 'Units', 'normalized');
    lpos = l1.Position;
    %     t1 = title(sprintf('Nodule %d', ipatient), 'Fontsize', 14);
    if ipatient == 1
        t1 = title("A", 'Fontsize', 18);
    else
        t1 = title("B", 'Fontsize', 18);
    end
    set(t1, 'Units', 'normalized');
    set(t1, 'Position', [-0.02, t1.Position(2), t1.Position(3)]);
    clear t1
    
    % next row
    row     = (rownr-1)*2 + 2;
    s2      = axes('Position', [hratio*(imarker-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + (vratio-hbottom) + upshift, w, hbottom], 'Color', 'none', 'XColor', 'none');
    set(s2, 'XGrid', 'off', 'box', 'off', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
 
    % raster
    cfg_raster              = [];
    cfg_raster.trialborders = 'no';
    cmap = repmat([linspace(0, 0.8, size(SpikeTrials{ipatient}{ipart}.SEIZURE.label, 2))], 3, 1)';
    cfg_raster.cmapneurons  = cmap;
%     cfg_raster.cmapneurons  = gray(size(SpikeTrials{ipatient}{ipart}.SEIZURE.label, 2));
    ft_spike_plot_raster(cfg_raster, SpikeTrials{ipatient}{ipart}.SEIZURE);
    set(gca,'Ydir', 'reverse', 'ticklength', [0.005 0.100]);
 
    y = ylim;
    yticks([y(2)]); ylabel('');
    %     xlim(cfg{1}.epoch.toi.SEIZURE);
    %     xticks([cfg{ipatient}.epoch.toi.SEIZURE(1), 0, cfg{ipatient}.epoch.toi.SEIZURE(2)]);
    xlim(xlimits);
    xticks([xlimits(1), 0, xlimits(2)]);
    yticks([]);
    l2 = ylabel('Seizures');
    set(l2, 'Units', 'normalized');   
    set(l2, 'Position', lpos);
  
    l3 = xlabel("Time from seizure onset (seconds)");
    set(l3, 'Position', [l3.Position(1), 7, l3.Position(3)]);
    
end


% write to figure
fname = fullfile(cfg{ipart}.imagesavedir, 'article', 'Figure4');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));
