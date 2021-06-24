function Figure2(cfg_fig)

if nargin == 0
    cfg_fig = [];
end

cfg_fig.colorscale          = ft_getopt(cfg_fig, 'colorscale', 0.2);

n = 100;
cmap              = parula(n);
cmap              = [repmat(cmap(1, :), round(n*cfg_fig.colorscale), 1); cmap; repmat(cmap(n,:), round(n*cfg_fig.colorscale, 1), 1)];

cfg = pnh_setparams;

cfg_fig.plot.name = {'PSW', 'FA', 'ES'};
cfg_fig.plot.channel{1}.PSW     = 'm1pNs_6';
cfg_fig.plot.channel{1}.FA      = 'm1pNs_4';
cfg_fig.plot.channel{1}.ES      = 'm1pNs_4';
cfg_fig.plot.channel{2}.PSW     = 'mCasd_2';
cfg_fig.plot.channel{2}.FA      = 'mCasd_2';
cfg_fig.plot.channel{2}.ES      = 'mCasd_2';
cfg_fig.plot.channel{3}.PSW     = 'mTNmi_3';
cfg_fig.plot.channel{3}.FA      = 'mTNmi_3';
cfg_fig.plot.channel{3}.ES      = 'mTNmi_3';
cfg_fig.plot.channel{4}.PSW     = 'mLMI1_7';
cfg_fig.plot.channel{4}.FA      = 'mLMI1_7';
cfg_fig.plot.channel{4}.ES      = 'mLMI1_7';

cfg_fig.FFT.toi.PSW{1}          = [0 0.35];
cfg_fig.FFT.toi.FA{1}           = [0 0.6];
cfg_fig.FFT.toi.ES{1}           = [-0.1 0.1];
cfg_fig.FFT.toi.PSW{2}          = [-0.05 0.25];
cfg_fig.FFT.toi.FA{2}           = [0 0.8];
cfg_fig.FFT.toi.ES{2}           = [-0.1 0.1];
cfg_fig.FFT.toi.FA{3}           = [0 0.4];
cfg_fig.FFT.toi.PSW{3}          = [0 0.4];
cfg_fig.FFT.toi.ES{3}           = [-0.1 0.1];
cfg_fig.FFT.toi.PSW{4}          = [0 0.6];
cfg_fig.FFT.toi.FA{4}           = [0 0.6];
cfg_fig.FFT.toi.ES{4}           = [-0.1 0.1];            
                
cfg_TFR                         = [];
cfg_TFR.channel                 = 1;
cfg_TFR.colorbar                = 'no';
cfg_TFR.zlim                    = 'maxabs';
cfg_TFR.title                   = ' ';
cfg_TFR.baselinetype            = 'relchange';
cfg_TFR.interactive             = 'no';
cfg_TFR.ylim                    = [60, 200];

ipart = 1;

for ipatient = 1 : 4
    LFP = readLFP(cfg{ipatient});
    TFR{ipatient} = TFRtrials(cfg{ipatient});
    
    for markername = string(fields(LFP{ipart}))'
        if strcmp(markername, 'SEIZURE')
            continue
        end
    	cfg_temp = [];
        cfg_temp.channel = cfg_fig.plot.channel{ipatient}.(markername);
        LFPavg{ipatient}.(markername) = ft_timelockanalysis(cfg_temp, LFP{ipart}.(markername));
        TFR{ipatient}{1}.(markername) = ft_selectdata(cfg_temp, TFR{ipatient}{ipart}.(markername));
    end
    clear LFP
end


% to move PFA of patient 3 under PSW instead of FA
TFR{3}{1}.PSW           = TFR{3}{1}.FA;
TFR{3}{1}               = rmfield(TFR{3}{1}, 'FA');
LFPavg{3}.PSW           = LFPavg{3}.FA;
LFPavg{3}               = rmfield(LFPavg{3}, 'FA');
cfg{3}.TFR.bl.PSW       = cfg{3}.TFR.bl.FA;
cfg{3}.epoch.toi.PSW    = cfg{3}.epoch.toi.FA;


% configure figure
papersize = 1000;
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'position', [100 100  papersize papersize*sqrt(2)]);
set(fig, 'position', [100 100  papersize papersize*sqrt(2)]);

% set(fig, 'PaperOrientation', 'portrait');
% set(fig, 'PaperUnits', 'normalized');
% set(fig, 'PaperPosition', [0 0 1 1]);
% set(fig, 'Renderer', 'Painters');
% orient(fig,'portrait')

ncols   = 3;
nrows   = 4;

% first print average LFP and TFR
for ipatient = 1 : 4
    
    maxrange = 0;
    for markername = string(cfg_fig.plot.name)
        if isfield(LFPavg{ipatient}, markername)
            newmax = max(abs(LFPavg{ipatient}.(markername).avg)) * 1.15;
            maxrange = max([newmax, maxrange]);
        end
    end
    
    imarker = 1;
    for markername = string(cfg_fig.plot.name)
        
        if ~isfield(LFPavg{ipatient}, markername) || ~isfield(TFR{ipatient}{ipart}, markername)
            imarker = imarker + 1;
            continue
        end
        
        % TFR
        s = subaxis(nrows, ncols, imarker + (ipatient-1) * ncols, 'SpacingVert', 0.02, 'SpacingHoriz', 0.01); hold;
        cfg_TFR.figure      = s;
        cfg_TFR.baseline    = cfg{ipatient}.TFR.bl.(markername);
        cfg_TFR.xlim        = cfg{ipatient}.epoch.toi.(markername);
        cfg_TFR.colormap    = cmap;        
        ft_singleplotTFR(cfg_TFR, TFR{ipatient}{ipart}.(markername));
        
        set(gca, 'XGrid', 'off', 'box', 'off',  'TickDir', 'out');
        xticks([]);         
        yticks(cfg_TFR.ylim);
        
        if imarker == 1 || (imarker == 2 && ipatient == 3)
            l = ylabel('Frequency (Hz)');
            set(l, 'Units', 'normalized');
            x = l.Position;
            x(1) = -0.025;            
            set(l, 'Position', x);
            set(gca, 'TickDir', 'out')            
        else
            set(gca, 'yticklabel', [], 'TickDir', 'out')
            yticks([]);
        end
        
        if ipatient ~= 4
            set(gca, 'xticklabel', []);
        end
        
        z = caxis; 
        z_new = [ceil(z(1)*100), floor(z(2)*100)];
        
        c = colorbar('TickDirection','out','TickLength', 0.03); 
        set(c, 'Location', 'west', 'color', [0 0 0], 'Xtick', [z(1), 0, z(2)], 'Xticklabel', {num2str(z_new(1)), 0, num2str(z_new(2))});
        x = get(c, 'Position'); x(1) = x(1) + 0.01; x(2) = x(2) + 0.1; x(4) = x(4) / 3 + 0.01; set(c, 'Position', x); 

        l = ylabel(c, 'Change (%)', 'Rotation', 90);
        x = get(l, 'Position'); x(1) = -0.8; set(l, 'Position', x);

        % baseline line
        y = ylim;
        line(cfg{ipatient}.TFR.bl.(markername), [y(1), y(1)] + 5, 'linewidth', 2, 'color', 'k');
        text(cfg{ipatient}.TFR.bl.(markername)(1) + 0.05, y(1) + 7, 'Baseline', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        
        % time line
        y = ylim;
        if strcmp(markername, "ES")
            duration = 0.5;
            t = sprintf('%.1f second', duration);            
        else
            duration = 1;
            t = sprintf('1 second'); 
        end
        line([cfg{ipatient}.epoch.toi.(markername)(end)-duration, cfg{ipatient}.epoch.toi.(markername)(end)], [y(1), y(1)] + 5, 'linewidth', 2, 'color', 'k');
        text(cfg{ipatient}.epoch.toi.(markername)(end) - 0.05, y(1) + 7, t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        
        xlim(cfg{ipatient}.epoch.toi.(markername));
        switch markername
            case "ES"
                title("Epileptic Spike", 'FontSize', 14);
                
            case "PSW"
                title("Periodic Slow Waves", 'FontSize', 14);
                if ipatient == 3
                    title("Periodic Fast Activity", 'FontSize', 14);
                end
            case "FA"
                title("Fast Activity", 'FontSize', 14);
        end
        
        % FFT toi line
        y = ylim;
        line(cfg_fig.FFT.toi.(markername){ipatient}, [y(1), y(1)] + 5, 'linewidth', 2, 'color', 'k');
        text(mean(cfg_fig.FFT.toi.(markername){ipatient}), y(1) + 7, 'Inset', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
 
        %% LFP
        yyaxis right
        x = LFPavg{ipatient}.(markername).time;
        y = LFPavg{ipatient}.(markername).avg;
        plot(x, y, 'k');
        ylim([-maxrange, maxrange]);
        yticks([ceil(-(maxrange)), 0, floor((maxrange))]);
        
        xlim(cfg{ipatient}.epoch.toi.(markername));
        ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
        if imarker == 3
            set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'off', 'box', 'off', 'TickDir', 'out');
            l = ylabel('Amplitude (microVolts)');
            set(l, 'Units', 'normalized');
            x = l.Position;
            x(1) = 1.075;            
            set(l, 'Position', x);
        else
            set(gca,'TickLabelInterpreter', 'none', 'XGrid', 'off', 'box', 'off', 'TickDir', 'out', 'yticklabel', []);
            yticks([]);
        end
        if ipatient ~= 4
            set(gca, 'xticklabel', []);
        end
        
        
        %% FFT INSET
        x = s.Position;
        
        if (ipatient == 2 || ipatient == 4) && strcmp(markername,"FA")
            x(1) = x(1) + 0.03;     % left 
            x(2) = x(2) + 0.03;     % bottom
            x(3) = x(3) / 3;        % width
            x(4) = x(4) / 3.5;      % height
        else
            x(1) = x(1) + x(3) * 0.62;
            x(2) = x(2) + x(4)/2 + 0.03;
            x(3) = x(3) / 3;        % width
            x(4) = x(4) / 3.5;      % height
        end
        
        inset = axes('position', x, 'color', 'none'); hold on;
        
        if imarker == 1 && ipatient == 1
            xlabel('Frequency (Hz)');
            l = ylabel('Change (%)');
            set(l, 'Units', 'normalized');
            x = l.Position;
            x(1) = -0.025; 
            x(2) = 0.4;            
            set(l, 'Position', x);
        end
        
        cfg_bl              = [];
        cfg_bl.baseline     = cfg{ipatient}.TFR.bl.(markername);
        TFR_bl              = ft_freqbaseline(cfg_bl, TFR{ipatient}{ipart}.(markername));
        
        cfg_bl              = [];
        cfg_bl.frequency    = [60, 200];
        cfg_bl.latency      =  cfg_fig.FFT.toi.(markername){ipatient};
        TFR_bl              = ft_selectdata(cfg_bl, TFR_bl);
        
        x = TFR_bl.freq;
        y = TFR_bl.powspctrm;
        y = nanmean(y, 1);
        y = nanmean(y, 3);
        y = squeeze(y);
        y = smooth(y, 10);
        
        [Ypk, Xpk, ~, ~] = findpeaks(y, TFR_bl.freq, 'NPeaks', 2, 'SortStr', 'descend');
        
        if ipatient == 2 && strcmp(markername, "PSW")
            Ypk = Ypk(2);
            Xpk = Xpk(2);
        else
            Ypk = Ypk(1);
            Xpk = Xpk(1);
        end
        
        % vertical line for peak
        plot(x, y, 'k');
        yl = ylim;
        for i = 1 : length(Xpk)
            plot([Xpk(i), Xpk(i)], [yl(1), Ypk(i)], 'k:');
        end
        xlim([60, 200]);
        xticks(round(sort(unique([60, Xpk, 200]))));
        
        % vertical line for peak
        for i = 1 : length(Xpk)
            plot([0, Xpk(i)], [Ypk(1), Ypk(i)], 'k:');
        end
        yticks(Ypk);
        yticklabels(sprintf('%d',round(Ypk*100)));
        set(inset,'color','none', 'box', 'off', 'TickDir', 'out');
        imarker = imarker + 1;
    end  
end

fname = fullfile(cfg{ipart}.imagesavedir, 'article', ['Figure2','_scale_', num2str(cfg_fig.colorscale )]);
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 600);
exportgraphics(fig, strcat(fname, '.pdf'));


