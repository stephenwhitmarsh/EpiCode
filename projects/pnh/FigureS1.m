function FigureS1

config = pnh_setparams; 

macro{1} = {'_1pNs_1', '_1pNs_2', '_1pNs_3', '_1pNs_4', '_1pNs_5', '_1pNs_6', '_1pNs_7', '_1pNs_8'};
macro{2} = {'_Casd_1', '_Casd_2', '_Casd_3', '_Casd_4', '_Casd_5', '_Casd_6', '_Casd_7', '_Casd_8'};
macro{3} = {'_TNmi_1', '_TNmi_2', '_TNmi_3', '_TNmi_4', '_TNmi_5', '_TNmi_6', '_TNmi_7', '_TNmi_8'};
macro{4} = {'_LMI1_1', '_LMI1_2', '_LMI1_3', '_LMI1_4', '_LMI1_5', '_LMI1_6', '_LMI1_7', '_LMI1_8'};

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


ipart = 1;
for ipatient = 1 : 4
    
    % read muse markers
    MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, false);
    
    % align markers
    MuseStruct{ipatient} = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
        
    % read LFP (micro/default)
    LFP = readLFP(config{ipatient});

    % read LFP (macro)
    config{ipatient}.LFP.postfix = "macro";
    config{ipatient}.LFP.channel = macro{ipatient};    
    LFP_macro = readLFP(config{ipatient}, MuseStruct{ipatient}, false);
    
    for markername = string(cfg_fig.plot.name)
        try
            LFP_avg{ipatient}.(markername)          = ft_timelockanalysis([], LFP{ipart}.(markername));
            LFP_macro_avg{ipatient}.(markername)    = ft_timelockanalysis([], LFP_macro{ipart}.(markername));
        catch
        end
    end
    
    clear LFP_macro LFP
end

for ipatient = 1 : 4
    for markername = string(cfg_fig.plot.name)
        try
            cfg = [];
            cfg.reref = 'yes';
            cfg.refmethod = 'bipolar';
            LFP_macro_avg_bipolar{ipatient}.(markername) = ft_preprocessing(cfg, LFP_macro_avg{ipatient}.(markername) );
        catch
        end
    end
end

clear cor
for ipatient = 1 : 4
    for markername = string(fields(LFP_macro_avg{ipatient}))'
        try     
            for ichan = 1 : 7
                
                cfg         = [];
                cfg.channel = cfg_fig.plot.channel{ipatient}.(markername);
                temp        = ft_selectdata(cfg, LFP_avg{ipatient}.(markername));
                cor{ipatient}.(markername)(ichan) = corr(LFP_macro_avg_bipolar{ipatient}.(markername).avg(ichan, :)', temp.avg');
            end
        catch
        end
    end
end

% configure figure

fig = figure; hold;
set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');
xi = 0;
labels = {};
cm = parula(2);
for ipatient =  1 : 4
    for markername = string(fields(cor{ipatient}))'
        xi = xi + 1;
        
        for ichan = 1 : 7
            xpos = (xi-1) * 3;
            ypos = (ichan-1) * 3;
            if cor{ipatient}.(markername)(ichan) >= 0
                col = cm(2, :);
            else
                col = cm(1, :);
            end
            p1 = scatter(xpos, ypos, abs(cor{ipatient}.(markername)(ichan) * 3000), col, 'filled', 'MarkerFaceAlpha', 0.8);
            if ipatient == 3 && strcmp(markername,'FA')
                labels{xi} = sprintf('N%d-%s', ipatient, "PFA");                
            else
                labels{xi} = sprintf('N%d-%s', ipatient, markername);
            end
            text(xpos, ypos, sprintf('%.2f', cor{ipatient}.(markername)(ichan)), 'HorizontalAlignment', 'center', 'VerticalAlignment','middle');
        end
    end
end
xticks((0:xi)*3);
xticklabels(labels);
yticks((0:6)*3);
yticklabels([1:7]);
xlim([-3 xi*3]);
ylim([-3 7*3]);
xlabel('Nodule-Pattern');
ylabel('Macroelectrode');
title('Correlation micro with macroelectrodes','Fontsize', 18);
box off;

% write to figure
fname = fullfile(config{1}.imagesavedir, 'article', 'FigureS1');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));
