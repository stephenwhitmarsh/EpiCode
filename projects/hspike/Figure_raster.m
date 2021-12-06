function Figure_raster

restoredefaultpath
if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/git/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip    
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
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
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epishare-master'));
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\SPIKY_apr_2021'));
    addpath          \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\MatlabImportExport_v6.0.0
end

config                                                          = hspike_setparams;
MuseStruct{ipatient}                                            = readMuseMarkers(config{ipatient}, false);
MuseStruct{ipatient}                                            = padHypnogram(MuseStruct{ipatient});
MuseStruct{ipatient}                                            = alignMuseMarkersXcorr(config{ipatient}, MuseStruct{ipatient}, false);
[clusterindx{ipatient}, LFP_cluster{ipatient}]                  = clusterLFP(config{ipatient}, MuseStruct{ipatient}, false);
[config{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}] = alignClusters(config{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6});
[MuseStruct{ipatient}, ~, LFP_cluster_detected{ipatient}]       = detectTemplate(config{ipatient}, MuseStruct{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
[config{ipatient}, MuseStruct{ipatient}]                        = addSlidingWindows(config{ipatient}, MuseStruct{ipatient});

% read spike data
SpikeRaw{ipatient}                                              = readSpikeRaw_Phy(config{ipatient}, true);

% epoch to IEDs and sliding windows
config{ipatient}.spike.name                                     = ["template1", "template2", "template3", "template4", "template5", "template6", "window"];
SpikeTrials{ipatient}                                           = readSpikeTrials(config{ipatient}, MuseStruct{ipatient}, SpikeRaw{ipatient}, true);
SpikeStats{ipatient}                                            = spikeTrialStats(config{ipatient}, SpikeTrials{ipatient}, true);

config{ipatient}.LFP.name                                       = ["template1", "template2", "template3", "template4", "template5", "template6"];
LFP{ipatient}                                                   = readLFP(config{ipatient}, MuseStruct{ipatient}, false);

% 
% 
% 
% for inodule = nodulenr
% 
% %     SpikeTrials{inodule}  = readSpikeTrials(cfg{inodule});
%     LFP{inodule}          = readLFP(cfg{inodule});
% end

% 
% for iunit = 1 : size(SpikeTrials{4}{ipart}.SEIZURE.trial, 2)
%     i2 = SpikeTrials{4}{ipart}.SEIZURE.trial{iunit} == 2;
%     i1 = SpikeTrials{4}{ipart}.SEIZURE.trial{iunit} == 1;
%     SpikeTrials{4}{ipart}.SEIZURE.trial{iunit}(i2) = 1;
%     SpikeTrials{4}{ipart}.SEIZURE.trial{iunit}(i1) = 2;
% end
%     
%     
% % remove artefacted and non-nodular seizure
% if isunix
%     addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/modified_fieldtrip_functions
% else
%     addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development\modified_fieldtrip_functions
% end
% cfg_temp = [];
% cfg_temp.trials = [1, 2, 4, 5, 6, 7];
% SpikeTrials{4}{ipart}.SEIZURE = ft_spike_select_rmfulltrials(cfg_temp, SpikeTrials{4}{ipart}.SEIZURE);


for ipatient = 1 : 8
    
    % rereference to bipolar
    if strcmp(config{ipatient}.cluster.reref, 'yes')
        
        for itemplate = 1 : 6
            
            if contains(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, '-')
                disp('Rereferencing already done');
                continue
            end
            disp('Rereferencing');
            labels_nonum    = regexprep(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, '[0-9_]', '');
            [~,~,indx]      = unique(labels_nonum);
            clear group
            for i = 1 : max(indx)
                cfgtemp             = [];
                cfgtemp.reref       = 'yes';
                cfgtemp.refmethod   = 'bipolar';
                cfgtemp.demean      = 'yes';
                cfgtemp.baselinewindow = [-0.3, -0.1];
                cfgtemp.channel     = LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label(indx==i);
                group{i}            = ft_preprocessing(cfgtemp, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
            end
            LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)) = ft_appenddata([], group{:});
        end
    end
end

% configure figure
fig = figure('visible', true);
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'position', get(0,'ScreenSize'));
set(fig, 'PaperOrientation', 'portrait');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer', 'Painters');

% parameters for subplots
ncols       = 6;
nrows       = 8 * 2;
w           = 0.90/ncols;
hratio      = 0.98/ncols;
vratio      = 1/nrows;
htop        = 1/nrows * 0.8;
hbottom     = 1/nrows * 0.8;
spacehor    = 0.3;
rshift      = 0;
upshift     = -0.01;
imarker     = 1;
xlimits     = [-0.2, 0.8];
rownr       = 0;

for ipatient = 1 : 1
    rownr = rownr + 1;

   
    for itemplate = 1 : 6
        
        cfg = [];
        cfg.avgoverrpt = 'yes';
        LFPavg{ipatient}{ipart} = ft_selectdata(cfg, LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)));
        maxchan = find(~cellfun(@isempty, strfind(LFP{ipatient}{ipart}.(sprintf('template%d', itemplate)).label, config{ipatient}.align.zerochannel)), 1, 'first');
        
        % LFP
        row     = (rownr-1)*2 + 1;
        s1      = axes('Position', [hratio*(itemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
        plot(LFPavg{ipatient}{ipart}.time{1}, LFPavg{ipatient}{ipart}.trial{1}(maxchan, :), 'k');
        set(s1, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);
        xlim(xlimits);
        xticks([xlimits(1), 0, xlimits(2)])
        y = ylim;
        ylim([-ceil(max(abs(y))), ceil(max(abs(y)))]);
        yticks([-ceil(max(abs(y))), 0, ceil(max(abs(y)))]);
        
        if itemplate == 1
            l1 = ylabel('Amplitude (μV)');
            set(l1, 'Units', 'normalized');
            lpos = l1.Position;
            set(l1, 'Position', [-0.1, 0.5, 0]);
        end
        
        
        % next row
        row     = (rownr-1)*2 + 2;
        
        s2      = axes('Position', [hratio*(itemplate-1) + w*(spacehor/2) + rshift, 1-(vratio*row) + upshift, w, htop]);
        set(s2, 'XGrid', 'off', 'box', 'off', 'xticklabel', [], 'XColor', 'none', 'tickdir', 'out', 'Color', 'none', 'ticklength', [0.005 0.100]);

        % raster
        cfg_raster              = [];
%         cfg_raster.spikechannel = 1;
        cfg_raster.trialborders = 'no';
        %     cmap = repmat([linspace(0, 0.8, size(SpikeTrials{ipatient}{ipart}.SEIZURE.label, 2))], 3, 1)';
        %     cfg_raster.cmapneurons  = cmap;
        cfg_raster.cmapneurons  = gray(size(SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)).label, 2));
        ft_spike_plot_raster(cfg_raster, SpikeTrials{ipatient}{ipart}.(sprintf('template%d',itemplate)));
        xlim(xlimits);
        xticks([xlimits(1), 0, xlimits(2)])
        y = ylim;
        ylim([1, ceil(max(abs(y)))]);
        yticks([1, ceil(max(abs(y)))]);
        set(gca,'ticklength', [0.005 0.100]);

        if itemplate == 1
            l1 = ylabel('Trial');
            set(l1, 'Units', 'normalized');
            lpos = l1.Position;
            set(l1, 'Position', [-0.1, 0.5, 0]);
        else
            ylabel([]);
        end

    end
    
end


% write to figure
fname = fullfile(cfg{ipart}.imagesavedir, 'article', 'Figure4');
exportgraphics(fig, strcat(fname, '.jpg'),  'Resolution', 300);
exportgraphics(fig, strcat(fname, '.pdf'));
