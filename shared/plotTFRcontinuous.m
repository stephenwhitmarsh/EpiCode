function plotTFRcontinuous(cfg,TFR)

for ipart = 1 : size(TFR,2)
    fig = figure;
    
    for ichan = 1 : size(TFR{ipart}.label,1)
        
        subplot(size(TFR{ipart}.label,1),1,ichan);
        
        cfgtemp                 = [];
        cfgtemp.channel         = TFR{ipart}.label{ichan};
        %         cfgtemp.ylim            = [1 40];
        cfgtemp.baseline        = [TFR{1}.time(1)+120 TFR{1}.time(end)-120];
        cfgtemp.baselinetype    = 'relative';
        cfgtemp.colorbar        = 'no';
        cfgtemp.colorbar        = 'yes';
        %     cfgtemp.zlim            = 'maxabs';
        %     cfgtemp.xlim            = config{ipatient}.;
        %     cfgtemp.title           = 'Relative change from Baseline';
        cfgtemp.parameter       = 'powspctrm';
        cfgtemp.colormap        = parula(5000);
        cfgtemp.renderer        = 'painters';
        ft_singleplotTFR(cfgtemp,TFR{ipart});
    
    end
    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-TFR_',TFR{ipart}.label{ichan},'.pdf']),'-r600');
    close all
end
