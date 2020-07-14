function rpv = plot_spike_quality(cfg,spikedata,spikewaveforms,force)
%plot 12*6, max pour 36 units (2 plots par unit)

%get default values
cfg.spike.RPV               = ft_getopt(cfg.spike, 'RPV', 0.003);
cfg.spikequal               = ft_getopt(cfg,'spikequal',[]);
cfg.spikequal.plot_std      = ft_getopt(cfg,'spikequal','yes');
cfg.spikequal.suffix        = ft_getopt(cfg,'spikequal',[]);
cfg.spikequal.part_list     = ft_getopt(cfg.spikequal,'part_list','all');
cfg.spikequal.label_list    = ft_getopt(cfg.spikequal,'label_list','all');

% load precomputed stats if required
fname = fullfile(cfg.datasavedir,[cfg.prefix,'spike_rpv',cfg.spikequal.suffix,'.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed spike rpv\n');
    load(fname,'rpv');
    return
end


if strcmp(cfg.spikequal.part_list, 'all')
    cfg.spikequal.part_list = 1:size(spikedata);
end

for ipart = cfg.spikequal.part_list
    
    if strcmp(cfg.spikequal.label_list, 'all')
        cfg.spikequal.label_list = 1:size(spikedata{ipart},2);
    end
    
    for ilabel = cfg.spikequal.label_list
        
        fig=figure;
        i=1;
        
        if size(spikedata{ipart}{ilabel}.label,2) > 36
            error('More then 36 units in part %d label %d. Adapt the nr of subplots', ipart, ilabel);
        end
        
        %compute isi
        cfgtemp         = [];
        cfgtemp.bins    = [0:0.001:0.025];
        isi             =  ft_spike_isi(cfgtemp,spikedata{ipart}{ilabel});
        
        for i_unit=1:size(spikedata{ipart}{ilabel}.label,2)
            if contains(spikedata{ipart}{ilabel}.cluster_group{i_unit}, 'noise')
                rpv{ipart}{ilabel}(i_unit) = NaN;
                continue
            end
            
            % compute rpv
            rpv{ipart}{ilabel}(i_unit) = ( sum(isi.isi{i_unit} < cfg.spike.RPV) / length(isi.isi{i_unit}) ) * 100;
                        
            %plot isi
            subplot(12,6,i)

            %color for mua
            if contains(spikedata{ipart}{ilabel}.cluster_group{i_unit}, 'mua')
                c = 'k';
            else
                %color for sua
                c = 'b';
            end
            
            bar(isi.time,isi.avg(i_unit,:), 'FaceColor', c, 'EdgeColor', c);
            xticklabels(xticks*1000); %convert in ms
            ax = axis;
            titlepos = title(sprintf('%s (E%d), RPV : %.2f%%, %d spikes',spikedata{ipart}{ilabel}.label{i_unit}, spikedata{ipart}{ilabel}.template_maxchan(i_unit), rpv, length(isi.isi{i_unit})),'Interpreter', 'none','Fontsize', 4);
            titlepos.Position(1) = ax(1);
            titlepos.HorizontalAlignment = 'left';
            setfig();
            
            subplot(12,6,i+1)
            cfgtemp                            = [];
            cfgtemp.morpho.mesurepeaktrough    = 'yes';
            cfgtemp.morpho.mesuretroughpeak    = 'no';
            cfgtemp.morpho.mesurehalfwidth     = 'yes';
            cfgtemp.morpho.blmethod            = 'min';
            cfgtemp.morpho.channame            = spikewaveforms{ipart}{ilabel}{i_unit}.label{1};
            cfgtemp.morpho.plotstd             = cfg.spikequal.plot_std;
            cfgtemp.morpho.plotraw             = 'no';
            cfgtemp.morpho.name                = cfg.name{ilabel};
            plot_morpho(cfgtemp,spikewaveforms{ipart}{ilabel}{i_unit});
            delete(findall(gca,'type','text'));
            
            xticklabels(xticks*1000); %convert in ms
            title([]);
            setfig();
            
            i=i+2; %nr of the subplot
            
        end
        %% print to file
        fprintf('Print to file and save data\n');
        
        if ~isfolder(fullfile(cfg.imagesavedir))
            mkdir(fullfile(cfg.imagesavedir));
        end
        
        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_',cfg.name{ilabel},'_spike_quality',cfg.spikequal.suffix,'.pdf']),'-r600');
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'_',cfg.name{ilabel},'_spike_quality',cfg.spikequal.suffix,'.png']),'-r600');
        
        close all
        
    end
end

save(fname, 'rpv');

end

% setfig = @() set(gca,'FontWeight','bold','TickDir','out');
function setfig()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set common features to all the subplots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
set(gca,'Fontsize',4);
end


