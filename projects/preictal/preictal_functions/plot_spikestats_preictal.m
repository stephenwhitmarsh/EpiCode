function [stats] = plot_spikestats_preictal(cfg, stats, spikedata)
% field added to stats for outputs : isi, waveform measurments, descriptive
% stats of the firing for all the windows of stats

% Mandatory input :

% Optional inputs
cfg.spikewaveforms          = ft_getopt(cfg, 'spikewaveforms', []);
cfg.statstime               = ft_getopt(cfg, 'statstime', []);
cfg.statstime.minbadtime    = ft_getopt(cfg.statstime, 'minbadtime', 0);
cfg.statstime.plot          = ft_getopt(cfg.statstime, 'plot', []);
cfg.statstime.plot.toi      = ft_getopt(cfg.statstime.plot,'toi',[-Inf Inf]);
cfg.statstime.plot.suffix   = ft_getopt(cfg.statstime.plot,'suffix',[]);
cfg.statstime.plot.marker   = ft_getopt(cfg.statstime.plot,'marker',[]);
cfg.statstime.plot.bad_start= ft_getopt(cfg.statstime.plot,'bad_start',[]);
cfg.statstime.plot.bad_end  = ft_getopt(cfg.statstime.plot,'bad_end',[]);
cfg.spike                   = ft_getopt(cfg, 'spike', []);
cfg.spike.ISIbins           = ft_getopt(cfg.spike, 'ISIbins', [0:0.003:0.150]);
cfg.spike.RPV               = ft_getopt(cfg.spike, 'RPV', 0.003);
cfg.circus                  = ft_getopt(cfg, 'circus', []);
cfg.circus.channel          = ft_getopt(cfg.circus, 'channel', cell(1,max(spikedata{1}{1}.template_maxchan)+1)); % name of the chans used by spyking circus (to set the same of the title). Default = empty cells


for ipart = 1:size(spikedata,2)
    for ilabel = 1:size(stats{ipart}, 2)
        for i_unit = 1:size(spikedata{ipart}.label, 2)
            
            %% prepare figure
            fig = figure;
            phy_channr      = spikedata{ipart}.template_maxchan(i_unit);
            nlx_channame    = cfg.circus.channel{phy_channr+1}; %+1 because phy_channr starts at zero
            sgtitle(sprintf('Electrode n°%d %s : %s %s',phy_channr,nlx_channame, spikedata{ipart}.label{i_unit}, cfg.statstime.plot.suffix), 'Fontsize', 22, 'Interpreter', 'none', 'FontWeight', 'bold');
            
            %% text with average values of the whole data
            subplot(4,4,1);hold;
            RPV             = ( sum(stats{ipart}{ilabel}.isi{i_unit} < cfg.spike.RPV) / length(stats{ipart}{ilabel}.isi{i_unit}) ) * 100;
            meanfreq        = nanmean(cell2mat(stats{ipart}{ilabel}.freq{i_unit}));
            meancv2         = nanmean(cell2mat(stats{ipart}{ilabel}.cv2{i_unit}));
            meanamplitude   = nanmean(cell2mat(stats{ipart}{ilabel}.amplitude{i_unit}));
            meanff          = nanmean(cell2mat(stats{ipart}{ilabel}.fanofactor{i_unit}));
            meanbi          = nanmean(cell2mat(stats{ipart}{ilabel}.burstindex{i_unit}));
            meancv          = nanmean(cell2mat(stats{ipart}{ilabel}.cv{i_unit}));
            
            disp_text = sprintf('Mean values : \n- Amplitude = %g µV \n- RPV = %g %% \n- Freq = %g Hz \n- CV2 = %g \n- Fano Factor = %g \n- CV = %g \n- Burst index = %g', meanamplitude, RPV, meanfreq, meancv2, meanff, meancv, meanbi);
            
            text(0,1,disp_text,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
            axis off
            
            stats{ipart}{ilabel}.meanvalues.rpv(i_unit)         = RPV;
            stats{ipart}{ilabel}.meanvalues.freq(i_unit)        = meanfreq;
            stats{ipart}{ilabel}.meanvalues.cv2(i_unit)         = meancv2;
            stats{ipart}{ilabel}.meanvalues.amplitude(i_unit)   = meanamplitude;
            stats{ipart}{ilabel}.meanvalues.fanofactor(i_unit)  = meanff;
            stats{ipart}{ilabel}.meanvalues.burstindex(i_unit)  = meanbi;
            stats{ipart}{ilabel}.meanvalues.cv(i_unit)          = meancv;
            
            %% plot ISI
            subplot(4,4,2);hold;
            stats{ipart}{ilabel}.isi{i_unit} = diff(double(spikedata{ipart}.samples{i_unit})) / spikedata{ipart}{ilabel}.hdr.Fs ;
            histogram(stats{ipart}{ilabel}.isi{i_unit},cfg.spike.ISIbins,'EdgeColor','none');
            
            yticklabels(yticks);
            xticklabels(xticks*1000); %convert in ms
            ylabel('Spike count');
            title('ISI (ms)');
            xlim([0 Inf]);
            set(gca,'FontWeight','bold','TickDir','out');
            
            
            %% plot template
            subplot(4,4,3);hold;
            tempsel = spikedata{ipart}.template{i_unit}(:,spikedata{ipart}.template_maxchan(i_unit)+1,:);%+1 because electrodes nr are zero-based
            temptime = ( (0 : size(spikedata{ipart}.template{i_unit},3) - 1) / spikedata{ipart}{ilabel}.hdr.Fs )';
            %tempsel dimensions : 1:ntemplates 2:maxchan 3:values
            
            % interpolate template
            temptime_int = linspace(temptime(1),temptime(end),10000);
            tempsel_int  = pchip(temptime,tempsel,temptime_int);
            
            %convert to Fieldtrip 'raw' type
            template = [];
            for itrial=1:size(tempsel,1)
                template.time{itrial}   = temptime_int;
                template.trial{itrial}  = tempsel_int(itrial,:,:);
                %permute to remove one of the 2 units dimensions in case of several templates
                if size(template.trial{itrial}, 2) == 1
                    template.trial{itrial} = permute(template.trial{itrial},[1 3 2]);
                end
            end
            template.label{1} = 'template';
            
            cfgtemp                     = [];
            cfgtemp.morpho.channame            = 'template';
            cfgtemp.morpho.mesurehalfwidth     = 'yes';
            cfgtemp.morpho.halfwidthmethod     = 'min';
            cfgtemp.morpho.mesurepeaktrough    = 'yes';
            cfgtemp.morpho.toiac               = 'all';
            cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.halfwidthmethod = 'min';
            [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,template);
            
            title([]);
            ylabel('Template (µV)');
            xlabel('Time (ms)');
            xticklabels(xticks*1000);
            set(gca,'FontWeight','bold','TickDir','out');
            
            stats{ipart}{ilabel}.template.halfwidth{i_unit}    = halfwidth;
            stats{ipart}{ilabel}.template.peaktrough{i_unit}   = peaktrough;
            stats{ipart}{ilabel}.template.troughpeak{i_unit}   = troughpeak;
            stats{ipart}{ilabel}.template.time{i_unit}         = template.time;
            stats{ipart}{ilabel}.template.values{i_unit}       = template.trial;
            
            %% plot waveform
            subplot(4,4,4);hold;
            if ~isempty(cfg.spikewaveforms)
                if ~isempty(cfg.spikewaveforms{ipart}{1}{i_unit}) %always 1 because raw data not cut into trials, so no 'labels'
                    cfgtemp                     = [];
                    cfgtemp.morpho.channame            = cfg.SpikeWaveforms{ipart}{iplot}{i_unit}.label{1};
                    cfgtemp.morpho.plotstd             = 'yes';
                    cfgtemp.morpho.removeoutliers      = 'yes'; %if big noise, impair seeing real data. Still present in avg and std.
                    cfgtemp.morpho.mesurehalfwidth     = 'yes';
                    cfgtemp.morpho.halfwidthmethod     = 'min';
                    cfgtemp.morpho.mesurepeaktrough    = 'yes';
                    cfgtemp.morpho.toibl               = []; %no need of bl if cfgtemp.halfwidthmethod     = 'min';
                    cfgtemp.morpho.toiac               = 'all';
                    cfgtemp.morpho.name                = cfg.name{iplot};
                    [halfwidth, peaktrough, troughpeak] = plot_morpho(cfgtemp,cfg.spikewaveforms{ipart}{iplot}{i_unit});
                    
                    xlabel('Time (ms)');
                    xticklabels(xticks*1000); %convert in ms
                    ylabel('Spike waveform (µV)');
                    title([]);
                    set(gca,'FontWeight','bold','TickDir','out');
                    
                    stats{ipart}{ilabel}.spikewaveform.halfwidth{i_unit}      = halfwidth;
                    stats{ipart}{ilabel}.spikewaveform.peaktrough{i_unit}     = peaktrough;
                    stats{ipart}{ilabel}.spikewaveform.troughpeak{i_unit}     = troughpeak;
                else
                    axis off;
                end
            else
                axis off;
            end
            
            %% plot values over time
            i=5;
            for ivalue = ["freq","amplitude","cv", "cv2", "fanofactor", "burstindex"]
                subplot(4,4,[i i+1]);hold;
                
                %plot values
                for itrial = 1:size(spikedata{ipart}{ilabel}.trialinfo, 1)
                    scatter(stats{ipart}{ilabel}.time{i_unit}{itrial}, stats{ipart}{ilabel}.(ivalue){i_unit}{itrial},'.b');
                end
                
                %add BAD patch
                xlim(cfg.statstime.plot.toi);
                ax = axis;
                for iBAD = 1:size(cfg.statstime.plot.bad_start,2)
                    %cfg.minbadtime
                    if cfg.statstime.plot.bad_end(iBAD) - cfg.statstime.plot.bad_start(iBAD) < cfg.statstime.minbadtime, continue, end
                    x = [cfg.statstime.plot.bad_start(iBAD) cfg.statstime.plot.bad_end(iBAD) cfg.statstime.plot.bad_end(iBAD) cfg.statstime.plot.bad_start(iBAD)];
                    y = [ax(3) ax(3) ax(4) ax(4)];
                    patch('XData',x,'YData',y,'facecolor',[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
                end
                
                %plot again values to overdraw bad patchs
                for itrial = 1:size(spikedata{ipart}{ilabel}.trialinfo, 1)
                    scatter(stats{ipart}{ilabel}.time{i_unit}{itrial}, stats{ipart}{ilabel}.(ivalue){i_unit}{itrial},'.b');
                end
                
                %add marker position if any
                for imarker = 1:size(cfg.statstime.plot.marker,2)
                    plot([cfg.statstime.plot.marker(imarker) cfg.statstime.plot.marker(imarker)], [ax(3) ax(4)], 'r');
                end
                
                ylabel(ivalue)
                xlim(cfg.statstime.plot.toi);
                set(gca,'FontWeight','bold','TickDir','out');
                
                if i == 17 || i == 19
                    xlabel('Time (s)');
                end
                
                i = i+2;
            end
            
            %% save figure
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',spikedata{ipart}.label{i_unit},'_spikestatsOverTime_raw_',strrep(num2str(i_toi'),'  ','_'),'_',cfg.statstime.plot.suffix,'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'-',spikedata{ipart}.label{i_unit},'_spikestatsOverTime_raw_',strrep(num2str(i_toi'),'  ','_'),'_',cfg.statstime.plot.suffix,'.png']),'-r600');
            close all
            
        end %i_unit
    end %ilabel
end %ipart

end