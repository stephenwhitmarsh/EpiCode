
function [spike_LFP] = spikeLFP(cfg, MuseStruct, SpikeRaw, force)

fname = fullfile(cfg.datasavedir,[cfg.prefix,'spike_LFP.mat']);

if exist(fname,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed spike-LFP data ***\n');
    fprintf('************************************\n\n');
    
    load(fname,'spike_LFP');
    
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing spike-LFP data ***\n');
    fprintf('********************************\n\n');
    

    for itemp = 1 : size(SpikeRaw.label,2)
        
        temp                = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.LFP.channel{ifile}, '.ncs']));
        fname{1}            = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
        dat                 = ft_read_neuralynx_interp(fname);
        
        
        hdr = ft_read_header(MuseStruct.directorylist{});
        
        trl = [];
        trl(:,1) = SpikeRaw.sample{itemp} - (cfg.spike.pre * hdr.Fs);
        trl(:,2) = SpikeRaw.sample{itemp} + (cfg.spike.post * hdr.Fs);
        trl(:,3) = ones(length(SpikeRaw.sample{itemp}),1) * -(cfg.spike.pre * hdr.Fs);
                
        cfgtemp                 = [];
%         cfgtemp.hpfilter        = 'yes';
%         cfgtemp.hpfreq          = 500;
        cfgtemp.demean          = 'yes';
        cfgtemp.baselinewindow  = cfg.spike.baseline; 
        cfgtemp.dataset         = fnames{ichan};
        cfgtemp.trl             = trl;
        spike_LFP{itemp}        = ft_preprocessing(cfgtemp);
    end
    
    save(fname,'spike_LFP','-v7.3');
    
    
end
% 
% figure;
% ft_singleplotER([],spike_LFP{itemp})
% 
% spike_LFP_trl{itemp}    = nan(size(spike_LFP{itemp}.trial,2),size(spike_LFP{itemp}.time{1},2));
% for trialnr = 1 : size(spike_LFP{itemp}.trial,2)
%     spike_LFP_trl{itemp}(trialnr,:) = spike_LFP{itemp}.trial{trialnr};
% end
%     
%     % remove the most extreme to avoid scaling problem
% %     spike_LFP_avg{ichan}    = spike_LFP_avg(rms(spike_LFP_avg,2) < (std(rms(spike_LFP_avg,2))*2),:);
%     
% %     subplot(size(cfg.circus.channel,2),2,(ichan*2)-1);
% %     hold;
%     
% figure; hold;
% rng('shuffle')
% trialnr = 1;
% while trialnr < 1000
%     r = randi(size(SpikeRaw.sample{itemp},1),1,1);
%     if size(find(SpikeRaw.sample{itemp} < SpikeRaw.sample{itemp}(r)+cfg.spike.width & SpikeRaw.sample{itemp} > SpikeRaw.sample{itemp}(r)-cfg.spike.width),1) == 1
%         plot(spike_LFP{itemp}.time{r},spike_LFP{itemp}.trial{r},'r');
%         trialnr = trialnr + 1;
%     end
% end
% plot(spike_LFP{itemp}.time{1},median(spike_LFP_trl{itemp}),'k','linewidth',2);
% 
% spike_LFP_std = std(spike_LFP_trl{itemp});
% spike_LFP_avg = median(spike_LFP_trl{itemp}); 
% 
% figure; hold;
% plot(spike_LFP{itemp}.time{1},spike_LFP_avg,'k-');
% plot(spike_LFP{itemp}.time{1},spike_LFP_avg+spike_LFP_std,'k:');
% plot(spike_LFP{itemp}.time{1},spike_LFP_avg-spike_LFP_std,'k:');
% set(gca,'XTick',(round(spike_LFP{itemp}.time{1}(1),4):0.001:round(spike_LFP{itemp}.time{1}(end),4)));
% 
% format long
% 
% 
% 
% 
%     %
%     %         subplot(size(cfg.circus.channel,2),2,(ichan*2)); hold;
%     %
%     %         if ichan == SpikeRaw.template_maxchan(itemp)
%     %
%     %             plot(spike_LFP{ichan}.time{1},median(spike_LFP_avg{ichan}),'r','linewidth',2);
%     %         else
%     %             plot(spike_LFP{ichan}.time{1},median(spike_LFP_avg{ichan}),'k','linewidth',2);
%     %         end
%     %     end
%     %
%     %     % print to file
%     %     set(fig,'PaperOrientation','landscape');
%     %     set(fig,'PaperUnits','normalized');
%     %     set(fig,'PaperPosition', [0 0 1 1]);
%     %     print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'spiketriggered_average_ofspike_',num2str(itemp),'.pdf']),'-r600');
%     %     set(fig,'PaperOrientation','portrait');
%     %     print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'spiketriggered_average_ofspike_',num2str(itemp),'.png']),'-r600');
%     %     close all
