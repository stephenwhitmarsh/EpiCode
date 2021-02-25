function spiketriggeredplots(cfg,dat_microFs,SpikeTrials,SpikeRaw)

%% read data on high-Fs

for ichan = 1 : size(cfg.fnames_ncs,2)
    cfgtemp = [];
    cfgtemp.dataset = cfg.fnames_ncs{ichan};
    dat_microFs_all{ichan} = ft_preprocessing(cfgtemp); 
end

dat_microFs_all_concatinated = ft_appenddata([],dat_microFs_all{:});
template_labels = ft_channelselection('temp*', dat_microFs_all_concatinated, 'all');

clear dat_microFs_all

%% Spike-triggered LFP

% combine spike and LFP
temp                        = dir(fullfile(cfg.datasavedir,[num2str(ipatient),'-all_data_',cfg.channel{imarker}(1:end-2),'_*.ncs']));
hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
dat_microFs_all_concatinated.hdr             = hdr;

data_all = ft_appendspike([],dat_microFs_all_concatinated,SpikeTrials{imarker});

template_labels = ft_channelselection('temp*', data_all, 'all');

for itemp = 1 : size(template_labels,1)
    
    %% interpolate LFP around spikes
    cfgtemp                         = [];
    cfgtemp.method                  = 'nan';
    cfgtemp.timwin                  = [-0.002 0.002];
    cfgtemp.spikechannel            = template_labels{itemp};
    cfgtemp.channel                 = cellstr(cfgtemp.micro_labels{ipatient,imarker});
    
    data_i                          = ft_spiketriggeredinterpolation(cfgtemp,data_all);
    
    cfgtemp                         = [];
    cfgtemp.timwin                  = [-1 1];
    cfgtemp.spikechannel            = template_labels{itemp};
    cfgtemp.channel                 = cellstr(cfgtemp.micro_labels{ipatient,imarker});
    cfgtemp.latency                 = [-cfg.prestim{ipatient}(imarker) cfg.poststim{ipatient}(imarker)];
    spikeavg                        = ft_spiketriggeredaverage(cfgtemp, data_i);
    
    fig = figure;
    plot(spikeavg.time,spikeavg.avg);
    legend
    title(['Spike-triggered average ',cfg.label{imarker}]);
    
    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'-spiketriggered_average_',cfg.label{imarker},'_template',num2str(itemp),'.pdf']),'-r600');
    set(fig,'PaperOrientation','portrait');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'-spiketriggered_average_',cfg.label{imarker},'_template',num2str(itemp),'.png']),'-r600');
end % itemplate


%% combining LFP with spikes

% remove trials that were somehow lost, i.e. fell off edge of file
[~, Ai, Bi] = intersect(dat_microFs.trialinfo(:,7),SpikeTrials{imarker}.trialinfo(:,7));

cfgtemp                     = [];
cfgtemp.trials              = Ai;
dat_microFs                 = ft_selectdata(cfgtemp,dat_microFs);

cfgtemp                     = [];
cfgtemp.trials              = Bi;
SpikeTrials{imarker}        = ft_spike_select(cfgtemp,SpikeTrials{imarker});

% add header - better to do it in readMicroFs
temp                        = dir(fullfile(cfg.datasavedir,[num2str(ipatient),'-all_data_',cfg.channel{imarker}(1:end-2),'_*.ncs']));
hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
dat_microFs.hdr             = hdr;

% recreate trl
dat_microFs.cfg.trl         = SpikeTrials{1}.trialinfo(Bi,3:4); % should be able to remove Bi after spike_select will also remove trials in trialinfo
dat_microFs.cfg.trl(:,3)    = -ones(size(dat_microFs.cfg.trl,1),1) * cfg.prestim(imarker) * hdr.Fs;

% combine spike and LFP
data_all = ft_appendspike([],dat_microFs,SpikeTrials{imarker});

%         data_all_sel = rmfield(data_all,'sampleinfo');
%                 % plot spikes together with lfp with databrowser
%                 cfg                 = [];
%                 cfg.viewmode        = 'vertical';
%                 cfg.continuous      = 'yes';
%                 cfg.mychan          = ft_channelselection('temp*', data_all, 'all');
%                 cfg.mychanscale     = repmat(max(max(data_all.trial{1}))*2,size(cfg.mychan,1));
%                 cfg.channelcolormap = [0,0,0;1,0,0];
%                 cfg.colorgroups     = [ones(1,size(data_all.label,1)-size(cfg.mychan,1)) ones(1,size(cfg.mychan,1))*2];
%                 cfg.linewidth       = 1;
%                 ft_databrowser(cfg,data_all_sel);

%% plot trials with LFP and spikes - no filter

for itrial = 1:10
    h = 300;
    
    fig = figure; hold;
    i1 = find(dat_microFs.time{itrial} >= -cfg.prestim(imarker),1,'first');
    i2 = find(dat_microFs.time{itrial} >=  cfg.poststim(imarker),1,'first');
    
    for ichannel = 1 : size(dat_microFs.label,1)
        plot(dat_microFs.time{itrial}(i1:i2),dat_microFs.trial{itrial}(ichannel,i1:i2) + ichannel*h,'color','k');
    end
    
    cmap = lines(size(SpikeTrials{imarker}.time,2));
    
    for itemp = 1:size(SpikeTrials{imarker}.time,2)
        
        maxchan     = SpikeTrials{1}.template_maxchan(itemp);
        indx        = find(SpikeTrials{imarker}.trial{itemp} == itrial);
        spiketime   = SpikeTrials{imarker}.time{itemp}(indx);
        if ~isempty(spiketime)
            for si = 1 : size(spiketime,2)
                
                xi  = find(dat_microFs.time{itrial} >= spiketime(si),1,'first');
                y   = dat_microFs.trial{itrial}(maxchan,xi) + maxchan*h;
                
                if spiketime(si) > -cfg.prestim(imarker) && spiketime(si) < cfg.poststim(imarker)
                    plot(spiketime(si),y+10,'v','markerfacecolor',cmap(itemp,:),'markeredgecolor',[1,1,1]);
                end
                
            end
        end
    end
    
    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'-LFPwithSpikes_',cfg.label{imarker},'_trial',cfg.num2str(itrial),'.pdf']),'-r300');
    set(fig,'PaperOrientation','portrait');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'-LFPwithSpikes_',cfg.label{imarker},'_trial',cfg.num2str(itrial),'.png']),'-r300');
end


%% plot trials with LFP and spikes - highpass filter
cfgtemp = [];
cfgtemp.hpfilter = 'yes';
cfgtemp.hpfreq = 500;
dat_microFs_hp = ft_preprocessing(cfgtemp,dat_microFs);

for itrial = 1:10
    h = 100;
    
    fig = figure; hold;
    i1 = find(dat_microFs_hp.time{itrial} >= -cfg.prestim{ipatient}(imarker),1,'first');
    i2 = find(dat_microFs_hp.time{itrial} >=  cfg.poststim{ipatient}(imarker),1,'first');
    
    for ichannel = 1 : size(dat_microFs_hp.label,1)
        plot(dat_microFs_hp.time{itrial}(i1:i2),dat_microFs_hp.trial{itrial}(ichannel,i1:i2) + ichannel*h,'color','k');
    end
    
    cmap = lines(size(SpikeTrials{imarker}.time,2));
    
    for itemp = 1:size(SpikeTrials{imarker}.time,2)
        
        maxchan     = SpikeTrials{1}.template_maxchan(itemp);
        indx        = find(SpikeTrials{imarker}.trial{itemp} == itrial);
        spiketime   = SpikeTrials{imarker}.time{itemp}(indx);
        if ~isempty(spiketime)
            for si = 1 : size(spiketime,2)
                
                xi  = find(dat_microFs_hp.time{itrial} >= spiketime(si),1,'first');
                y   = dat_microFs_hp.trial{itrial}(maxchan,xi) + maxchan*h;
                
                if spiketime(si) > -cfg.prestim{ipatient}(imarker) && spiketime(si) < cfg.poststim{ipatient}(imarker)
                    plot(spiketime(si),y+10,'v','markerfacecolor',cmap(itemp,:),'markeredgecolor',[1,1,1]);
                end
                
            end
        end
    end
    
    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'-LFPwithSpikes_HP_',cfg.label{imarker},'_trial',num2str(itrial),'.pdf']),'-r300');
    set(fig,'PaperOrientation','portrait');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'-LFPwithSpikes_HP_',cfg.label{imarker},'_trial',num2str(itrial),'.png']),'-r300');
end



% 
% %% Spike triggered spectrum
% cfg = [];
% cfg.method = 'mtmconvol';
% cfg.foi = 20:10:100;
% cfg.t_ftimwin = 5./cfg.foi;
% cfg.taper = 'hanning';
% %     stsConvol = ft_spiketriggeredspectrum(cfg,data_all); % better to use third argument
% stsConvol = ft_spiketriggeredspectrum(cfg,dat_microFs,SpikeTrials{imarker}); % better to use third argument
% 
% for k = 1:length(stsConvol.label)
%     
%     % compute the statistics on the phases
%     cfg               = [];
%     cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
%     %         excludeChan       = str2num(stsConvol.label{k}(6)); % exclude the same channel
%     %    chan              = true(1,length(stsConvol.label));
%     %         chan(k)             = false;
%     cfg.spikechannel  = 2;  % stsConvol.label{k};
%     cfg.channel       = 1;    %stsConvol.lfplabel(chan); bug% selected LFP channels
%     cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
%     cfg.timwin        = 'all'; % compute over all available spikes in the window
%     %         cfg.latency       = [0.3 nanmax(stsConvol.trialtime(:))]; % sustained visual stimulation period
%     statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
%     
%     % plot the results
%     figure
%     plot(statSts.freq,statSts.ppc0')
%     xlabel('frequency')
%     ylabel('PPC')
% end
% 
% 
% % compute the statistics on the phases
% cfg               = [];
% cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
% %         excludeChan       = str2num(stsConvol.label{k}(6)); % exclude the same channel
% %    chan              = true(1,length(stsConvol.label));
% %         chan(k)             = false;
% cfg.spikechannel  = 5;  % stsConvol.label{k};
% cfg.channel       = 2;    %stsConvol.lfplabel(chan); bug% selected LFP channels
% cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
% cfg.timwin        = 'all'; % compute over all available spikes in the window
% %         cfg.latency       = [0.3 nanmax(stsConvol.trialtime(:))]; % sustained visual stimulation period
% cfg.avgoverchan = 'unweighted';
% cfg.winstepsize = 0.01;
% cfg.timwin = 0.5;
% cfg.latency = [-prestim{ipatient}(imarker) poststim{ipatient}(imarker)];
% %                 cfg.latency = [0 poststim{ipatient}(imarker)];
% 
% statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);
% 
% %         % plot the results
% %         figure
% %         plot(statSts.freq,statSts.ppc0')
% %         xlabel('frequency')
% %         ylabel('PPC')
% 
% 
% figure
% cfg            = [];
% cfg.parameter  = 'ppc0';
% cfg.refchannel = statSts.labelcmb{1,1};
% cfg.channel = statSts.labelcmb{1,2};
% ft_singleplotTFR(cfg, statSts)

