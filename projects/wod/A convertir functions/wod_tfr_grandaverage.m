function wod_grandaverage(cfg)


%load data
analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_recovery','timefreq_wor','timefreq_wor_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_recovery_blcorrected','timefreq_wor_blcorrected','timefreq_wor_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_recovery','log_timefreq_wor','log_timefreq_wor_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected', 'log_timefreq_recovery_blcorrected','log_timefreq_wor_blcorrected','log_timefreq_wor_timenorm_blcorrected', 'log_timefreq_baseline_blcorrected'};

%go through each data structure
for idata = 1:size( analysis_names,2)
    
    fprintf('For %s\n', analysis_names{idata});
    
    
    cfgcommon = cfg{4}; %same common config for all rats
    
    %initialize data structures
    data_allrats = [];
    count_trials                    = 0;
    
    % load timefreq data from all rats
    for irat = 1:size(cfg,2)
        %load data
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(cfg,2));
        data_temp = load(fullfile(cfg{irat}.datasavedir,[cfg{irat}.prefix,analysis_names{idata},'.mat']));
        
        %pool all rat data in the same structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});
            for ichan = 1:numel(chan_list)
                chan_name = chan_list{ichan};
                data_allrats.(chan_name){count_trials}     = data_temp.(analysis_names{idata}){itrial}.(chan_name);
            end
        end
        
        clear data_temp
    end %irat
    
    % pool data
    % if not the same number of points, keep only the number of points of
    % the shorter trial
    chan_list                                   = fieldnames(data_allrats);
    
    for ichan = 1:numel(chan_list)
        chan_name                               = chan_list{ichan};
        %replace empty channels by nan channels
        hasdata = find(~cellfun(@isempty,data_allrats.(chan_name)));
        if isempty(hasdata)
            continue
        end
        for iwod = 1:size(data_allrats.(chan_name),2)
            %replace empty channel by nan
            if isempty(data_allrats.(chan_name){iwod})
                data_allrats.(chan_name){iwod}                       = data_allrats.(chan_name){hasdata(1)};
                data_allrats.(chan_name){iwod}.powspctrm             = nan(size(data_allrats.(chan_name){iwod}.powspctrm));
            end
        end
        
        %data separated for each channel
        cfgtemp                 = [];
        cfgtemp.keepindividual  = 'yes';
        data.(chan_name)        = ft_freqgrandaverage(cfgtemp, data_allrats.(chan_name){:});
    end
    chan_list                                   = fieldnames(data);
    
    clear data_allrats
    
    
    
    
    
    %% TFR average by chan and over rats
    
    for ichan= 1:numel(chan_list)%chan_name = string(fieldnames(data)')
        
        chan_name                               = chan_list{ichan};
        %average over protocols
        data_plot = data.(chan_name);
        data_plot.label = {'channel'};
        data_plot.powspctrm = permute(nanmean(data.(chan_name).powspctrm,1), [2 3 4 1]);
        data_plot.dimord = 'chan_freq_time';
        
        fig=figure;
        
        sgtitle([analysis_names{idata},chan_name, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        subplot(2,1,1)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfgtemp.colormap= 'jet';
        cfg.fontsize = 12;
        cfgtemp.ylim= [51 100];
        cfgtemp.masknans    = 'yes';
        ft_singleplotTFR(cfgtemp, data_plot);
        ft_pimpplot(fig, jet(5000))
        
        
        title('Time-Frequency [50 : 100]');
        %figure settings (fontsize,fontweight, ticks dir, renderer etc.)
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time (s)');
        ylabel(sprintf('Frequency (Hz)'));
        axis tight;
        
        subplot(2,1,2)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfgtemp.masknans    = 'yes';
        cfgtemp.colormap= 'jet';
        cfgtemp.ylim= [0 50];
        cfg.fontsize = 12;
        ft_singleplotTFR(cfgtemp, data_plot);
        ft_pimpplot(fig, jet(5000))
        
        
        title('Time-Frequency [1 : 50]');
        %figure settings (fontsize,fontweight, ticks dir, renderer etc.)
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time (s)');
        ylim([1 50]);
        ylabel(sprintf('Frequency (Hz)'));
        axis tight;
        
        if ~isfolder(cfg{4}.imagesavedir)
            fprintf('Creating directory %s\n', cfg{4}.imagesavedir);
            mkdir(cfg{4}.imagesavedir)
        end
        
        
        if idata >= 7
            cfg{4}.imagesavedir_data{1} = fullfile(cfg{4}.imagesavedir_data{1}, 'baseline_corrected');
        end
        
        if idata > 12 & idata < 19
            cfg{4}.imagesavedir_data{1} = fullfile(cfg{4}.imagesavedir_data{3});
        end
        
        if idata >= 19
            cfg{4}.imagesavedir_data{1} = fullfile(cfg{4}.imagesavedir_data{3},'baseline_corrected');
        end
        
        if ~isfolder(cfg{4}.imagesavedir_data{1})
            fprintf('Creating directory %s\n', cfg{4}.imagesavedir);
            mkdir(cfg{4}.imagesavedir_data{1})
        end
        
        
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s.pdf', analysis_names{idata}, chan_name)),'-r600');
        print(fig, '-dpng', fullfile(cfg{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s.png', analysis_names{idata}, chan_name)),'-r600');
        savefig(fig,fullfile(cfg{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s', analysis_names{idata}, chan_name)));
        close all
        
        cfg{4}.imagesavedir_data{1}= fullfile(cfg{4}.imagesavedir,'TFR');
        
        clear data_plot
        
        
        
        
        
    end %ichan
    
    %% plot frequency bands over time for all rats
    
    fig = figure;
    sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
    
    for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
        
        subplot(2,2,ifreq);hold;
        C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
        
        for chan_name = string(fieldnames(data))'
            %go trhough each electrode in the ascending order
            chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
            legend_name{chan_nr} = chan_name;
            %                 if isempty(data{ifreq}.(chan_name))
            %                     continue
            %                 end
            color   = C(chan_nr,:);
            
            x = data.(chan_name).time;
            clear pow_band y
            
            %find frequencies of the band (comparer des valeurs : > <
            %>= <= ==)
            freq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
            pow_band= data.(chan_name).powspctrm(:,:,freq_band_idx,:);
            %average over frequencies and over protocols
            pow_band = nanmean(pow_band,1); %average over protocols
            pow_band = nanmean(pow_band,3); %average over frequencies
            y = squeeze(pow_band);
            
            if ismember(idata, [3 4 5]) %only for recovery, wor and wor timenorm
                y = movmean(y,cfg{4}.timefreq.movmeanwin(idata),'omitnan');
            end
            
            leg{chan_nr} = plot(x, y, 'Color', color);
        end
        
        %set figure display :
        %set(gca, 'YScale', 'log');
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time');
        
        if idata <= 7
            ylabel(sprintf('%g-%gHz\nPower (mV²)',cfg{irat}.timefreq.foi_band{ifreq}(1),cfg{irat}.timefreq.foi_band{ifreq}(end)));
        else
            ylabel(sprintf('%g-%gHz\nNormalized power',cfg{irat}.timefreq.foi_band{ifreq}(1),cfg{irat}.timefreq.foi_band{ifreq}(end)));
        end
        leg = flip(leg);
        legend_name = flip(legend_name); %set E16 on top
        legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
        legend('boxoff');
        axis tight;
    end
    
    if ~isfolder(cfg{4}.imagesavedir)
        fprintf('Creating directory %s\n', cfg{4}.imagesavedir);
        mkdir(cfg{4}.imagesavedir)
    end
    
    if idata >= 7
        cfg{4}.imagesavedir_data{2}=fullfile(cfg{4}.imagesavedir_data{2},'baseline_corrected');
    end
    
    if idata > 12 & idata < 19
        cfg{4}.imagesavedir_data{2} = fullfile(cfg{4}.imagesavedir_data{2},'Log');
    end
    
    if idata >= 19
        cfg{4}.imagesavedir_data{2} = fullfile(cfg{4}.imagesavedir_data{2},'Log','baseline_corrected');
    end
    
    
    if ~isfolder(cfg{4}.imagesavedir_data{2})
        fprintf('Creating directory %s\n', cfg{4}.imagesavedir_data{2});
        mkdir(cfg{4}.imagesavedir_data{2})
    end
    
    %save figure :
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata},'.png']),'-r600');
    savefig(fig,fullfile(cfg{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata}]));
    close all
    cfg{4}.imagesavedir_data{2}= fullfile(cfg{4}.imagesavedir,'band_depth');
  
    
    %% plot low frequency over high frequency ratio over time for all rats
    
    fig = figure;hold;
    sgtitle([analysis_names{idata}, ' All Rats'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
    
    
        
      
        C        	= colormap(winter(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
        
        for ichan = 1:size(chan_list,1)%chan_name = string(fieldnames(data))'
            chan_name= chan_list{ichan};
            
            %go trhough each electrode in the ascending order
            chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
            legend_name{chan_nr} = chan_name;
            %                 if isempty(data{ifreq}.(chan_name))
            %                     continue
            %                 end
            color   = C(chan_nr,:);
            
            x = data.(chan_name).time;
            clear pow_band y
            
            %find frequencies in the Low range
            lffreq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{1}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{1}(end);
            lfpow_band= data.(chan_name).powspctrm(:,:,lffreq_band_idx,:);
            %average over frequencies and over protocols
            lfpow_band = nanmean(lfpow_band,1); %average over protocols
            lfpow_band = nanmean(lfpow_band,3); %average over frequencies
            A = squeeze(lfpow_band);
            
            %find frequencies in the High range
            hffreq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{3}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{4}(end);
            hfpow_band= data.(chan_name).powspctrm(:,:,hffreq_band_idx,:);
            %average over frequencies and over protocols
            hfpow_band = nanmean(hfpow_band,1); %average over protocols
            hfpow_band = nanmean(hfpow_band,3); %average over frequencies
            B = squeeze(hfpow_band);
            
            y= A ./ B;
            
            clear A B
            
            if ismember(idata, [3 4 5]) %only for recovery, wor and wor timenorm
                y = movmean(y,cfg{4}.timefreq.movmeanwin(idata),'omitnan');
            end
            
            leg{chan_nr} = plot(x, y, 'Color', color);
        end %ichan
        
        
        %set figure display :
        %set(gca, 'YScale', 'log');
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Time');

        ylabel(sprintf('LF/HF ratio'));

        legend_name = flip(legend_name); %set E16 on top
        legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
        legend('boxoff');
        axis tight;

    
    if ~isfolder(cfg{4}.imagesavedir)
        fprintf('Creating directory %s\n', cfg{4}.imagesavedir);
        mkdir(cfg{4}.imagesavedir)
    end
    
    if ~isfolder(cfg{4}.imagesavedir_data{4})
        fprintf('Creating directory %s\n', cfg{4}.imagesavedir);
        mkdir(cfg{4}.imagesavedir_data{4})
    end
    
    if idata >= 7 & idata <= 12
        cfg{4}.imagesavedir_data{4}=fullfile(cfg{4}.imagesavedir_data{4},'baseline_corrected');
    end
    
    if idata > 12 & idata < 19
        cfg{4}.imagesavedir_data{4} = fullfile(cfg{4}.imagesavedir_data{4},'Log');
    end
    
    if idata >= 19
        cfg{4}.imagesavedir_data{4} = fullfile(cfg{4}.imagesavedir_data{2},'Log','baseline_corrected');
    end
    
    
    if ~isfolder(cfg{4}.imagesavedir_data{2})
        fprintf('Creating directory %s\n', cfg{4}.imagesavedir_data{2});
        mkdir(cfg{4}.imagesavedir_data{2})
    end
    
    %save figure :
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(cfg{4}.imagesavedir_data{4},['AllRats_',analysis_names{idata},'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg{4}.imagesavedir_data{4},['AllRats_',analysis_names{idata},'.png']),'-r600');
    savefig(fig,fullfile(cfg{4}.imagesavedir_data{4},['AllRats_',analysis_names{idata}]));
    close all
    cfg{4}.imagesavedir_data{4}= fullfile(cfg{4}.imagesavedir,'LFHF_ratio');  
    
end %idata
end %wod_grandaverage


