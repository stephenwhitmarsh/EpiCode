function wod_tfr_plotrat(cfg)

analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_recovery','timefreq_wor','timefreq_wor_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_recovery_blcorrected','timefreq_wor_blcorrected','timefreq_wor_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_recovery','log_timefreq_wor','log_timefreq_wor_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected', 'log_timefreq_recovery_blcorrected','log_timefreq_wor_blcorrected','log_timefreq_wor_timenorm_blcorrected', 'log_timefreq_baseline_blcorrected'};

%go through each data structure
for idata = 1:size( analysis_names,2)
    %% TFR by chan by rat
      %load data
        data_temp = load(fullfile(cfg.datasavedir,[cfg.prefix,analysis_names{idata},'.mat']));
        
        data_rat = [];
        count_trials                    = 0;
        %get data from one rat in a structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list = fieldnames(data_temp.(analysis_names{idata}){itrial});
            
            data_rat= data_temp.(analysis_names{idata}){itrial};
            
            for ichan= 1: numel(chan_list)
                chan_name = chan_list{ichan};
                %replace empty channels by nans
                hasdata = find(~structfun(@isempty,data_rat));
                if isempty(hasdata)
                    continue
                end
                
                if isempty(data_rat.(chan_name))
                    data_rat.(chan_name)                       = data_rat.(chan_list{hasdata(1)});
                    data_rat.(chan_name).powspctrm             = ones(size(data_rat.(chan_name).powspctrm));
                end
                
                data_plot = data_rat.(chan_name);
                
                
                fig=figure;
                
                sgtitle([analysis_names{idata},chan_name,sprintf('Rat %s',cfg.prefix),' ',sprintf('WOD %d/%d\n',wod_rat(count_trials),size(data_temp.(analysis_names{idata}),2))], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
                subplot(2,1,1)
                %plot TFR
                %voir les paramètres optionnels dans le descriptifs de la fonction pour
                %modifier l'aspect du TFR. Avec les paramètres par défaut :
                cfgtemp         = [];
                cfgtemp.channel = 'all';
                cfgtemp.interactive = 'no';
                cfgtemp.colormap= 'jet';
                cfgtemp.fontsize = 12;
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
                cfgtemp.fontsize = 12;
                ft_singleplotTFR(cfgtemp, data_plot);
                ft_pimpplot(fig, jet(5000))
                
                
                title('Time-Frequency [1 : 50]');
                %figure settings (fontsize,fontweight, ticks dir, renderer etc.)
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time (s)');
                ylim([1 50]);
                ylabel(sprintf('Frequency (Hz)'));
                axis tight;
                
                if ~isfolder(cfg.imagesavedir)
                    fprintf('Creating directory %s\n', cfg.imagesavedir);
                    mkdir(cfg.imagesavedir)
                end
                
                rat_imagedir= fullfile(cfg.imagesavedir_data{1},'by_rat',sprintf(cfg.prefix));
                
                
                if idata >= 7
                    rat_imagedir = fullfile(rat_imagedir, 'baseline_corrected');
                end
                
                if idata > 12 & idata < 19
                    rat_imagedir = fullfile(rat_imagedir, 'Log');
                end
                
                if idata >= 19
                    rat_imagedir = fullfile(rat_imagedir, 'Log','baseline corrected');
                end
                
                
                
                if ~isfolder(rat_imagedir)
                    fprintf('Creating directory %s\n',  rat_imagedir);
                    mkdir(rat_imagedir)
                end
                
                
                fname= fullfile(rat_imagedir,sprintf('Rat_%s_WoD_%d_%s_%s.pdf',cfg.prefix,wod_rat(count_trials),analysis_names{idata}, chan_name));
                dtx_savefigure(fig,fname,'pdf','png','close');
                
                
            end %ichan
        end %itrial
        clear data_plot fname
    
    
    %% Plot LF/HF ratio by rat

        data_rat = [];
        count_trials                    = 0;
        %get data from one rat in a structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});
            
            data_rat= data_temp.(analysis_names{idata}){itrial};

            
            fig = figure;hold;
            C        	= colormap(winter(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels

            sgtitle([analysis_names{idata},sprintf('Rat %s',cfg.prefix),' ',sprintf('WOD %d/%d\n',wod_rat(count_trials),size(data_temp.(analysis_names{idata}),2))], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
            
            for ichan = 1:size(chan_list,1)%chan_name = string(fieldnames(data))'

                chan_name= chan_list{ichan};
                
                %go trhough each electrode in the ascending order
                chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
                legend_name{chan_nr} = chan_name;
                %                 if isempty(data{ifreq}.(chan_name))
                %                     continue
                %                 end
                color   = C(chan_nr,:);
                
                if isempty(data_rat.(chan_name))
                    continue
                end
                
                x = data_rat.(chan_name).time;
                
                %find frequencies in the Low range
                lffreq_band_idx = data_rat.(chan_name).freq >= cfg.timefreq.foi_band{1}(1) & data_rat.(chan_name).freq <= cfg.timefreq.foi_band{1}(end);
                lfpow_band= data_rat.(chan_name).powspctrm(:,lffreq_band_idx,:);
                %average over frequencies
                lfpow_band = nanmean(lfpow_band,2); %average over frequencies
                A = squeeze(lfpow_band);
                
                %find frequencies in the High range
                hffreq_band_idx = data_rat.(chan_name).freq >= cfg.timefreq.foi_band{4}(1) & data_rat.(chan_name).freq <= cfg.timefreq.foi_band{4}(end);
                hfpow_band= data_rat.(chan_name).powspctrm(:,hffreq_band_idx,:);
                %average over frequencies
                hfpow_band = nanmean(hfpow_band,2); %average over frequencies
                B = squeeze(hfpow_band);
                
                y= A ./ B;
                
                clear A B
                
                if ismember(idata, [3 4 5]) %only for recovery, wor and wor timenorm
                    y = movmean(y,cfg.timefreq.movmeanwin(idata),'omitnan');
                end
                
                leg{chan_nr} = plot(x, y, 'Color', color);
            end %ichan
            
            
            %set figure display :
            %set(gca, 'YScale', 'log');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            
            ylabel(sprintf('LF/HF ratio'));
            
            legend_name = flip(legend_name); %set E16 on top
            legend( legend_name{:}, 'Fontsize',8,'location','eastoutside');
            legend('boxoff');
            axis tight;
            
            imagessavedir=fullfile(cfg.imagesavedir_data{4},'by rat',sprintf('%s',cfg.prefix));
            
            
            if ~isfolder(cfg.imagesavedir)
                fprintf('Creating directory %s\n', cfg.imagesavedir);
                mkdir(cfg.imagesavedir)
            end
            
            if ~isfolder(imagessavedir)
                fprintf('Creating directory %s\n', imagessavedir);
                mkdir(imagessavedir)
            end
            
            if idata >= 7 & idata <= 12
                imagessavedir=fullfile(imagessavedir,'baseline_corrected');
            end
            
            if idata > 12 & idata < 19
                imagessavedir = fullfile(imagessavedir,'Log');
            end
            
            if idata >= 19
                imagessavedir = fullfile(imagessavedir,'Log','baseline_corrected');
            end
            
            
            if ~isfolder(cfg.imagesavedir_data{2})
                fprintf('Creating directory %s\n', cfg.imagesavedir_data{2});
                mkdir(cfg.imagesavedir_data{2})
            end
            
            %save figure :
            
            fname=fullfile(imagessavedir,sprintf('%sWOD_%i/%i_%s',cfg.prefix,itrial,numel(wod_rat),analysis_names{idata}));
            dtx_savefigure(fig,fname,'pdf','png','close');

        end %itrial
end %idata
