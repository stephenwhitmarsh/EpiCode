function wod_plot

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/development'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
 
end

ft_defaults

config = wod_setparams;




%load data
analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_recovery','timefreq_wor','timefreq_wor_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_recovery_blcorrected','timefreq_wor_blcorrected','timefreq_wor_timenorm_blcorrected', 'timefreq_baseline_blcorrected'};

%go through each data structure
for idata = 1:size( analysis_names,2)
    
    fprintf('For %s\n', analysis_names{idata});
    
    
    cfgcommon = config{4}; %same common config for all rats
    
    %initialize data structures
    data_allrats = [];
    count_trials                    = 0;
    
    % load timefreq data from all rats
    for irat = 4:size(config,2)
        %load data
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        data_temp = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis_names{idata},'.mat']));
        
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
    
    clear data_allrats
    
    
    
    
    
    %% TFR average by chan and over rats
    
    %gather powspctrm data for each channel
    %             powspctrm_allchans = nan(size(data.(chan_name).powspctrm));
    %             ichan=0;
    %             for chan_name = string(fieldnames(data))'
    %                 ichan=ichan+1;
    %                 temp = data.(chan_name).powspctrm;
    %                 if size(data.(chan_name).powspctrm,4) > size(powspctrm_allchans,4)
    %                     temp = data.(chan_name).powspctrm(:,:,:,1:size(powspctrm_allchans,4));
    %                 elseif size(data.(chan_name).powspctrm,4) < size(powspctrm_allchans,4)
    %                     powspctrm_allchans = powspctrm_allchans(:,:,:,1:size(data.(chan_name).powspctrm,4));
    %                 end
    %                 powspctrm_allchans(:,ichan,:,:) = temp;
    %                 clear temp
    %             end
    %
    %             %re-create fieldtrip structure
    %             datafreq_allchans.label = fieldnames(data)';
    %             datafreq_allchans.powspctrm = powspctrm_allchans;
    %             datafreq_allchans.freq = data.(chan_name).freq;
    %             datafreq_allchans.time = data.(chan_name).time(1:size(powspctrm_allchans,4));
    %             datafreq_allchans.dimord = 'subj_chan_freq_time';
    %             clear powspctrm_allchans
    
    
    
    
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
        
        if ~isfolder(config{4}.imagesavedir)
            fprintf('Creating directory %s\n', config{4}.imagesavedir);
            mkdir(config{4}.imagesavedir)
        end
        
        
        if idata >= 7
            config{4}.imagesavedir_data{1} = fullfile(config{4}.imagesavedir_data{1}, 'baseline_corrected');
        end
        
        if ~isfolder(config{4}.imagesavedir_data{1})
            fprintf('Creating directory %s\n', config{4}.imagesavedir);
            mkdir(config{4}.imagesavedir_data{1})
        end
        
        
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(config{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s.pdf', analysis_names{idata}, chan_name)),'-r600');
        print(fig, '-dpng', fullfile(config{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s.png', analysis_names{idata}, chan_name)),'-r600');
        savefig(fig,fullfile(config{4}.imagesavedir_data{1},sprintf('AllRats_TFR_%s_%s', analysis_names{idata}, chan_name)));
        close all
        
        config{4}.imagesavedir_data{1}= fullfile(config{4}.imagesavedir,'TFR');
        
        clear data_plot
        
        
        
        
        
    end %ichan
    
    %% TFR by chan by rat
    
    
    % load timefreq data from all rats
    for irat = 4:size(config,2)
        %load data
        
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(config,2));
        data_temp = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis_names{idata},'.mat']));
        
        data_rat = [];
        count_trials                    = 0;
        %get data from one rat in a structure
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});

            data_rat= data_temp.(analysis_names{idata}){itrial};
            
            for ichan= 1: numel(chan_list)
                chan_name = chan_list{ichan};
                %replace empty channels by nans
                hasdata = find(~structfun(@isempty,data_rat));
                if isempty(data_rat.(chan_name))
                    data_rat.(chan_name)                       = data_rat.(chan_list{hasdata(1)});
                    data_rat.(chan_name).powspctrm             = zeros(size(data_rat.(chan_name).powspctrm));
                end
                
                data_plot = data_rat.(chan_name);
                
               
                fig=figure;
                
                sgtitle([analysis_names{idata},chan_name,sprintf('Rat %d/%d',irat,size(config,2)),' ',sprintf('WOD %d/%d\n',wod_rat(count_trials),size(data_temp.(analysis_names{idata}),2))], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
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
                
                if ~isfolder(config{4}.imagesavedir)
                    fprintf('Creating directory %s\n', config{4}.imagesavedir);
                    mkdir(config{4}.imagesavedir)
                end
                
                rat_imagedir= fullfile(config{4}.imagesavedir_data{1},'by_rat',sprintf(config{irat}.prefix));
                
                
                if idata >= 7
                    rat_imagedir = fullfile(rat_imagedir, 'baseline_corrected');
                end
                
                if ~isfolder(rat_imagedir)
                    fprintf('Creating directory %s\n',  rat_imagedir);
                    mkdir(rat_imagedir)
                end
                
                
                
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.pdf',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)),'-r600');
                print(fig, '-dpng', fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.png',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)),'-r600');
                %savefig(fig,fullfile(rat_imagedir,sprintf('Rat_%d_WoD_%d_%s_%s.pdf',irat,wod_rat(count_trials),analysis_names{idata}, chan_name)));
                close all
                
                
                
                
                
                
            end %ichan
            
            
            
            
        end %itrial
        clear data_temp data_plot
    end %irat

            %
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
                        y = movmean(y,config{4}.timefreq.movmeanwin(idata),'omitnan');
                    end

                    %plot std
                    %                 std = nanstd(squeeze(data{ifreq}.(chan_name).powspctrm(:,1,1,:)));
                    %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
                    %                 filled_SD = area(x,y_area);
                    %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
                    %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
                    %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
                    %                 filled_SD(1).ShowBaseLine = 'off';

                    %plot avg

                    leg{chan_nr} = plot(x, y, 'Color', color);
                end

                %set figure display :
                %set(gca, 'YScale', 'log');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel('Time');
                
                if idata <= 7
                    ylabel(sprintf('%g-%gHz\nPower (mV²)',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(end)));
                else
                    ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(end)));
                end
                leg = flip(leg); legend_name = flip(legend_name); %set E16 on top
                legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
                legend('boxoff');
                axis tight;
            end
            
            if ~isfolder(config{4}.imagesavedir)
                    fprintf('Creating directory %s\n', config{4}.imagesavedir);
                    mkdir(config{4}.imagesavedir)
                end
            
            if idata >= 7
                config{4}.imagesavedir_data{2}=fullfile(config{4}.imagesavedir_data{2},'baseline_corrected');
            end
            if ~isfolder(config{4}.imagesavedir_data{2})
                    fprintf('Creating directory %s\n', config{4}.imagesavedir_data{2});
                    mkdir(config{4}.imagesavedir_data{2})
                end
            
            %save figure :
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(config{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata},'.pdf']),'-r600');
            print(fig, '-dpng', fullfile(config{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata},'.png']),'-r600');
            savefig(fig,fullfile(config{4}.imagesavedir_data{2},['AllRats_',analysis_names{idata}]));
            close all
            config{4}.imagesavedir_data{2}= fullfile(config{4}.imagesavedir,'band_depth');

    
end %idata
end %wod_plot
%             config{4}.imagesavedir_data{1}='\\lexport\iss01.charpier\analyses\wod\TFR';
%
%             %% TFR average over chans by rat
%             %get wod_rat and wod_rat_number for titles
%             %gather powspctrm data for each channel
%
%
%             %% plot data by frequency band by protocol
%             for iwod = 1 :size(data.(chan_name).powspctrm,1)
%
%                 fig = figure;
%                 rat_wod_sum = sum(wod_rat == wod_rat(iwod));
%                 rat_wod_nr = 1;
%                 if iwod >1
%                     if wod_rat(iwod-1) == wod_rat(iwod)
%                         rat_wod_nr = 2;
%                     end
%                 end
%                 sgtitle(sprintf('%s Rat %d (%d/%d)',analysis_names{idata},wod_rat(iwod), rat_wod_nr, rat_wod_sum), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
%
%                 for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
%
%                     subplot(2,2,ifreq);hold;
%                     C        	= colormap(autumn(numel(chan_list)));%fais un dégradé codé en RGB avec autant de couleurs que le channels
%
%                     for chan_name = string(fieldnames(data))'
%                         %go trhough each electrode in the ascending order
%                         chan_nr          = str2num(cell2mat(regexp(chan_name,'\d*','Match')))-7;%-7 because E8 is 1st and E16 is 9th
%                         legend_name{chan_nr} = chan_name;
%                         color   = C(chan_nr,:);
%
%                         x = data.(chan_name).time;
%                         %select wod and frequency band
%                         freq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
%                         pow_band= data.(chan_name).powspctrm(iwod,:,freq_band_idx,:);
%                         %average over frequencies and over protocols
%                         pow_band = nanmean(pow_band,3); %average over frequencies
%                         y = squeeze(pow_band);
%                         y = movmean(y,config{4}.timefreq.movmeanwin(idata),'omitnan');
%
%
%                         %plot std
%                         %                 std = nanstd(squeeze(data{ifreq}.(chan_name).powspctrm(:,1,1,:)));
%                         %                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
%                         %                 filled_SD = area(x,y_area);
%                         %                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%                         %                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%                         %                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
%                         %                 filled_SD(1).ShowBaseLine = 'off';
%
%                         %find maximum of the first half of data to adapt y scale
%                         leg{chan_nr} = plot(x, y, 'Color', color);
%                     end
%
%                     %set figure display :
%                     %set(gca, 'YScale', 'log');
%                     set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
%                     xlabel('Time');
%                     ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(2)));
%                     leg = flip(leg); legend_name = flip(legend_name); %set E16 on top
%                     legend([leg{:}], legend_name{:}, 'Fontsize',8,'location','eastoutside');
%                     legend('boxoff');
%                     axis tight;
%
%
%                 end
%
%                 if ~isfolder(fullfile(config{4}.imagesavedir_data{2},'by_rat'))
%                               mkdir(fullfile(config{4}.imagesavedir_data{2},'by_rat'));
%                               analydir= fullfile(config{4}.imagesavedir_data{2},'by_rat');
%                 end
%                           if dobaseline
%                               analydir=fullfile(cfgcommon.imagesavedir,'band_depth','by_rat');
%                               if ~isfolder(analydir)
%                                   mkdir(analydir);
%                               end
%                           end
%
%                 %save figure :
%                 set(fig,'PaperOrientation','landscape');
%                 set(fig,'PaperUnits','normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(analydir,sprintf('%s_Rat%d_%d.pdf',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
%                 print(fig, '-dpng', fullfile(analydir,sprintf('%s_Rat%d_%d.png',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
%                 savefig(fig,fullfile(analydir,sprintf('%s_Rat%d_%d',analysis_names{idata},wod_rat(iwod), rat_wod_nr)));
%
%                 close all
%             end
%
%             % compute statistiques : y a-t-il des profondeurs qui récupèrent plus vite ? Lesquelles ?
%
%             % trouver période stable de recovery
%
%             %% plot frequency band averaged over channels over rats
%             fig = figure;
%             sgtitle([analysis_names{idata}, ' All Rats, average over cortical channels'], 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
%
%             for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
%
%                 subplot(2,2,ifreq);hold;
%                 color = 'b';
%
%                 ichan = 0;
%                 for chan_name = string(fieldnames(data))'
%                     ichan = ichan+1;
%
%                     %find frequencies of the band, averaged over all wod (comparer des valeurs : > <
%                     %>= <= ==)
%                     freq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
%                     pow_band= data.(chan_name).powspctrm(:,:,freq_band_idx,:);
%                     %average over frequencies and over protocols
%                     pow_band = nanmean(pow_band,1); %average over protocols
%                     pow_band = nanmean(pow_band,3); %average over frequencies
%                     y = squeeze(pow_band);
%
%                     if ismember(idata, [3 4 5]) %only for recovery
%                         y = movmean(y,config{4}.timefreq.movmeanwin(idata),'omitnan');
%                     end
%
%                     data_chan.time{ichan} = data.(chan_name).time;
%                     data_chan.trial{ichan} = y';
%                     data_chan.label{1} = 'dummy';
%                 end
%
%                 data_avgchan = ft_timelockanalysis([],data_chan);
%
%                 x = data_avgchan.time;
%                 y = data_avgchan.avg;
%
%                 %plot std
%                 y_std = sqrt(data_avgchan.var);
%                 y_area = [y - y_std; y_std; y_std]'; %FIXME tester avec 2*std
%                 filled_SD = area(x,y_area);
%                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
%                 filled_SD(1).ShowBaseLine = 'off';
%
%                 %plot avg
%                 plot(x, y, 'Color', color);
%
%
%
%                 %set figure display :
%                 %set(gca, 'YScale', 'log');
%                 set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
%                 xlabel('Time (s)');
%                 ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(end)));
%                 axis tight;
%
%
%             end
%
%             if ~isfolder(fullfile(config{4}.imagesavedir_data{3}))
%                               mkdir(fullfile(config{4}.imagesavedir_data{3}));
%                               analydir= fullfile(config{4}.imagesavedir_data{3});
%                           end
%
%             %save figure :
%             set(fig,'PaperOrientation','landscape');
%             set(fig,'PaperUnits','normalized');
%             set(fig,'PaperPosition', [0 0 1 1]);
%             print(fig, '-dpdf', fullfile(analydir,['AllRats_AllCortex_',analysis_names{idata},'.pdf']),'-r600');
%             print(fig, '-dpng', fullfile(analydir,['AllRats_AllCortex_',analysis_names{idata},'.png']),'-r600');
%             savefig(fig,fullfile(analydir,['AllRats_AllCortex_',analysis_names{idata}]));
%             close(fig);
%
%             %% plot frequency band averaged over channels by rat
%             for iwod= 1:size(data.(chan_name).powspctrm,1)
%                 fig = figure;
%                 sgtitle(sprintf('%s Rat %d (%d/%d) average over all cortical channels',analysis_names{idata},wod_rat(iwod), rat_wod_nr, rat_wod_sum), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
%
%                 for ifreq = 1:size(cfgcommon.timefreq.foi_band,2)
%
%                     subplot(2,2,ifreq);hold;
%                     color = 'b';
%
%                     ichan = 0;
%                     for chan_name = string(fieldnames(data))'
%                         ichan = ichan+1;
%
%                         %find frequencies of the band, averaged over all wod (comparer des valeurs : > <
%                         %>= <= ==)
%                         freq_band_idx = data.(chan_name).freq >= cfgcommon.timefreq.foi_band{ifreq}(1) & data.(chan_name).freq <= cfgcommon.timefreq.foi_band{ifreq}(2);
%                         pow_band= data.(chan_name).powspctrm(iwod,:,freq_band_idx,:);
%                         %average over frequencies
%                         pow_band = nanmean(pow_band,3);
%                         y = squeeze(pow_band);
%
%                         if ismember(idata, [3 4 5]) %only for recovery
%                             y = movmean(y,config{4}.timefreq.movmeanwin(idata),'omitnan');
%                         end
%
%                         data_chan.time{ichan} = data.(chan_name).time;
%                         data_chan.trial{ichan} = y';
%                         data_chan.label{1} = 'dummy';
%                     end
%
%                     data_avgchan = ft_timelockanalysis([],data_chan);
%
%                     x = data_avgchan.time;
%                     y = data_avgchan.avg;
%
%                     %plot std
%                     y_std = sqrt(data_avgchan.var);
%                     y_area = [y - y_std; y_std; y_std]'; %FIXME tester avec 2*std
%                     filled_SD = area(x,y_area);
%                     filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%                     filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%                     filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
%                     filled_SD(1).ShowBaseLine = 'off';
%
%                     %plot avg
%                     plot(x, y, 'Color', color);
%
%
%
%                     %set figure display :
%                     %set(gca, 'YScale', 'log');
%                     set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
%                     xlabel('Time (s)');
%                     legend
%                     ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi_band{ifreq}(1),config{irat}.timefreq.foi_band{ifreq}(end)));
%                     axis tight;
%
%                 end
%
%                 if ~isfolder(fullfile(config{4}.imagesavedir_data{3},'by_rat'))
%                               mkdir(fullfile(config{4}.imagesavedir_data{3},'by_rat'));
%                               analydir= fullfile(config{4}.imagesavedir_data{3},'by_rat');
%                           end
%
%                 %save figure :
%                 set(fig,'PaperOrientation','landscape');
%                 set(fig,'PaperUnits','normalized');
%                 set(fig,'PaperPosition', [0 0 1 1]);
%                 print(fig, '-dpdf', fullfile(analydir,sprintf('%s_Allcortex_Rat%d_%d.pdf',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
%                 print(fig, '-dpng', fullfile(analydir,sprintf('%s_Allcortex_Rat%d_%d.png',analysis_names{idata},wod_rat(iwod), rat_wod_nr)),'-r600');
%                 savefig(fig,fullfile(analydir,sprintf('%s_Allcortex_Rat%d_%d',analysis_names{idata},wod_rat(iwod), rat_wod_nr)));
%                 close(fig);
%             end %iwod
%             clear data
%         end %dobaseline
%     end %idata
% end %slurm task id == 0

