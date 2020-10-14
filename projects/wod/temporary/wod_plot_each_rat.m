%% plot all data
    if ~isfolder(cfgcommon.imagesavedir)
        fprintf('Creating directory %s\n', cfgcommon.imagesavedir);
        mkdir(cfgcommon.imagesavedir)
    end
    
    analysis_names = {'Timefreq_WOD', 'Timefreq_WOD_timenorm', 'Timefreq_recovery'};
    data           = {timefreq_wod_grandavg, timefreq_wod_timenorm_grandavg, timefreq_recovery_grandavg};
    
    for idata = 1:size(data,2)
        
        fig = figure;
        sgtitle(analysis_names{idata}, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(config{irat}.timefreq.foi,2)
            subplot(2,2,ifreq);hold;
            
            chan_list 	= fieldnames(timefreq_wod_timenorm_allrats{ifreq});
            C        	= colormap(autumn(numel(chan_list)));
            for ichan = 1:numel(chan_list)
                
                chan_nr          = numel(chan_list)-ichan+1;
                chan_label{ichan} = ['E', num2str(chan_nr)];
                if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                    continue
                end        
       for idata = 1:size(data,2)
        
        fig = figure;
        sgtitle(analysis_names{idata}, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(config{irat}.timefreq.foi,2)
            subplot(2,2,ifreq);hold;
            
            chan_list 	= fieldnames(timefreq_wod_timenorm_allrats{ifreq});
            C        	= colormap(autumn(numel(chan_list)));
            for ichan = 1:numel(chan_list)
                
                chan_nr          = numel(chan_list)-ichan+1;
                chan_label{ichan} = ['E', num2str(chan_nr)];
                if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                    continue
                end        
                %find color for this chan. deep chans are darker
%                 chan_nr = split(chan_label{ichan}, 'E');
%                 chan_nr = str2double(chan_nr{2});
                color   = C(chan_nr,:);
                
                x = data{idata}{ifreq}.(chan_label{ichan}).time;
                y = squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(1,1,1,:));
                
                if idata == 3 %only for recovery
                    y = movmean(y,10,'omitnan');
                end
                
                
                %plot std
%                 std = nanstd(squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(:,1,1,:)));
%                 y_area = [y - std; std; std]'; %FIXME tester avec 2*std
%                 filled_SD = area(x,y_area);
%                 filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.4; filled_SD(3).FaceAlpha = 0.4;
%                 filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
%                 filled_SD(2).FaceColor = color; filled_SD(3).FaceColor = color;
%                 filled_SD(1).ShowBaseLine = 'off';
                
                %find maximum of the first half of data to adapt y scale
                data_max(ichan) = max(y(x>0 & x<x(end)/2)+std(x>0 & x<x(end)/2)); 
                
                %plot avg
                leg{ichan} = plot(x, y, 'Color', color);
            end
            
            %set figure display : 
            %set(gca, 'YScale', 'log');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('Time');
            ylabel(sprintf('%g-%gHz\nNormalized power',config{irat}.timefreq.foi{ifreq}(1),config{irat}.timefreq.foi{ifreq}(end)));
            legend([leg{:}], chan_label{:}, 'Fontsize',8,'location','eastoutside');
            legend('boxoff');
            axis tight;
            if contains(analysis_names{idata}, 'WOD') 
                ylim([0 max(data_max)]);
            end
            
        end
        
        %save figure : 
%         set(fig,'PaperOrientation','landscape');
%         set(fig,'PaperUnits','normalized');
%         set(fig,'PaperPosition', [0 0 1 1]);
%         print(fig, '-dpdf', fullfile(cfgcommon.imagesavedir,['AllRats_',analysis_names{idata},'.pdf']),'-r600');
%         print(fig, '-dpng', fullfile(cfgcommon.imagesavedir,['AllRats_',analysis_names{idata},'.png']),'-r600');
        
       end
       
      for idata = 1:size(data,2)
       for iwod = 1 : size(data{idata}{ifreq}.(chan_label{ichan}).powspctrm,1)
            
           fig = figure;
        sgtitle(analysis_names{idata}, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(config{irat}.timefreq.foi,2)
            subplot(2,2,ifreq);hold;
            
            chan_list 	= fieldnames(timefreq_wod_timenorm_allrats{ifreq});
            C        	= colormap(autumn(numel(chan_list)));
            for ichan = 1:numel(chan_list)
                
                chan_nr          = numel(chan_list)-ichan+1;
                chan_label{ichan} = ['E', num2str(chan_nr)];
                if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                    continue
                end        
       for idata = 1:size(data,2)
        
        fig = figure;
        sgtitle(analysis_names{idata}, 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 23);
        
        for ifreq = 1:size(config{irat}.timefreq.foi,2)
            subplot(2,2,ifreq);hold;
            
            chan_list 	= fieldnames(timefreq_wod_timenorm_allrats{ifreq});
            C        	= colormap(autumn(numel(chan_list)));
            for ichan = 1:numel(chan_list)
                
                chan_nr          = numel(chan_list)-ichan+1;
                chan_label{ichan} = ['E', num2str(chan_nr)];
                if isempty(data{idata}{ifreq}.(chan_label{ichan}))
                    continue
                end        
                %find color for this chan. deep chans are darker
%                 chan_nr = split(chan_label{ichan}, 'E');
%                 chan_nr = str2double(chan_nr{2});
                color   = C(chan_nr,:);
                
                x = data{idata}{ifreq}.(chan_label{ichan}).time;
                y = squeeze(data{idata}{ifreq}.(chan_label{ichan}).powspctrm(iwod,1,1,:));
            end
        end
       end
            end
        end
       end
      end
            
        
    
    
    
    close all

    
    % compute statistiques : y a-t-il des profondeurs qui récupèrent plus vite ? Lesquelles ?
    
    % trouver période stable de recovery
    
end %slurm_task_id == 0

        end
    end