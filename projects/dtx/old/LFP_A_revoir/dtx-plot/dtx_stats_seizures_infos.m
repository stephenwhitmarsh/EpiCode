function dtx_stats_seizures_infos(cfg,Seizure_Infos)
%plot seizures info

%Bug : need to be fixed
if length(Seizure_Infos) > 1
    if isempty(Seizure_Infos{1})
        Seizure_Infos = Seizure_Infos{end};
    end
end

pat_idx = 0;
for ipatient = 1:size(Seizure_Infos, 2)
    if ~isempty(Seizure_Infos{ipatient})
        pat_idx = pat_idx +1;
        pat_name{pat_idx} = Seizure_Infos{ipatient}.ID;
        pat_name{pat_idx} = strrep(pat_name{pat_idx},'_',' ');
        pat_name{pat_idx} = strrep(pat_name{pat_idx},'-',' ');
        pat_name{pat_idx} = strrep(pat_name{pat_idx},' MERGED',[]);
        if length(pat_name{pat_idx}) >=12
            pat_name{pat_idx} = pat_name{pat_idx}(1:12);
        end
    end
end

%Temps entre 2 crises
for ivalue_name = ["timeBetween2SlowWaves","timeBetween2SlowWaves_R","timeBetween2SlowWaves_L"]
    if isfield(Seizure_Infos{1}, ivalue_name)

        pat_idx = 0;
        values = [];
        cv2=[];
        for ipatient = 1:size(Seizure_Infos,2)
            %gel values of all patients
            if ~isempty(Seizure_Infos{ipatient})
                pat_idx                = pat_idx + 1;
                values{pat_idx}        = minutes(Seizure_Infos{ipatient}.(ivalue_name));
                x_values{pat_idx}      = minutes(Seizure_Infos{ipatient}.(sprintf('x_%s',ivalue_name)));
                meanvalues{pat_idx}    = nanmean(values{pat_idx});
                
                %CV2
                if length(values{pat_idx}) > 2
                    cv2.hascv2(pat_idx)         = true;
                    for i = 1:length(values{pat_idx})-1
                        cv2.values{pat_idx}(i)           = 2*abs(values{pat_idx}(i)-values{pat_idx}(i+1))/(values{pat_idx}(i)+values{pat_idx}(i+1));
                    end
                    cv2.x_values{pat_idx}       = minutes(Seizure_Infos{ipatient}.(sprintf('x_%s',ivalue_name))(2:end));
                else
                    cv2.values{pat_idx}         = NaN;
                    cv2.x_values{pat_idx}       = NaN;
                    cv2.hascv2(pat_idx)         = false;
                end
                cv2.meancv2{pat_idx}        = nanmean(cv2.values{pat_idx});
                
                
            end
        end
        
        %scatter plot time between 2 SW
        fig1 = figure;
        hold;
        for ipat_idx = 0:size(values,2)-1
            meandata = mean(values{ipat_idx+1});
            stddata = std(values{ipat_idx+1});
            xdata = (rand(1,size(values{ipat_idx+1},2)))*0.2+0.9+0.5*ipat_idx; %random values between 0.9 and 1.1
            %plot data
            scatter(xdata,values{ipat_idx+1},'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
            %plot mean
            plot([0.8+0.5*ipat_idx, 1.2+0.5*ipat_idx],[meandata, meandata], 'k');
            %plot vertical error bars
            plot([1+0.5*ipat_idx, 1+0.5*ipat_idx], [meandata, meandata+stddata],'k', [1+0.5*ipat_idx, 1+0.5*ipat_idx], [meandata, meandata-stddata], 'k');
            %plot horizontal error bars
            plot([0.9+0.5*ipat_idx, 1.1+0.5*ipat_idx], [meandata+stddata, meandata+stddata],'k', [0.9+0.5*ipat_idx, 1.1+0.5*ipat_idx], [meandata-stddata, meandata-stddata], 'k');
            xlim([0.7, 1.3+0.5*ipat_idx]);
        end
        y = ylim;
        ylim([0 y(2)]);
        xticks(1:0.5:size(values,2)*2);
        ylabel('Time (minutes)');
        set(gca,'Fontsize',12);
        set(gca, 'XTickLabel',pat_name);
        title(sprintf('%s : \n%g +/- %g (avg of avg)',ivalue_name, nanmean([meanvalues{:}]), nanstd([meanvalues{:}])), 'Fontsize', 18, 'Interpreter', 'none');
        set(gca,'FontWeight','bold' );
        set(gca,'TickDir','out');

        %scatter plot CV2 of time between 2 SW
        fig2 = figure;
        hold;
        for ipat_idx = 0:size(cv2.values,2)-1
            meandata = mean(cv2.values{ipat_idx+1});
            stddata = std(cv2.values{ipat_idx+1});
            xdata = (rand(1,size(cv2.values{ipat_idx+1},2)))*0.2+0.9+0.5*ipat_idx; %random values between 0.9 and 1.1
            %plot data
            scatter(xdata,cv2.values{ipat_idx+1},'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
            %plot mean
            plot([0.8+0.5*ipat_idx, 1.2+0.5*ipat_idx],[meandata, meandata], 'k');
            %plot vertical error bars
            plot([1+0.5*ipat_idx, 1+0.5*ipat_idx], [meandata, meandata+stddata],'k', [1+0.5*ipat_idx, 1+0.5*ipat_idx], [meandata, meandata-stddata], 'k');
            %plot horizontal error bars
            plot([0.9+0.5*ipat_idx, 1.1+0.5*ipat_idx], [meandata+stddata, meandata+stddata],'k', [0.9+0.5*ipat_idx, 1.1+0.5*ipat_idx], [meandata-stddata, meandata-stddata], 'k');
            xlim([0.7, 1.3+0.5*ipat_idx]);
        end
        ylim([0 2]);
        xticks(1:0.5:size(values,2)*2);
        ylabel('Time (minutes)');
        set(gca,'Fontsize',12);
        set(gca, 'XTickLabel',pat_name);
        title(sprintf('CV2 of %s : \n%g +/- %g (avg of avg)',ivalue_name, nanmean([cv2.meancv2{:}]), nanstd([cv2.meancv2{:}])), 'Fontsize', 18, 'Interpreter', 'none');
        set(gca,'FontWeight','bold' );
        set(gca,'TickDir','out');
        
        %Interseizure-time over time
        fig3 = figure;
        hold;
        for i = 1:size(values, 2)
            plot(x_values{i}, values{i},'o', 'LineWidth', 2);
        end
        ylabel('Time between 2 seizures (minuts)');
        xlabel('Time from begining of recording (minuts)');
        set(gca,'Fontsize',12);
        title(sprintf('%s over time',ivalue_name), 'Fontsize', 18, 'Interpreter', 'none');
        set(gca,'FontWeight','bold' );
        set(gca,'TickDir','out');
        legend(pat_name{:},'location','northeast','fontsize',10);
        
        %CV2 over time
        fig4 = figure;
        hold;
        for ipatient = 1:size(values, 2)
            if length(cv2.x_values{ipatient}) > 2
                plot(cv2.x_values{ipatient}, cv2.values{ipatient},'o', 'LineWidth', 2);
            end
        end
        ylim([0 2]);
        ylabel('CV2 of time between 2 seizures');
        xlabel('Time from begining of recording (minuts)');
        set(gca,'Fontsize',12);
        title(sprintf('CV2 of %s over time',ivalue_name), 'Fontsize', 18, 'Interpreter', 'none');
        set(gca,'FontWeight','bold' );
        set(gca,'TickDir','out');
        legend(pat_name{cv2.hascv2},'location','northeast','fontsize',10);
        
        
        
        %save figures
        if ~(exist(cfg.imagesavedir)==7)
            mkdir(cfg.imagesavedir);
            fprintf('Create forlder %s',cfg.imagesavedir);
        end
        
        set(fig1,'PaperOrientation','landscape');
        set(fig1,'PaperUnits','normalized');
        set(fig1,'PaperPosition', [0 0 1 1]);
        set(fig1,'Renderer','Painters');
        print(fig1, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_seizures_timings']),'-r600');
        print(fig1, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_seizures_timings']),'-r600');
        
        set(fig2,'PaperOrientation','landscape');
        set(fig2,'PaperUnits','normalized');
        set(fig2,'PaperPosition', [0 0 1 1]);
        set(fig2,'Renderer','Painters');
        print(fig2, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_seizures_timings_CV2']),'-r600');
        print(fig2, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_seizures_timings_CV2']),'-r600');
        
        set(fig3,'PaperOrientation','landscape');
        set(fig3,'PaperUnits','normalized');
        set(fig3,'PaperPosition', [0 0 1 1]);
        set(fig3,'Renderer','Painters');
        print(fig3, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_interseizure_overtime']),'-r600');
        print(fig3, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_interseizure_overtime']),'-r600');
        
        set(fig4,'PaperOrientation','landscape');
        set(fig4,'PaperUnits','normalized');
        set(fig4,'PaperPosition', [0 0 1 1]);
        set(fig4,'Renderer','Painters');
        print(fig4, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_interseizure_overtime_CV2']),'-r600');
        print(fig4, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,char(ivalue_name),'_interseizure_overtime_CV2']),'-r600');
        close all
    end
        
end

