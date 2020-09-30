timings = Seizure_Infos{2};

%get values
for irat = 1:5
    max_idx = min([length(timings{irat}.timeBetween2SlowWaves), length(timings{irat}.seizureLength)]);
    interictal_time{irat} = timings{irat}.timeBetween2SlowWaves(1:max_idx) - timings{irat}.seizureLength(1:max_idx);
end

% get patient names
pat_idx = 0;
for ipatient = 1:size(Seizure_Infos{2}, 2)
    if ~isempty(Seizure_Infos{2}{ipatient})
        pat_idx = pat_idx +1;
        pat_name{pat_idx} = Seizure_Infos{2}{ipatient}.ID;
        pat_name{pat_idx} = strrep(pat_name{pat_idx},'_',' ');
        pat_name{pat_idx} = strrep(pat_name{pat_idx},'-',' ');
        pat_name{pat_idx} = strrep(pat_name{pat_idx},' MERGED',[]);
        if length(pat_name{pat_idx}) >=12
            pat_name{pat_idx} = pat_name{pat_idx}(1:12);
        end
    end
end


pat_idx = 0;
values = [];
cv2=[];
for ipatient = 1:size(Seizure_Infos{2},2)
    %gel values of all patients
    if ~isempty(Seizure_Infos{2}{ipatient})
        pat_idx                = pat_idx + 1;
        values{pat_idx}        = minutes(interictal_time{ipatient});
        meanvalues{pat_idx}    = nanmean(values{pat_idx});
        
        %CV2
        if length(values{pat_idx}) > 2
            cv2.hascv2(pat_idx)         = true;
            for i = 1:length(values{pat_idx})-1
                cv2.values{pat_idx}(i)           = 2*abs(values{pat_idx}(i)-values{pat_idx}(i+1))/(values{pat_idx}(i)+values{pat_idx}(i+1));
            end
%             cv2.x_values{pat_idx}       = minutes(Seizure_Infos{ipatient}.(sprintf('x_%s',ivalue_name))(2:end));
        else
            cv2.values{pat_idx}         = NaN;
%             cv2.x_values{pat_idx}       = NaN;
            cv2.hascv2(pat_idx)         = false;
        end
        cv2.meancv2{pat_idx}        = nanmean(cv2.values{pat_idx});
        
        
    end
end

%scatter plot interictal time
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
set(gca, 'XTickLabel', pat_name);
title(sprintf('Interictal durations : \n%g +/- %g (avg of avg)', nanmean([meanvalues{:}]), nanstd([meanvalues{:}])), 'Fontsize', 18, 'Interpreter', 'none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');