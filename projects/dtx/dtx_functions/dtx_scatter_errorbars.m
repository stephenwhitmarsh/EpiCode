function fig = dtx_scatter_errorbars(values, xlabel_str, ylabel_str, title_str)
% values is an array with a line per plot to do

fig = figure;
hold;

for iplot = 0:size(values,1)-1
%     data = Seizure_Infos{2}{1}.timeBetween2SlowWaves;
%     data = seconds(data);
    meandata = mean(values(iplot+1,:));
    stddata = std(values(iplot+1,:));
    xdata = (rand(1,size(values(iplot+1,:),2)))*0.2+0.9+0.5*iplot; %random values between 0.9 and 1.1
    %plot data
    scatter(xdata,values(iplot+1,:),'filled','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]);
    %plot mean
    plot([0.8+0.5*iplot, 1.2+0.5*iplot],[meandata, meandata], 'k');
    %plot vertical error bars
    plot([1+0.5*iplot, 1+0.5*iplot], [meandata, meandata+stddata],'k', [1+0.5*iplot, 1+0.5*iplot], [meandata, meandata-stddata], 'k');
    %plot horizontal error bars
    plot([0.9+0.5*iplot, 1.1+0.5*iplot], [meandata+stddata, meandata+stddata],'k', [0.9+0.5*iplot, 1.1+0.5*iplot], [meandata-stddata, meandata-stddata], 'k');
    xlim([0.7, 1.3+0.5*iplot]);
    
end
xticks(1:size(values,1));
ylabel(ylabel_str);
xlabel(xlabel_str);
title(title_str, 'Fontsize', 18); 
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');
set(gca,'Fontsize',15);
set(gca, 'XTickLabel',[]);



% 
% dmean = mean(data);                                                 % Mean
% dci  = std(data)*tinv(0.975,size(data,1)-1);                        % Confidence Intervals
% xt = 1;                                                         % X-Ticks
% xtd = repmat(xt, size(data,1), 1);                                  % X-Ticks For Data
% sb = [xt'-ones(size(data,2),1)*0.1,  xt'+ones(size(data,2),1)*0.1]; % Short Bar X
% lb = [xt'-ones(size(data,2),1)*0.2,  xt'+ones(size(data,2),1)*0.2]; % Long Bar X
% figure(1)
% plot(xt, data, '+')
% hold on
% for k1 = 1:size(data,2)
%     plot(lb(k1,:), [1 1]*dmean(k1), sb(k1,:), [1 1]*(dmean(k1)-dci(k1)), sb(k1,:), [1 1]*(dmean(k1)+dci(k1)), '-k')
% end
% hold off
% set(gca, 'XTick', xt, 'XTickLabel', {'Left','Middle','Right'})
% xlabel('Group')
% ylabel('Velocity (Furlongs/Fortnight)')