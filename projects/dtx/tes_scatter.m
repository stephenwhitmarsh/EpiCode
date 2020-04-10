% for patients : one subplot per patient, and one figure for all patients
% for each parameter


figure;
hold;

for i = 0:1
    data = Seizure_Infos{2}{1}.timeBetween2SlowWaves;
    data = seconds(data);
    meandata = mean(data);
    stddata = std(data);
    xdata = (rand(1,size(data,2)))*0.2+0.9+0.5*i; %random values between 0.9 and 1.1
    scatter(xdata,data);
    plot([0.8+0.5*i, 1.2+0.5*i],[meandata, meandata], 'k');
    plot([1+0.5*i, 1+0.5*i], [meandata, meandata+stddata],'k', [1+0.5*i, 1+0.5*i], [meandata, meandata-stddata], 'k');
    plot([0.9+0.5*i, 1.1+0.5*i], [meandata+stddata, meandata+stddata],'k', [0.9+0.5*i, 1.1+0.5*i], [meandata-stddata, meandata-stddata], 'k');
    xlim([0.7, 1.3+0.5*i]);
    xticks([]);
end

ylabel('Time between 2 SlowWaves (seconds)');
xlabel('DTX5');
title('Time between 2 SlowWaves : 

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