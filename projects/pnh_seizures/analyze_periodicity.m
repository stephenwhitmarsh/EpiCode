imagesavedir = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/images';

dat = readtable('pattern_quantifications.xlsx');

% distance between discharges - absolute time
dat.interdischarge = [nan; diff(dat.start)];

% add counter for every discharge within pattern
nr = 0;
for i = 1 : find(~isnan(dat.diff),1,'last')
    if i > 1
        if dat.count(i) ~= dat.count(i-1)
            nr = 1;
            dat.interdischarge(i) = nan;
        else
            nr = nr + 1;
        end
        dat.nr(i) = nr;
    end 
end

% average distance between discharges per pattern
for inodule = unique(dat.nodule)'
    for icount = unique(dat.count(dat.nodule == inodule))'
        sel = dat.nodule == inodule & dat.count == icount;
        dat.avginterdischarge(sel) = nanmean(dat.interdischarge(sel));
    end
end

% boxplot duration discharges
fig = figure; 
boxplot(dat.diff*1000,dat.nodule)
title('Duration of discharges');
xlabel('Nodule');
ylabel('Duration discharge (ms)');
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(imagesavedir,'pattern_boxplot_duration_discharges.pdf'),'-r600');


% distance between discharges percentage of pattern
dat.interdischargeperc = ((dat.interdischarge ./ dat.avginterdischarge) * 100)-100;

avg_duration_discharge_1 = mean(dat.diff(dat.nodule == 1)*1000)
avg_duration_discharge_2 = mean(dat.diff(dat.nodule == 2)*1000)
avg_duration_discharge_3 = mean(dat.diff(dat.nodule == 3)*1000)

std_duration_discharge_1 = std(dat.diff(dat.nodule == 1)*1000)
std_duration_discharge_2 = std(dat.diff(dat.nodule == 2)*1000)
std_duration_discharge_3 = std(dat.diff(dat.nodule == 3)*1000)

avg_duration_interdischarge_1 = nanmean(dat.interdischarge(dat.nodule == 1)*1000)
avg_duration_interdischarge_2 = nanmean(dat.interdischarge(dat.nodule == 2)*1000)
avg_duration_interdischarge_3 = nanmean(dat.interdischarge(dat.nodule == 3)*1000)

std_duration_interdischarge_1 = nanstd(dat.interdischarge(dat.nodule == 1)*1000)
std_duration_interdischarge_2 = nanstd(dat.interdischarge(dat.nodule == 2)*1000)
std_duration_interdischarge_3 = nanstd(dat.interdischarge(dat.nodule == 3)*1000)

std_duration_interdischarge_perc_1 = nanstd(dat.interdischargeperc(dat.nodule == 1))
std_duration_interdischarge_perc_2 = nanstd(dat.interdischargeperc(dat.nodule == 2))
std_duration_interdischarge_perc_3 = nanstd(dat.interdischargeperc(dat.nodule == 3))


fig = figure;
title('Interval between discharges');
hold on
boxplot(dat.interdischarge*1000,dat.nodule);
xlabel('Nodule');
ylabel('Duration interval between discharge (ms)');

fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(imagesavedir,'pattern_interdischarge_interval_boxplot.pdf'),'-r600');


fig = figure;
% title('Interval between deflections in percentage per mean interval per pattern');
hold on
h1 = histogram(dat.interdischargeperc(dat.nodule == 1),30);
h2 = histogram(dat.interdischargeperc(dat.nodule == 2),30);
h3 = histogram(dat.interdischargeperc(dat.nodule == 3),30);
xlim([-100, 100]);
xlabel('% change in duration between deflections compared to mean interval');
ylabel('Nr. of discharges');
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
xticks([-100, -75, -50, -25, 0, 25, 50, 75, 100]);
ax = axis;
p1 = plot([-25, -25],[ax(3) ax(4)],'--k');
p2 = plot([25, 25],  [ax(3) ax(4)],'--k');
p3 = plot([50, 50],  [ax(3) ax(4)],':k');
p4 = plot([-50, -50],[ax(3) ax(4)],':k');

legend([h1, h2, h3, p1, p3],{'Nodule 1','Nodule 2','Nodule 3','25%','50%'});

print(fig, '-dpdf', fullfile(imagesavedir,'pattern_interval_histogram.pdf'),'-r600');



size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc))) <= 25),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc))) > 25),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc))) > 25 & abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc))) < 50),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc))) > 50 ),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 1 & ~isnan(dat.interdischargeperc)))),1)

size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc))) <= 25),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc))) > 25),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc))) > 25 & abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc)) < 50)),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc))) > 50),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 2 & ~isnan(dat.interdischargeperc)))),1)

size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc))) <= 25),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc))) > 25),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc))) > 25 & abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc))) < 50),1) / size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc)))),1)
size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc))) > 50),1)  / size(find(abs(dat.interdischargeperc(dat.nodule == 3 & ~isnan(dat.interdischargeperc)))),1)


% dat.nodule = categorical(dat.nodule,[1,2,3],{'one','two','three'});
% figure; gscatter(dat.diff,dat.count,dat.nodule,'bgr')
% xlabel('duration');
% ylabel('count');
fitted = fitlm(dat,'diff ~ count*nodule*nr','RobustOpts','on');
anova(fitted)

fitted = fitlme(dat,'diff ~ nodule + nr + (nodule|count) + (nr|count)');
anova(fitted)

fitted = fitlme(dat(~isnan(dat.interdischarge),:),'interdischarge ~ nodule + nr + (count|nodule)');
anova(fitted)


[h,p,ci,stats]  = ttest2(dat.diff(dat.nodule == 1),dat.diff(dat.nodule == 2)) %no
[h,p,ci,stats]  = ttest2(dat.diff(dat.nodule == 1),dat.diff(dat.nodule == 3)) %yes
[h,p,ci,stats]  = ttest2(dat.diff(dat.nodule == 2),dat.diff(dat.nodule == 3)) %yes

[h,p,ci,stats]  = ttest2(dat.interdischarge(dat.nodule == 1),dat.interdischarge(dat.nodule == 2)) %no
[h,p,ci,stats]  = ttest2(dat.interdischarge(dat.nodule == 1),dat.interdischarge(dat.nodule == 3)) %yes
[h,p,ci,stats]  = ttest2(dat.interdischarge(dat.nodule == 2),dat.interdischarge(dat.nodule == 3)) %yes


[h,p,ci,stats]  = ttest2(dat.interdischargeperc(dat.nodule == 1),dat.interdischargeperc(dat.nodule == 2)) %no
[h,p,ci,stats]  = ttest2(dat.interdischargeperc(dat.nodule == 1),dat.interdischargeperc(dat.nodule == 3)) %yes
[h,p,ci,stats]  = ttest2(dat.interdischargeperc(dat.nodule == 2),dat.interdischargeperc(dat.nodule == 3)) %yes


nanmean(dat.interdischarge(dat.nodule == 1)*1000)
nanmean(dat.interdischarge(dat.nodule == 2)*1000)
nanmean(dat.interdischarge(dat.nodule == 3)*1000)
nanstd(dat.interdischarge(dat.nodule == 1)*1000)
nanstd(dat.interdischarge(dat.nodule == 2)*1000)
nanstd(dat.interdischarge(dat.nodule == 3)*1000)


% 
% 
% 
% d = linspace(min(dat.diff),max(dat.diff));
% 
% figure; hold;
% gscatter(dat.diff,dat.count,dat.nodule,'bgr')
% line(d,feval(fitted,d,'one'),'Color','b');
% fit(d,'one')
% 
% 
% mdl = fitlm(dat,'diff ~ nodule*count*nr');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% 
% figure; boxplot(dat.diff,dat.nodule);
% mdl = fitlm(dat,'diff ~ nodule*count*nr');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% mdl = fitlm(dat,'diff ~ nodule');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% mdl = fitlm(dat(dat.nodule==1,:),'diff ~ count*nr');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% mdl = fitlm(dat(dat.nodule==2,:),'diff ~ count*nr');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% mdl = fitlm(dat(dat.nodule==3,:),'diff ~ count*nr');
% figure; plotEffects(mdl)
% anova(mdl,'summary')
% anova(mdl)
% 
% 
% 
% 
% nr = 0;
% for i = 1 : find(~isnan(dat.diff1),1,'last')
%     if i > 1
%         if dat.count1(i) > dat.count1(i-1)
%             nr = 1;
%         else
%             nr = nr + 1;
%         end
%     end
%     dat.nr1(i) = nr;
% end
% 
% nr = 0;
% for i = 1 : find(~isnan(dat.diff2),1,'last')
%     if i > 1
%         if dat.count2(i) > dat.count2(i-1)
%             nr = 1;
%         else
%             nr = nr + 1;
%         end
%     end
%     dat.nr2(i) = nr;
% end
% 
% nr = 0;
% for i = 1 : find(~isnan(dat.diff3),1,'last')
%     if i > 1
%         if dat.count3(i) > dat.count3(i-1)
%             nr = 1;
%         else
%             nr = nr + 1;
%         end
%     end
%     dat.nr3(i) = nr;
% end
% 
% clear pat
% pat = table
% ii = 1;
% for i = 1 : max(dat.count1)
%     pat.avg(ii) = mean(dat.diff1(dat.count1 == i));
%     pat.nr(ii)  = 1;
%     ii = ii + 1;
% end
% for i = 1 : max(dat.count2)
%     pat.avg(ii) = mean(dat.diff2(dat.count2 == i));
%     pat.nr(ii)  = 2;    
%     ii = ii + 1;
% end
% for i = 1 : max(dat.count3)
%     pat.avg(ii) = mean(dat.diff3(dat.count3 == i));
%     pat.nr(ii) = 3; 
%     ii = ii + 1;
% end
% 
% figure; boxplot(pat.avg,pat.nr);
% mdl = fitlm(pat,'avg ~ nr');
% 
% ttest(pat.avg
% 
% 


