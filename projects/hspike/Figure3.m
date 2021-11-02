
for ipatient = 1 : 8
    [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient});
end

clear totaldur totalsum
hypTable = table;
mrkTable = table;

for ipatient = 1 : size(hypnogram, 2)
    
    % plotting IEDs x Sleep stage for all marker (e.g. all templates)
    for ipart = 1:3
        for stage = unique(hypnogram{ipatient}.hyplabel)
            totaldur(ipatient, ipart, i+1) = seconds(sum(hypnogram{ipatient}.duration(hypnogram{ipatient}.stage == i & hypnogram{ipatient}.part == ipart)));
        end
        for i = 0 : 4
            totalsum(ipatient, ipart, i+1) = sum(marker{ipatient}.stage == i & marker{ipatient}.ipart == ipart);
        end
    end
%     hypnogram{ipatient}.ipatient = ipatient * ones(height(hypnogram{ipatient}),1);
%     marker{ipatient}.ipatient    = ipatient * ones(height(marker{ipatient}),1);
%     
    hypTable = [hypTable; hypnogram{ipatient}];
    temp = sortrows(marker{ipatient},'clock');
    temp.rate(2:height(temp)) = seconds(diff(temp.clock));
    mrkTable = [mrkTable; temp];
end

stgTable = table;
for ipatient = 1 : size(hypnogram, 2)
    
    % plotting IEDs x Sleep stage for all marker (e.g. all templates)
    for ipart = unique(hypnogram{ipatient}.part)'
        for i = 0 : 4
            temp            = table;
            temp.duration   = totaldur(ipatient, ipart, i+1);
            temp.number     = totalsum(ipatient, ipart, i+1);
            temp.rate       = totalsum(ipatient, ipart, i+1) / totaldur(ipatient, ipart, i+1);            
            temp.ipatient   = ipatient;
            temp.ipart      = ipart;
            temp.istage     = i;
            stgTable        = [stgTable; temp];
        end
    end
end
stgTable.istage = categorical(stgTable.istage);




filename = fullfile(config{ipatient}.datasavedir,'IEDs.csv');
writetable(mrkTable, filename)





figure;
boxplot(stgTable.rate,stgTable.istage)
lme = fitlme(stgTable,'rate~istage+(1|ipatient)+(1|ipart)');



mrkTableSel = mrkTable(mrkTable.rate > 0, :);
figure;
histogram(mrkTableSel.rate, 2000,'BinLimits',[0, 100]);

boxplot(mrkTableSel.rate, mrkTableSel.stage)





totaldur = round(squeeze(sum(totaldur, 2)), 0);
totalsum = round(squeeze(sum(totalsum, 2)), 0);
IEDrate  = totalsum ./ totaldur * 60;

totalsum_sum = round(sum(totalsum, 1), 0);
totaldur_sum = round(sum(totaldur, 1), 0);
IEDrate_sum  = round(sum(IEDrate, 1), 1);

format shortE
%% 
cfg.visible = 'on';
fig = figure('visible', cfg.visible);

% Number of IEDs per night / stage
subplot(5,2,1);
hb = bar(0:4, totalsum_sum, 1);
for ib = 1 : numel(hb)
    xData = hb(ib).XData + hb(ib).XOffset;
    text(xData, totalsum_sum(ib,:)', num2str(totalsum_sum(ib,:)'), 'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'});
ylabel('Total number of IEDs');
box off

% Number of IEDs per night / stage
subplot(5,2,2);
hb = bar(0:4, totalsum_sum, 1);
for ib = 1 : numel(hb)
    xData = hb(ib).XData + hb(ib).XOffset;
    text(xData, totalsum_sum(ib,:)', num2str(totalsum_sum(ib,:)'), 'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'});
ylabel('Total number of IEDs');
box off

% hours spend in sleep stage
subplot(5,2,3);
hb = bar(0 : 4, totaldur / 60 / 60, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,totaldur(ib,:)' /60 /60, num2str(round(totaldur(ib,:)' /60 /60, 1)),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
ylabel('Total duration (hrs)');
box off
% legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
ylim([0, max(max(totaldur) /60 /60) * 1.3]);
set(gca, 'TickLength', [0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% hours spend in sleep stage
subplot(5,2,4);
hb = bar(0 : 4, totaldur_sum / 60 / 60, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,totaldur_sum(ib,:)' /60 /60, num2str(round(totaldur_sum(ib,:)' /60 /60, 1)),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
ylabel('Total duration (hrs)');
box off
% legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
ylim([0, max(max(totaldur_avg) /60 /60) * 1.3]);
set(gca, 'TickLength', [0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% number of IEDs per minute
subplot(5,2,5);
hb = bar(0 : 4, IEDrate, 1); hold;
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData, IEDrate(ib,:)', num2str(round(IEDrate(ib,:),1)'), 'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'})
bar(0 : 4, IEDrate, 1);
ylabel('IEDs per minute');
box off

% number of IEDs per minute
subplot(5,2,6); hold;
hb = bar(0 : 4, IEDrate_sum, 1, 'stacked');
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData, IEDrate_sum(ib,:)', num2str(round(IEDrate_sum(ib,:),1)'), 'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'})
bar(0 : 4, IEDrate, 1, 'stacked');
ylabel('IEDs per minute');
box off

% legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% IED rate normalized to wake
subplot(5,2,7);
hb = bar(0 : 4, IEDrateNorm, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,IEDrateNorm(ib,:)', num2str(round(IEDrateNorm(ib,:),1)'),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
ylabel('Summed IED rate vs. wake');
box off
if max(max(IEDrateNorm)) > 0
    ylim([0, max(max(IEDrateNorm)) * 1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% IED rate normalized to wake
IEDrateNorm_sum = IEDrate_sum ./ IEDrate_sum(:,1); 
subplot(5,2,8);
hb = bar(0 : 4, IEDrateNorm_sum, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,IEDrateNorm_sum(ib,:)', num2str(round(IEDrateNorm_sum(ib,:),1)'),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
ylabel('Summed IED rate vs. wake');
box off
if max(max(IEDrateNorm_sum)) > 0
    ylim([0, max(max(IEDrateNorm_sum)) * 1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';


%
IEDrateNorm = IEDrate ./ IEDrate(:,1); 
subplot(5,2,10);
hb = bar(0 : 4, mean(IEDrateNorm), 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData, mean(IEDrateNorm), num2str(round(mean(IEDrateNorm),2)'),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
ylabel('Avg IED rate vs. wake');
box off
if max(max( mean(IEDrateNorm))) > 0
    ylim([0, max(max( mean(IEDrateNorm))) * 1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';


