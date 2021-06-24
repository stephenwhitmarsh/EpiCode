function T = analyze_periodicity(cfg, MuseStruct)

fname = fullfile(cfg{1}.datasavedir, 'Manual periodicity annotation.xlsx');
dat = readtable(fname);
fname_out = fullfile(cfg{1}.datasavedir, 'Periodicity.xlsx');

for ipatient = 1 : size(cfg, 2)
    intervals{ipatient} = inter_trial_intervals(cfg{ipatient}, MuseStruct{ipatient}, false);
end

% distance between discharges - absolute time
dat.interdischarge = [nan; diff(dat.Start)];

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

subplot(2,2,1);
boxplot(dat.diff*1000,dat.nodule)
title('Duration of discharges');
xlabel('Nodule');
ylabel('Duration discharge (ms)');
fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);


% distance between discharges percentage of pattern
dat.interdischargeperc = ((dat.interdischarge ./ dat.avginterdischarge) * 100) - 100;

for inodule = 1 : 4
    avg_duration_discharge{inodule} = mean(dat.diff(dat.nodule == inodule)*1000);
    std_duration_discharge{inodule} = std(dat.diff(dat.nodule == inodule)*1000);
    avg_duration_interdischarge{inodule} = nanmean(dat.interdischarge(dat.nodule == inodule)*1000);
    std_duration_interdischarge{inodule} = nanstd(dat.interdischarge(dat.nodule == inodule)*1000);
    std_duration_interdischarge{inodule} = nanstd(dat.interdischargeperc(dat.nodule == inodule));
    
    deviation{inodule}      = (dat.interdischarge(dat.nodule == inodule) ./ dat.avginterdischarge(dat.nodule == inodule) * 100) - 100;
    deviation_abs{inodule}  = abs(deviation{inodule});
    deviation{inodule}      = deviation{inodule}(~isnan(deviation{inodule} ));
    deviation_abs{inodule}  = deviation_abs{inodule}(~isnan(deviation_abs{inodule} ));
     
end


subplot(2,2,2);
title('Interval between discharges');
hold on
boxplot(dat.interdischarge, dat.nodule);
xlabel('Nodule');
ylabel('Duration interval between discharge (ms)');

fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig, '-dpdf', fullfile(cfg{1}.imagesavedir,'interval_boxplot.pdf'),'-r600');

subplot(2,2,3);
title('Interval between deflections in percentage per mean interval per pattern');
hold on
for inodule = 1 : 4
    h(inodule) = histogram(dat.interdischargeperc(dat.nodule == inodule), 'BinEdges', [-100:2.5:100]);
end
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

legend([h(:)', p1, p3],{'Nodule 1','Nodule 2','Nodule 3','Nodule 4','25%','50%'});

print(fig, '-dpdf', fullfile(cfg{1}.imagesavedir,'pattern_interval_histogram.pdf'),'-r600');

for inodule = 1 : 4
    fprintf('Nodule %d: %0.1f%% (n=%d/%d) of the waves varied less than 25%% from the mean interval between waves (%0.0f ms).\n', inodule, mean(deviation_abs{inodule} <= 25)*100, sum(deviation_abs{inodule} <= 25), size(deviation_abs{inodule}, 1), avg_duration_interdischarge{inodule});
end
% summary of numbers and mode ISI of each pattern

t = table;
clear totaltime
for ipatient = 1 : 4
    totaltime(ipatient) = seconds(0);
    for idir = 1 : size(MuseStruct{ipatient}{1}, 2)
        totaltime(ipatient) = totaltime(ipatient) + MuseStruct{ipatient}{1}{idir}.endtime - MuseStruct{ipatient}{1}{idir}.starttime;
    end
    t.nodule(ipatient) = ipatient;    
    t.totaltime(ipatient) = hours(totaltime(ipatient));
end

for ipatient = 1 : 4
    for markername = ["PSW", "FA", "ES"]
        if isfield(intervals{ipatient}.table, markername)
            fprintf('%s, Nodule %d: n=%d, λ=%0.0fms\n', markername, ipatient, size(intervals{ipatient}.table.(markername), 1), intervals{ipatient}.mode.(markername)*1000);

            fn = sprintf('%s_n', markername);
            t.(fn)(ipatient) = size(intervals{ipatient}.table.(markername), 1);
           
            fn = sprintf('%s_λ (ms)', markername);
            t.(fn)(ipatient) = intervals{ipatient}.mode.(markername)*1000;
        end
    end
end

for ipatient = 1 : 4
    fprintf('Nodule %d: %0.1f%% (n=%d/%d) of the waves varied less than 25%% from the mean interval between waves (%0.0f ms).\n', ipatient, mean(deviation_abs{ipatient} <= 25)*100, sum(deviation_abs{ipatient} <= 25), size(deviation_abs{ipatient}, 1), avg_duration_interdischarge{ipatient});
    t.mean_inter_wave_interval(ipatient)    = avg_duration_interdischarge{ipatient};
    t.std_inter_wave_interval(ipatient)     = std_duration_interdischarge{ipatient};
%     t.sel_n_dev(ipatient)                   = sum(deviation_abs{ipatient} <= 25);
%     t.sel_n(ipatient)                       = size(deviation_abs{ipatient}, 1);
%     t.perc_dev_from_mean(ipatient)          = mean(deviation_abs{ipatient} <= 25)*100;
    t.("proportion deviation")(ipatient) = sprintf("%d/%d (%0.1f%%)", sum(deviation_abs{ipatient} <= 25), size(deviation_abs{ipatient}, 1), mean(deviation_abs{ipatient} <= 25)*100);
end

writetable(t, fname_out);

