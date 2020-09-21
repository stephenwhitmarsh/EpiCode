%
%
% figure;hold;
%                         c_greys = 0.9 : -0.9 / size(cfg.SpikeTrials{ipart}{ievent}.trialinfo,1) : 0;
%
% %         plot(psth_event.time, psth_event.avg);
% %         plot(sdf1.time, sdf1.avg);
%         for itrial = 1:size(sdfdata1.trial,2)
%             plot(sdfdata1.time{itrial}, sdfdata1.trial{itrial}(i_unit,:), 'Color', [c_greys(itrial) c_greys(itrial) c_greys(itrial)]);
%         end
%
%             plot(psth_event.time, normalize(psth_event.avg,2,'range'));
%             plot(sdf1.time, sdf1.avg);
%             plot(sdf1.time, normalize(sdf1.avg,2,'range'));
%             plot(sdf1.time, log10(sdf1.avg));
%             plot(sdfdata1.time{itrial}, normalize(sdfdata1.trial{itrial},2,'range'));
%             plot(sdfdata1.time{itrial}, log10(sdfdata1.trial{itrial}));
%             xlim([-0.5 0.5]);
%
%             testxcorr = xcorr2(sdf1.avg);
%             figure;
%             plot(testxcorr);
for irat = 1:5
    sdf = stats{irat}{ipart}.SlowWave.sdfavg;
    %             figure;hold;
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        h{irat}(i_unit) = stats{irat}{ipart}.Interictal.maxchan{i_unit};
        [peak{irat}(i_unit), loc{irat}(i_unit)] = findpeaks(sdf.avg(i_unit,:), sdf.time,'NPeaks',1,'SortStr','descend'); %Npeaks : max nr of peaks/ SortStr : peak sorting : descend = from largest to smallest
        %                 plot(sdf.time, sdf.avg(i_unit,:), 'k');
        %                 scatter(loc(i_unit), peak(i_unit), 'xr');
    end
end
    
%plot profondeur en fonction du délai
figure;hold;
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'noise')
            if strcmp(stats{irat}{ipart}.Interictal.celltype{i_unit}, 'pn')
                plottype = '^b';
            else
                plottype = 'or';
            end
            if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'mua')
                scatter(loc{irat}(i_unit),h{irat}(i_unit),plottype, 'filled', 'MarkerFaceColor', plottype(2), 'MarkerEdgeColor', plottype(2));
            else
                scatter(loc{irat}(i_unit),h{irat}(i_unit),plottype, 'MarkerEdgeColor', plottype(2));
            end
        end
    end
end
setfig();
ylabel('deepness');
xlabel('time');
xlim([-0.2 0.2]);
yticklabels((yticks-16).*250);

%plot pic sdf en fonction du délai
figure;hold;
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'noise')
            if strcmp(stats{irat}{ipart}.Interictal.celltype{i_unit}, 'pn')
                plottype = '^b';
            else
                plottype = 'or';
            end
            if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'mua')
                scatter(loc{irat}(i_unit),peak{irat}(i_unit),plottype, 'filled', 'MarkerFaceColor', plottype(2), 'MarkerEdgeColor', plottype(2));
            else
                scatter(loc{irat}(i_unit),peak{irat}(i_unit),plottype, 'MarkerEdgeColor', plottype(2));
            end
        end
    end
end
setfig();
ylabel('maxfreq (Hz)');
xlabel('time');
xlim([-0.2 0.2]);
ylim([0 250]);

%plot du type en fonction du délai (voirproject spikes)
figure;hold;
for irat = 1:5
    for i_unit = 1:size(stats{irat}{ipart}.SlowWave.label,2)
        if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'noise')
            if strcmp(stats{irat}{ipart}.Interictal.celltype{i_unit}, 'pn')
                plottype = '^b';
            else
                plottype = 'or';
            end
            if  ~contains(stats{irat}{ipart}.Interictal.group{i_unit}, 'mua')
                scatter(loc{irat}(i_unit),stats{irat}{ipart}.SlowWave.code_slowwave{i_unit},plottype, 'filled', 'MarkerFaceColor', plottype(2), 'MarkerEdgeColor', plottype(2));
            else
                scatter(loc{irat}(i_unit),stats{irat}{ipart}.SlowWave.code_slowwave{i_unit},plottype, 'MarkerEdgeColor', plottype(2));
            end
        end
    end
end
setfig();
ylabel('type of behaviour');
xlabel('delay');
xlim([-0.2 0.2]);
ylim([0.5 5.5]);
