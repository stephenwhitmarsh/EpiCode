fig = figure;

subplot(1,2,1);

hold;
h = 250;
n = 1;
for itrial = 1:53
    i1 = find(dat_micro_chansel.time{itrial} >= -cfg.prestim,1,'first');
    i2 = find(dat_micro_chansel.time{itrial} >= cfg.poststim,1,'first');
    plot(dat_micro_chansel.time{itrial}(i1:i2),dat_micro_chansel.trial{itrial}(i1:i2) + n*h,'color','k');
    n = n + 1;
end
ylabel('Trials');
xlabel('Time (s)');
axis tight
set(gca, 'YTickLabel', '');

subplot(1,2,2);

% plot grid of trial voltages (micro)
imagesc(trialgrid_micro*255);
set(gca,'YDir','normal')
% colormap hot(255)
set(gca, 'XTickLabel', '');
% title('Raw trial amplitudes over time');
ylabel('Trials');
axis tight
ax = axis;

      % print to file
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(cfg.imagesavedir,['alltrials_overview_timecourses_.pdf']),'-r600');

