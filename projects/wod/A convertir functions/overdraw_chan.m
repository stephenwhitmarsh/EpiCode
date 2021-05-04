
fig = figure;hold on
color = 'k';%color of trials when plotted
h = config{ipatient}.plotseizure.h /2;%5000;% FIXME à ajuster

%plot each trial
for itrial = 1:size(data_plot.trialinfo,1)
    ichan = 0;
    for channame = string(data_plot.label)'
        
        ichan = ichan+1;
        chan_idx = strcmp(channame,(data_plot.label));
        
        p = plot(data_1chan.time{itrial},data_1chan.trial{itrial}(chan_idx,:)+(numel(data_plot.label)+1)*h-h*ichan,'Color', color); %first on top
        p.Color(4) = 0.2;
        
        %plot average after the last trial
        if itrial == size(data_plot.trialinfo,1)
            chan_idx = strcmp(channame, data_plot_avg.label);
            plot(data_plot_avg.time, data_plot_avg.avg(chan_idx,:)+(numel(data_plot.label)+1)*h-h*ichan, 'b', 'LineWidth', 2);
        end
        
    end
end


%set figure display
axis tight
ylim([-h (numel(data_plot.label)+2)*h]);
xlim([-2 2]);
xlabel('Time (s)', 'Fontsize',15);
ylabel('Channel name', 'Fontsize',15);
tick = h;
yticks(h : tick : numel(data_plot.label)*h);
chan_name_all = string(data_plot.label)';
set(gca, 'YTickLabel', flip(chan_name_all));
set(gca, 'FontWeight','bold', 'Fontsize',15, 'TickDir','out', 'YTickLabel',flip(chan_name_all));


