cfgtemp                 = [];
cfgtemp.spikechannel    = SpikeTrials{ipart}{baseline_index}.label{i_unit};
cfgtemp.normtime        = 'no';
cfgtemp.normvalues      = 'no';
cfgtemp.cutlength       = 20; % if normtime = 'yes', must be between 0 and 1. Otherwise, it is in seconds
cfgtemp.method          = 'freq';
cfgtemp.removebursts    = 'no';
cfgtemp.timelock        = 'no';
cfgtemp.removeempty     = 'no';
cfgtemp.removeoutlier   = 'no';
cfgtemp.plot            = 'trialavg';
cfgtemp.color           = [];
cfgtemp.saveplot        = 'no';
cfgtemp.name            = [];
cfgtemp.prefix          = [];
cfgtemp.imagesavedir    = [];

[stats{ipart}.(cfg.name{baseline_index}).stats_over_time.freq{i_unit}, legend] = spikestatsOverTime(cfgtemp,SpikeTrials{ipart}{baseline_index});

axis tight
ax = axis;
% xticks(0:3600:ax(2));
% xticklabels(xticks/3600); %convert seconds to hours
% xlim([0 Inf]);
% xlabel('Time (hours)');
% set(gca, 'YColor', 'k');
% ylabel(sprintf('Mean firing rate \nof each trial'));
% yticklabels(10.^yticks);
setfig();