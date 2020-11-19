SpikeWaveforms{ipart}.(markername){icluster}

itemp =1 ;

cm = parula(100);

figure; hold;
y = cat(1, SpikeWaveforms{ipart}.window{itemp}.trial{:});
avg_waveshape = mean(y, 1);
for itrial = 1 : size(SpikeWaveforms{ipart}.window{itemp}.trial, 2)
    lh = plot(SpikeWaveforms{ipart}.window{itemp}.time{itrial}*1000, SpikeWaveforms{ipart}.window{itemp}.trial{itrial}, 'color', [0.3 0.3 0.3]);
    lh.Color = [cm(itrial,:),1];
end
plot(SpikeWaveforms{ipart}.window{itemp}.time{itrial}*1000, avg_waveshape, 'color', [1 1 1 ], 'linewidth', 1);
axis tight
% ylim( [-max(abs(avg_waveshape)) max(abs(avg_waveshape))]*1.5);
xlabel('Time (ms)'); ylabel('uV');
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'box', 'off', 'TickDir', 'out');