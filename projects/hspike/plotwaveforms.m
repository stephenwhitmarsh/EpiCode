figure;
ipart = 1;

for icluster = 1 : size(SpikeWaveforms{ipatient}{ipart}, 2)
    
    q = cat(1, SpikeWaveforms{ipatient}{ipart}{icluster}.trial{:});
    subplot(6,6,icluster); 
    plot(mean(q)); 
    
    labelW = SpikeWaveforms{ipatient}{ipart}{icluster}.label;
    labelRAW = SpikeRaw{ipatient}{ipart}.label{icluster};
    electrode = SpikeRaw{ipatient}{ipart}.channelname{icluster};
    maxchan = SpikeRaw{ipatient}{ipart}.template_maxchan(icluster);
    group = SpikeRaw{ipatient}{ipart}.cluster_group(icluster);
 
    title(strcat(num2str(icluster)," ", group," ",labelW," ",labelRAW," ",electrode, " ", num2str(maxchan)),'interpreter', 'none')
end
    
    
