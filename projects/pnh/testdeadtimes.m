
size(dat.trial{1}(1,:))

deadtimes

maxsamples = dat.fsample*26;
% maxsamples = size(dat.trial{1},2);

figure; hold;
plot(dat.time{1}(1:maxsamples),dat.trial{1}(1,1:maxsamples))

mindat = min(dat.trial{1}(1,1:maxsamples));
maxdat = max(dat.trial{1}(1,1:maxsamples));


idead = 1;
while deadtimes(idead,1) <= maxsamples && idead <= size(deadtimes,1)
    plot([deadtimes(idead,1)/dat.fsample,deadtimes(idead,1)/dat.fsample],[mindat,maxdat],'k:');
    plot([deadtimes(idead,2)/dat.fsample,deadtimes(idead,2)/dat.fsample],[mindat,maxdat],'k');
    idead = idead + 1;
end

ftemp = fullfile(cfg.datasavedir,[num2str(stimrate),'Hz_concatinated_',dat.label{1},'.ncs']);

cfgtemp = [];
cfgtemp.dataset = ftemp;
dattemp = ft_preprocessing(cfgtemp);


figure; 
subplot(2,1,1);
plot(dat.time{1}(1:maxsamples),dat.trial{1}(1,1:maxsamples))
subplot(2,1,2);
plot(dattemp.time{1}(1:maxsamples),dattemp.trial{1}(1,1:maxsamples),'r')

