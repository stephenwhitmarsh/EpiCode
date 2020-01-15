function [dat_micro_clean,dat_macro_clean,c] = removeartefacts_corr(cfg,dat_micro, dat_macro)

chanindx = 0;
for ichan = 1 : size(dat_micro.label,1)
    if strcmp(cfg.channel,dat_micro.label{ichan})
        chanindx = ichan;
    end
end
if chanindx == 0
    fprintf('Can not find channel %s! \n',cfg.channel);
end
datrms = zeros(size(dat_micro.trial,2),1);
for itrial = 1 : size(dat_micro.trial,2)
    datrms(itrial) = rms(dat_micro.trial{itrial}(chanindx,:));
end
fprintf('removing %d artefacted trials out of %d (> 3 STD(RMS)) \n',sum((datrms-mean(datrms)) >= std(datrms)*3),length(datrms));

cfgtemp = [];
cfgtemp.trials = find(datrms-mean(datrms) < std(datrms)*3);
dat_micro_clean = ft_selectdata(cfgtemp,dat_micro);
dat_macro_clean = ft_selectdata(cfgtemp,dat_macro);

cfgtemp = [];
cfgtemp.vartrllength = 2;
avg = ft_timelockanalysis(cfgtemp,dat_micro_clean);
clear c
for itrial = 1 : size(dat_micro_clean.trial,2)
    fprintf('correlating trial %d out of %d \n',itrial,size(dat_micro_clean.trial,2));
%     s = find(avg.time > -cfg.prestim,1,'first');
    s = size(dat_micro_clean.trial{itrial},2);
    c(itrial) = corr(avg.avg(chanindx,1:s)',dat_micro_clean.trial{itrial}(chanindx,1:s)');
end

cfgtemp = [];
cfgtemp.trials = find(c > 0);
fprintf('removing %d trials out of %d (c <= 0) \n',sum(c<=0),size(c,2));
dat_micro_clean = ft_selectdata(cfgtemp,dat_micro_clean);
dat_macro_clean = ft_selectdata(cfgtemp,dat_macro_clean);
c = c(c > 0);
