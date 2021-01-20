itrial = cfg.representtrials(1)


temp = rmfield(TFR_macro_single,'time');
temp.time{1} = TFR_macro_single.time;

cfg = [];
cfg.time = temp.time;
cfg.trials = itrial;
cfg.method = 'pchip';
dat_micro_resampled = ft_resampledata(cfg,dat_micro);

A = squeeze(TFR_macro_single.powspctrm(1,:,:))';
B = squeeze(dat_micro_resampled.trial{1}(2,:)');


q = corr(A,B,'rows','pairwise');
figure; plot(TFR_macro_single.freq,q);