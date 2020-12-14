
ipatient = 3;
ipart = 1;

datafile{1} = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2660\p1\mHa2g\2660-p1-multifile-mHa2g_3.ncs';
datafile{2} = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2660\p1\mHa2g\2660-p1-multifile-mHa2g_4.ncs';
datafile{3} = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2660\p1\mHa2g\2660-p1-multifile-mHa2g_7.ncs';
datafile{4} = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2660\p1\mHa2g\2660-p1-multifile-mHa2g_8.ncs';

for ichan = 1 : 4
    cfgtemp                         = [];
    cfgtemp.dataset                 = datafile{ichan};
    dat_chan{ichan}                 = ft_preprocessing(cfgtemp);
end


load '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data\hspike\2660-SpikeRaw_Phy.mat'

icluster = 6;

spikes_idx_sel = sort(randperm(length(SpikeRaw{ipart}.sample{icluster}), 1000));
hdr = SpikeRaw{ipart}.hdr;

% define Fieldtrip trials
trialcount  = 0;
Startsample = [];
Endsample   = [];
Offset      = [];
Trialnr     = [];

% note that the sample field has to been added by ft_spike_maketrials, so use the adapted version in /external/fieldtrip:
% note that the indexing of the channels can be different.
for ispike = spikes_idx_sel
    trialcount   = trialcount + 1;
    %                     Startsample  = [Startsample; round(SpikeRaw{ipart}.sample{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(1) * hdr.Fs)];
    %                     Endsample    = [Endsample;   round(SpikeRaw{ipart}.sample{icluster}(ispike) / hdr.TimeStampPerSample  + cfg.spikewaveform.toi(2) * hdr.Fs)];
    Startsample  = [Startsample; double(round(SpikeRaw{ipart}.sample{icluster}(ispike) - 0.001 * hdr.Fs))];
    Endsample    = [Endsample;   double(round(SpikeRaw{ipart}.sample{icluster}(ispike) + 0.001 * hdr.Fs))];
    Offset       = [Offset; -round(0.01 * hdr.Fs)];
    Trialnr      = [Trialnr; trialcount];
end



cfgtemp                         = [];
cfgtemp.trl                     = [Startsample, Endsample, Offset];
cfgtemp.trlunit                 = 'samples';
trials                          = ft_redefinetrial(cfgtemp, dat);
        
figure; plot(trials.time{1}, mean(cat(1, trials.trial{:})));
