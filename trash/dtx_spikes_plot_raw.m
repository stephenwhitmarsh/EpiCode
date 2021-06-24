
config = dtx_spikes_setparams;
ipart = 1;
irat = 1;

%% plot example of raw intercritic
SpikeTrials   = readSpikeTrials_MuseMarkers(config{irat}, [], [], false);
%same trial length as the intracellular example : 
trialsel      = find((SpikeTrials{ipart}.Interictal.trialtime(:,2) - SpikeTrials{ipart}.Interictal.trialtime(:,1)) > 170 & (SpikeTrials{ipart}.Interictal.trialtime(:,2) - SpikeTrials{ipart}.Interictal.trialtime(:,1)) < 220);

itrial = trialsel(3);
idir = SpikeTrials{ipart}.Interictal.trialinfo.idir(itrial);

channame = 'E14';%E07, E10, E11, E12
temp = dir(fullfile(config{irat}.rawdir, config{irat}.directorylist{ipart}{idir}, sprintf('*%s.ncs', channame)));
datapath = fullfile(temp.folder, temp.name);

cfgtemp             = [];
cfgtemp.trl(1)      = SpikeTrials{ipart}.Interictal.trialinfo.begsample(itrial) - SpikeTrials{ipart}.Interictal.trialinfo.fileoffset(itrial);
cfgtemp.trl(2)      = SpikeTrials{ipart}.Interictal.trialinfo.endsample(itrial) - SpikeTrials{ipart}.Interictal.trialinfo.fileoffset(itrial);
cfgtemp.trl(3)      = SpikeTrials{ipart}.Interictal.trialinfo.offset(itrial);
cfgtemp.dataset     = datapath;
cfgtemp.hpfilter    = 'yes';
cfgtemp.hpfreq      = 300;
cfgtemp.hpfiltord   = 3;
dat                 = ft_preprocessing(cfgtemp);

%remove artefacts
Fs = SpikeTrials{ipart}.Interictal.hdr.Fs;
[~,loc] = findpeaks(-dat.trial{1}, 'MinPeakHeight', 200);
for iart = 1:size(loc, 2)
    toremove = loc(iart)- round(0.02*Fs) : loc(iart) + round(0.02*Fs);
    dat.trial{1}(toremove) = nan;
end

fig = figure; hold on;
plot(dat.time{1}, dat.trial{1}, 'k', 'linewidth', 1);

set(gca, 'tickdir', 'out', 'fontsize', 25);

unitssel = find(contains(SpikeTrials{ipart}.Interictal.label, channame));
 
spikes=[];
for i_unit = unitssel(2)
     idx = SpikeTrials{ipart}.Interictal.trial{i_unit} == itrial;
     spikes = SpikeTrials{ipart}.Interictal.time{i_unit}(idx);
     for ispike = 1 : size(spikes, 2)
         sel = spikes(ispike) - round(0.02 * Fs) : spikes(ispike) + round(0.02 * Fs);
         sel = sel - SpikeTrials{ipart}.Interictal.trialinfo.begsample(itrial);
         plot(dat.time{1}(sel), dat.trial{1}(sel), 'color', [1 0.6 0.2]);
     end
 end

set(gca, 'tickdir', 'out', 'fontsize', 25);

figname = fullfile(config{irat}.imagesavedir, '..', 'raw_spikes', sprintf('%sIntercritic_trial%d', config{irat}.prefix, itrial));
dtx_savefigure(fig, figname, 'pdf', 'png', 'close');
    

%% PLOT RAW SLOWWAVE :  EEG and MUA example : same seizure for each group

irat    = 1;
SpikeTrials_timelocked{irat} = readSpikeTrials_MuseMarkers(config{irat}, [], [], false);


for itrial = 98%1:size(SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo, 1)
    idir    = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.idir(itrial);
    
    %load surface eeg
    temp = dir(fullfile(config{irat}.rawdir, config{irat}.directorylist{ipart}{idir}, '*E12.ncs'));
    cfgtemp         = [];
    cfgtemp.trl(1)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.begsample(itrial) - SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.fileoffset(itrial);
    cfgtemp.trl(2)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.endsample(itrial) - SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.fileoffset(itrial);
    cfgtemp.trl(3)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.offset(itrial);
    cfgtemp.dataset = fullfile(temp.folder, temp.name);
    cfgtemp.lpfilter = 'yes';
    cfgtemp.lpfreq   = 100;
    eeg             = ft_preprocessing(cfgtemp);
    
    fig = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(sprintf('%d', itrial));
    subplot(4,1,1)
    plot(eeg.time{1}, -eeg.trial{1}, 'k', 'LineWidth', 2);
    set(gca, 'TickDir', 'out', 'FontSize', 25, 'LineWidth', 2);
    ylim([-4000 1500])
    iplot = 1;
	
    for unitname = ["cluster_7_E14","cluster_0_E08", "cluster_0_E10"]
               
        %load raw MUA
        chan = split(unitname, '_');
        chan = chan(end);
        temp = dir(fullfile(config{irat}.rawdir, config{irat}.directorylist{ipart}{idir}, sprintf('*%s.ncs', chan)));
        cfgtemp         = [];
        cfgtemp.trl(1)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.begsample(itrial) - SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.fileoffset(itrial);
        cfgtemp.trl(2)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.endsample(itrial) - SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.fileoffset(itrial);
        cfgtemp.trl(3)  = SpikeTrials_timelocked{irat}{ipart}.SlowWave.trialinfo.offset(itrial);
        cfgtemp.dataset = fullfile(temp.folder, temp.name);
        cfgtemp.hpfilter = 'yes';
        cfgtemp.hpfilttype  = 'but';
        cfgtemp.hpfiltord   = 3;
        cfgtemp.hpfreq      = 300;
        mua             = ft_preprocessing(cfgtemp);
        
        i_unit = find(strcmp(unitname, SpikeTrials_timelocked{irat}{ipart}.SlowWave.label));
        iplot = iplot+1;
        subplot(4,1,iplot); hold on;
        plot(mua.time{1}, mua.trial{1}, 'k', 'linewidth', 2);

             for ispike = find(SpikeTrials_timelocked{irat}{ipart}.SlowWave.trial{i_unit} == itrial)
                 t = SpikeTrials_timelocked{irat}{ipart}.SlowWave.time{i_unit}(ispike);
                 %plot([t t], [-1 1], 'k', 'LineWidth', 2);
                  sel = mua.time{1}> t-0.001 & mua.time{1} < t+0.001;
                  plot(mua.time{1}(sel), mua.trial{1}(sel), 'r');
             end
		
        xlim([-2 2]);
        set(gca, 'TickDir', 'out', 'FontSize', 25, 'LineWidth', 2);
    end
    figname = fullfile(config{irat}.imagesavedir, '..', 'raw_spikes', sprintf('%s_slowwave_raw_onetrial-%d', config{irat}.prefix, itrial));
    dtx_savefigure(fig, figname, 'png', 'pdf', 'close');
    
end

