function preictal_project_SW

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))
%     addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_Dec_2019'))
    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
%     addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_Dec_2019'))
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
end
ft_defaults

config = preictal_setparams;

%% prepare data for Spyking Circus
% write a new joblist for the cluster
preictal_spikes_slurm_joblist

for ielec = 3 % à définir
    
    % read muse markers
    MuseStruct = readMuseMarkers(config{ielec}, false);
    
    % remove post ictal from the whole analysis,
    % according to config (some 'patients' will have
    % shorter postictal kept because of noise, see setparams)
    MuseStruct = addMuseBAD(config{ielec}, MuseStruct);
    
    % write artefacts to file
    writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresNotRemoved');
    
%     %%
%     % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresNotRemoved
%     %%
    
    % add arteafct marker from seizure start to seizure end  
    cfgtemp                       = [];
    cfgtemp.bad.markerStart       = 'CriseStart';
    cfgtemp.bad.markerEnd         = 'CriseEnd';
    MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);
    
    writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresRemoved');
   
%     %%
%     % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresRemoved
%     %%
    
    % write parameters file for spyking circus 
    writeSpykingCircusParameters(config{ielec});
    
    % write file list for spyking circus LINUX
    writeSpykingCircusFileList(config{ielec}, true);
    
    %%
    % Now do your spike sorting
    %%
    
    %% perform the analysis after spike sorting 
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ielec}, false);
    
    %read spike waveforms
%     SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, false);
  
    % create sliding timewindows
    % overwrite settings
    for ipart = 1 : size(config{ielec}.directorylist, 2)
        for idir = 1 : size(config{ielec}.directorylist{ipart}, 2)
            temp = dir(fullfile(config{ielec}.rawdir, config{ielec}.directorylist{ipart}{idir}, ['*', config{ielec}.circus.channel{1}, '.ncs']));
            hdr  = ft_read_header(fullfile(config{ielec}.rawdir, config{ielec}.directorylist{ipart}{idir}, temp.name));
            MuseStruct{ipart}{idir}.markers.window__START__.synctime = 1 : (config{ielec}.spikewin.windowsize - config{ielec}.spikewin.windowsize * config{ielec}.spikewin.windowoverlap) : hdr.nSamples/hdr.Fs - (config{ielec}.spikewin.windowsize);
            MuseStruct{ipart}{idir}.markers.window__START__.clock    = seconds(MuseStruct{ipart}{idir}.markers.window__START__.synctime) + MuseStruct{ipart}{idir}.starttime;
            MuseStruct{ipart}{idir}.markers.window__START__.events   = size(MuseStruct{ipart}{idir}.markers.window__START__.synctime, 2);
            MuseStruct{ipart}{idir}.markers.window__END__.synctime   = MuseStruct{ipart}{idir}.markers.window__START__.synctime + config{ielec}.spikewin.windowsize;
            MuseStruct{ipart}{idir}.markers.window__END__.clock      = MuseStruct{ipart}{idir}.markers.window__START__.clock + seconds(config{ielec}.spikewin.windowsize);
            MuseStruct{ipart}{idir}.markers.window__END__.events     = size(MuseStruct{ipart}{idir}.markers.window__END__.synctime, 2);  
        end
    end
    config{ielec}.muse.startmarker.window    = 'window__START__';
    config{ielec}.muse.endmarker.window      = 'window__END__';
    config{ielec}.spike.toi.window           = [0 0];
    config{ielec}.spike.pad.window           = 0;
    config{ielec}.spike.name                 = "window";
    config{ielec}.spike.postfix              = '-windowed';
    
    SpikeTrials_windowed                    = readSpikeTrials_MuseMarkers(config{ielec}, MuseStruct, SpikeRaw, false);
    
    SpikeStats_windowed                     = spikeTrialStats(config{ielec}, SpikeTrials_windowed, false);
  
    Spiky_windowed = compute_synchrony_spiky(config{ielec}, SpikeTrials_windowed, false);

    % create time-axis relative to seizures
    time = seconds(SpikeTrials_windowed{1}.window.trialinfo.starttime - MuseStruct{1}{2}.markers.CriseStart.clock);
    
    for iunit = 1 : size(SpikeStats_windowed{1}.window, 2)
        
        unit_t = 1:size(SpikeStats_windowed{1}.window, 2);
        unit_t(iunit) = [];
        
        clear spike_dist
        for iwindow = 1 : size(Spiky_windowed{1}.window.trials, 2)
            spike_dist(iwindow, :) = Spiky_windowed{1}.window.trials{iwindow}.SPIKE.matrix(iunit, unit_t);
        end
        
        clear CV2_trial FR

        % plot selection
        figure
        subplot(4,1,1);
        plot(time, SpikeStats_windowed{1}.window{iunit}.trialfreq); title(sprintf('FR unit %d %s', iunit, deblank(SpikeTrials_windowed{1}.window.cluster_group{iunit})));                  xlim([time(1), time(end)]);  ylim([0, 10]);
        
        subplot(4,1,2); hold;
        plot([time(1), time(end)], [1 1],':k');
        plot(time, SpikeStats_windowed{1}.window{iunit}.CV2_trial); title('CV2');          xlim([time(1), time(end)]);  ylim([.2, 1.8]);
        
        subplot(4,1,3);
        plot(time, SpikeStats_windowed{1}.window{iunit}.burst_trialsum); title('Nr. of Bursts');   xlim([time(1), time(end)]);  %ylim([.2, 1.8]);
        
        subplot(4,1,4);
        plot(time, spike_dist); title('Spike distance');  xlim([time(1), time(end)]);  ylim([0, 1]);
        
    end
end
