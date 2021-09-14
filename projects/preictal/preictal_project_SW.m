function preictal_project_SW

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))   
    
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
end
ft_defaults

config = preictal_setparams_SW;

%% prepare data for Spyking Circus
% write a new joblist for the cluster
preictal_spikes_slurm_joblist

for ielec = 3 % à définir
    
    ipart = 1; 

    % read muse markers
    MuseStruct = readMuseMarkers(config{ielec}, false);
    
    % remove post ictal from the whole analysis,
    % according to config (some 'patients' will have
    % shorter postictal kept because of noise, see setparams)
    MuseStruct = addMuseBAD(config{ielec}, MuseStruct);
    
    % write artefacts to file
    writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresNotRemoved');
    
    %%
    % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresNotRemoved
    %%
    
    % add arteafct marker from seizure start to seizure end  
    cfgtemp                       = [];
    cfgtemp.bad.markerStart       = 'CriseStart';
    cfgtemp.bad.markerEnd         = 'CriseEnd';
    MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);
    
    writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true, '_SeizuresRemoved');
   
    %%
    % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresRemoved
    %%
    
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
    SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, false);
    
    % add (sliding) timewindow
    [config{ielec}, MuseStruct]     = addSlidingWindows(config{ielec}, MuseStruct);
    
    % epoch data into windows
    SpikeTrials                    = readSpikeTrials(config{ielec}, MuseStruct, SpikeRaw, false);
    
    % calculate statistics per window
    SpikeStats                    = spikeTrialStats(config{ielec}, SpikeTrials, false);
    
    % plot overview 
        
    
end

% plotting
for ipatient  = 2 : 3
    
    % read settings
    config = preictal_setparams_SW;

    % read muse markers
    MuseStruct = readMuseMarkers(config{ipatient}, false);
    
    % add (sliding) timewindow
    [config{ipatient}, MuseStruct] = addSlidingWindows(config{ipatient}, MuseStruct);
  
%     % template LFP
    config{ipatient}.LFP.name   = {'window'};    
    LFP               = readLFP(config{ipatient}, MuseStruct, true);
    
    % calculate FFT on sliding timewindow
    config{ipatient}.FFT.name  = {'window'};
    FFT              = FFTtrials(config{ipatient}, true);
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ipatient}, false);
    
    % epoch spike data into windows
    config{ipatient}.spike.name  = {'window'};    
    SpikeTrials                    = readSpikeTrials(config{ipatient}, MuseStruct, SpikeRaw, false);
    
    % calculate statistics per window
    SpikeStats                    = spikeTrialStats(config{ipatient}, SpikeTrials, false);
    
    % you can plot any metric shown here: SpikeStats{1}.window{1}
    iunit = 1:size(SpikeStats{ipart}.window, 2); % spike to plot
    tunit = 3; % for spike distance
    ipart = 1;

    cfg                 = [];
    cfg.ipart           = ipart;
    cfg.filename        = fullfile(config{ipatient}.imagesavedir, sprintf('windowed_s%sp%d.jpg', config{ipatient}.prefix, ipart));
    cfg.xlim            = seconds([MuseStruct{1}{1}.starttime, MuseStruct{1}{2}.endtime] - MuseStruct{1}{2}.markers.CriseStart.clock);
    cfg.orientation     = 'landscape';
    cfg.offset          = MuseStruct{1}{2}.markers.CriseStart.clock;
    
    cfg.type{1}         = 'power';
    cfg.frequency{1}    = [1, 7];
    cfg.channel{1}      = 'all';
    cfg.title{1}        = sprintf('Power %d-%dHz', cfg.frequency{1}(1), cfg.frequency{1}(2));
    cfg.plotart(1)      = false;
    cfg.log(1)          = true;
    cfg.hideart(1)      = false;
    cfg.marker_indx{1}  = FFT{ipart}.window.trialinfo.BAD_cnt>0;
    cfg.marker_sign{1}  = '.r';
    cfg.marker_label{1} = 'artefact';
    
    cfg.type{2}         = 'power';
    cfg.frequency{2}    = [8, 14];
    cfg.channel{2}      = 'all';
    cfg.title{2}        = sprintf('Power %d-%dHz', cfg.frequency{2}(1), cfg.frequency{2}(2));
    cfg.plotart(2)      = false;
    cfg.log(2)          = true;
    cfg.hideart(2)      = true;
    
    cfg.type{3}         = 'relpower';
    cfg.frequency1{3}   = [1, 7];
    cfg.frequency2{3}   = [8, 14];
    cfg.channel{3}      = 'all';
    cfg.title{3}        = sprintf('Power (%d-%dHz)/(%d-%dHz)', cfg.frequency1{3}(1), cfg.frequency1{3}(2), cfg.frequency2{3}(1), cfg.frequency2{3}(2));
    cfg.plotart(3)      = true;
    cfg.log(3)          = false;
    
    cfg.type{4}        = 'spike';
    cfg.title{4}       = sprintf('CV2 unit %d-%d', iunit(1), iunit(end));
    cfg.metric{4}      = 'CV2_trial';
    cfg.unit{4}        = iunit;
    cfg.plotart(4)     = true;
    cfg.log(4)         = false;
    
    cfg.type{5}         = 'trialinfo';
    cfg.title{5}        = sprintf('BAD duration');
    cfg.metric{5}       = 'BAD_sec';
    
    cfg.type{6}         = 'spike';
    cfg.title{6}        = sprintf('Firing rate (corrected) unit %d-%d', iunit(1), iunit(end));
    cfg.metric{6}       = 'trialfreq_corrected';
    cfg.unit{6}         = iunit;
    cfg.plotart(6)      = true;
    cfg.log(6)          = true;
    
    cfg.type{7}         = 'spike';
    cfg.title{7}        = sprintf('Nr. bursts unit %d-%d', iunit(1), iunit(end));
    cfg.metric{7}       = 'burst_trialsum';
    cfg.unit{7}         = iunit;
    cfg.plotart(7)      = true;
    
    cfg.type{8}         = 'spike';
    cfg.title{8}        = sprintf('CV unit %d-%d', iunit(1), iunit(end));
    cfg.metric{8}       = 'CV_trial';
    cfg.unit{8}         = iunit;
    cfg.plotart(8)      = true;

    cfg.type{9}         = 'spike';
    cfg.title{9}        = sprintf('Distance %d to %d-%d', iunit(1), iunit(1), iunit(end));
    cfg.metric{9}       = 'dist';
    cfg.unit{9}         = iunit(1);
    cfg.index{9}        = 1:size(iunit,2)-1;
    cfg.plotart(9)      = true;
    
    plotWindowedData(cfg, MuseStruct, SpikeTrials, SpikeStats, FFT);
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    cfg                 = [];
    cfg.offset          = MuseStruct{1}{2}.markers.CriseStart.clock;
    cfg.ipart           = 1;
    cfg.filename        = '\\lexport\iss01.charpier\analyses\vn_preictal\images\test.jpg';
    
    cfg.type{1}         = 'spike';
    cfg.title{1}        = sprintf('firing rate unit %d', iunit);
    cfg.metric{1}       = 'trialfreq';
    cfg.unit(1)         = iunit;
    cfg.plotart(1)      = true;
    cfg.marker_indx{1}  = SpikeTrials{1}.window.trialinfo.BAD_sec > 5;
    cfg.marker_sign{1}  = '.b';
    
    cfg.type{2}         = 'spike';
    cfg.title{2}        = sprintf('firing rate (corrected) unit %d', iunit);
    cfg.metric{2}       = 'trialfreq_corrected';
    cfg.unit(2)         = iunit;
    cfg.plotart(2)      = true;
    cfg.marker_indx{2}  = SpikeTrials{1}.window.trialinfo.BAD_sec > 5;
    cfg.marker_sign{2}  = 'dk';
    
    cfg.type{3}         = 'spike';
    cfg.title{3}        = sprintf('Bursts unit %d', iunit);
    cfg.metric{3}       = 'burst_trialsum';
    cfg.unit(3)         = iunit;
    cfg.plotart(3)      = false;

    cfg.type{4}         = 'spike';
    cfg.title{4}        = sprintf('Spike distance between unit %d and %d', iunit, tunit);
    cfg.metric{4}       = 'dist';
    cfg.unit(4)         = iunit;
    cfg.index(4)        = tunit;
    cfg.plotart(4)      = false;
    cfg.marker_indx{4}  = SpikeTrials{1}.window.trialinfo.BAD_sec > 5;
    cfg.marker_sign{4}  = '.b';
    plotWindowedData(cfg, MuseStruct, SpikeStats, SpikeTrials);
    
end
