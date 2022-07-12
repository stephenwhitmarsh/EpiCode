function preictal_project_LC

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))
    %     addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_Dec_2019'))

    addpath \\lexport\iss01.charpier\analyses\vn_preictal\scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
    %     addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_Dec_2019'))
    addpath /network/lustre/iss01/charpier/analyses/vn_preictal/scripts/fieldtrip-20200607
end
ft_defaults

config = preictal_setparams;

%% prepare data for Spyking Circus
% write a new joblist for the cluster
preictal_spikes_slurm_joblist

for ielec = 12 % à définir
    
    ipart = 1; 
    
    % read muse markers
    MuseStruct = readMuseMarkers(config{ielec}, true);
    
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
    
    % write file list for spyking circus ---> LINUX
    writeSpykingCircusFileList(config{ielec}, true);

    %%
    % Now do your spike sorting
    %%
    
    %% perform the analysis after spike sorting 
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ielec}, true);
    
    %read spike waveforms
%     SpikeWaveforms = readSpikeWaveforms(config{ielec}, SpikeRaw, false);
    config{ielec}.spikewin.windowoverlap = 0;

    % create sliding timewindows
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
    config{ielec}.muse.startmarker.window   = 'window__START__';
    config{ielec}.muse.endmarker.window     = 'window__END__';
    config{ielec}.spike.toi.window          = [0 0];
    config{ielec}.spike.pad.window          = 0;
    config{ielec}.spike.name                = "window";
    config{ielec}.spike.postfix             = '-windowed';
    
    % epoch data into windows
    SpikeTrials_windowed                    = readSpikeTrials_MuseMarkers(config{ielec}, MuseStruct, SpikeRaw, true);
    
    % calculate statistics per window
    SpikeStats_windowed                     = spikeTrialStats(config{ielec}, SpikeTrials_windowed, true);
    
    % plot overview
    config{ielec}.plot.artthresh            = 1; % threshold artefacted based on total time artefacted per window (in seconds)
    config{ielec}.imagesavedir              = '\\lexport\iss01.charpier\analyses\vn_preictal\images'; % 
    plotOverviewPreictal(config{ielec}, MuseStruct, SpikeTrials_windowed, SpikeStats_windowed)
    
end
