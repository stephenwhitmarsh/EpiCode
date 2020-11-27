function hspike_cluster(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
 

%% Add path

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/development/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/fieldtrip/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/DBSCAN/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/altmany-export_fig-8b0ba13   
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\development
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\fieldtrip\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\DBSCAN\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0 % to read neuralynx files faster
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses
config                                                                  = hspike_setparams;
[MuseStruct_orig{ipatient}]                                             = readMuseMarkers(config{ipatient}, false);
[MuseStruct_aligned{ipatient}]                                          = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
[clusterindx{ipatient}, LFP_cluster{ipatient}]                          = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
[MuseStruct_template{ipatient}, ~,~, LFP_cluster_detected{ipatient}]    = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);

% update to any new artefacts
[MuseStruct_orig{ipatient}]                                             = readMuseMarkers(config{ipatient}, true);
for ipart = 1 : 3
    for idir = 1 : size(MuseStruct_template{ipatient}{ipart}, 2)
        try
            MuseStruct_template{ipatient}{ipart}{idir}.markers.BAD__START__ = MuseStruct_orig{ipatient}{ipart}{idir}.markers.BAD__START__;
            MuseStruct_template{ipatient}{ipart}{idir}.markers.BAD__END__ = MuseStruct_orig{ipatient}{ipart}{idir}.markers.BAD__END__;
        catch
        end
    end
end

switch ipatient
    case 1
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined1'; 'template5', 'combined2'; 'template6', 'combined2'};
    case 2
        markernames = {'combined1'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template4', 'combined1'; 'template6', 'combined1'};
    case 3
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined2'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined1'};
    case 4
        markernames = {'combined1', 'combined2', 'combined3'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined2'; 'template4', 'combined2'; 'template5', 'combined2'; 'template6', 'combined3'};
    case 5
        markernames = {'combined1', 'combined2', 'combined3'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined2'; 'template3', 'combined3'; 'template5', 'combined2'; };
    case 6
        markernames = {'combined1'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined1'; 'template4', 'combined1'; 'template5', 'combined1'; 'template6', 'combined1'};
    case 7
        markernames = {'combined1', 'combined2'};
        config{ipatient}.editmarkerfile.torename = {'template1', 'combined1'; 'template2', 'combined1'; 'template3', 'combined2'; 'template4', 'combined2'; 'template5', 'combined1'; 'template6', 'combined2'};
end

MuseStruct_combined{ipatient} = editMuseMarkers(config{ipatient}, MuseStruct_template{ipatient});

% focus time period a bit more
config{ipatient}.epoch.toi.Hspike           = [-0.2  0.8];
config{ipatient}.epoch.pad.Hspike           = 0.5;
config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];
config{ipatient}.LFP.baselinewindow.Hspike  = [-0.2  0.8];

% from now on work on manual and combined templates
itemp = 1;
for markername = string(markernames)
    config{ipatient}.name{itemp}                      = markername;
    config{ipatient}.muse.startmarker.(markername)    = markername;
    config{ipatient}.muse.endmarker.(markername)      = markername;
    config{ipatient}.epoch.toi.(markername)           = [-0.2  0.8];
    config{ipatient}.epoch.pad.(markername)           = 0.5;
    config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
    config{ipatient}.LFP.baselinewindow.(markername)  = [-0.2  0.8];
    config{ipatient}.LFP.name{itemp}                  = markername;
    config{ipatient}.hyp.markers{itemp}               = markername;
    itemp = itemp + 1;
end

[marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_combined{ipatient}, false);
[LFP{ipatient}]                                                = readLFP(config{ipatient}, MuseStruct_combined{ipatient}, false);
[TFR{ipatient}]                                                = TFRtrials(config{ipatient}, LFP{ipatient}, false);

% trim files to only those within a hypnogram
MuseStruct_trimmed  = MuseStruct_combined;
config_trimmed      = config;
for ipart = 1 : 3
    sel     = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
    first   = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
    last    = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
    config_trimmed{ipatient}.directorylist{ipart}   = config{ipatient}.directorylist{ipart}(first:last);
    MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct_combined{ipatient}{ipart}(first:last);

    % if still more than 7, cut off the beginning
    if size(config_trimmed{ipatient}.directorylist{ipart}, 2) > 7
        config_trimmed{ipatient}.directorylist{ipart}   = config_trimmed{ipatient}.directorylist{ipart}(end-6:end);
        MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct_trimmed{ipatient}{ipart}(end-6:end);
    end
end

% read spike data from Phy as one continuous trial
SpikeRaw{ipatient}                    = readSpikeRaw_Phy(config_trimmed{ipatient}, true);

% segment into trials based on IED markers
SpikeTrials_timelocked{ipatient}      = readSpikeTrials_MuseMarkers(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, false);
SpikeDensity_timelocked{ipatient}     = spikeTrialDensity(config_trimmed{ipatient}, SpikeTrials_timelocked{ipatient}, false);

% segment into equal periods
SpikeTrials_windowed{ipatient}        = readSpikeTrials_windowed(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, false);
SpikeStats_windowed{ipatient}         = spikeTrialStats(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, false, 'windowed');
SpikeWaveforms{ipatient}              = readSpikeWaveforms(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);

% plotOverviewHspike(config{ipatient}, marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}, ...
%     SpikeTrials_timelocked{ipatient}, SpikeTrials_windowed{ipatient}, SpikeStats_windowed{ipatient}, ...
%     SpikeDensity_timelocked{ipatient}, LFP{ipatient}, TFR{ipatient}, SpikeWaveforms{ipatient});

plotOverviewHspike(config{ipatient}, marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}, ...
    SpikeTrials_timelocked{ipatient}, SpikeTrials_windowed{ipatient}, SpikeStats_windowed{ipatient}, ...
    SpikeDensity_timelocked{ipatient}, LFP{ipatient}, TFR{ipatient}, SpikeWaveforms{ipatient});
