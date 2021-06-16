function hspike_cluster(ipatient)

%% Analysis script for SLURM cluster
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%

%% Add path

restoredefaultpath

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/hspike/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared/utilities
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/releaseDec2015/'));
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/subaxis
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external/sigstar-master
    addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/epishare-master'));
end

if ispc
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared\utilities
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\altmany-export_fig-8b0ba13\
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\subaxis
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\MatlabImportExport_v6.0.0
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external\sigstar-master
    addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\epishare-master'));
end

ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
labels = ["BAD__START__", "BAD__END__", "PHASE_1__START__", "PHASE_1__END__", "PHASE_2__START__", "PHASE_2__END__", "PHASE_3__START__", "PHASE_3__END__", "REM__START__", "REM__END__", "AWAKE__START__", "AWAKE__END__", "NO_SCORE__START__", "NO_SCORE__END__"];

%% General analyses
config                                                                  = hspike_setparams;
[MuseStruct_orig{ipatient}]                                             = readMuseMarkers(config{ipatient}, false);
[MuseStruct_aligned{ipatient}]                                          = alignMuseMarkersXcorr(config{ipatient}, MuseStruct_orig{ipatient}, false);
[~, LFP_cluster{ipatient}]                                              = clusterLFP(config{ipatient}, MuseStruct_aligned{ipatient}, false);
[MuseStruct_template{ipatient}, ~,~, ~]                                 = detectTemplate(config{ipatient}, MuseStruct_aligned{ipatient}, LFP_cluster{ipatient}{1}.Hspike.kmedoids{6}, false);
MuseStruct_template{ipatient}                                           = updateMarkers(config{ipatient}, MuseStruct_template{ipatient}, labels);

writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct_template{ipatient}, true);
writeSpykingCircusFileList(config{ipatient}, true);
writeSpykingCircusParameters_new(config{ipatient});



% %% t-zero LFPs
% fig = figure;
% MuseStruct_template_zeroed{ipatient} = MuseStruct_template{ipatient};
% for ipart = 1 : 3
%     for imarker = 1 : size(LFP_cluster_detected{ipatient}{1}, 2)
%         if isempty(LFP_cluster_detected{ipatient}{ipart}{imarker})
%             continue
%         end
%         chani       = find(strcmp(LFP_cluster_detected{ipatient}{ipart}{imarker}.label, config{ipatient}.align.zerochannel));
%         subplot(3,6,imarker + (ipart-1) * 6); hold;
%         plot(LFP_cluster_detected{ipatient}{ipart}{imarker}.time, LFP_cluster_detected{ipatient}{ipart}{imarker}.avg', 'k');
%         plot(LFP_cluster_detected{ipatient}{ipart}{imarker}.time, LFP_cluster_detected{ipatient}{ipart}{imarker}.avg(chani, :)', 'r');
%         timeindx    = LFP_cluster_detected{ipatient}{ipart}{imarker}.time > -0.1 & LFP_cluster_detected{ipatient}{ipart}{imarker}.time < 0.1;
%         [PKS, LOCS] = findpeaks(-LFP_cluster_detected{ipatient}{ipart}{imarker}.avg(chani, timeindx), 'SortStr', 'descend');
%         if isempty(LOCS)
%             continue
%         end
%         i0          = find(LFP_cluster_detected{ipatient}{ipart}{imarker}.time > -0.1, 1, 'first') + LOCS(1) - 1;
%         t0          = LFP_cluster_detected{ipatient}{ipart}{imarker}.time(i0);
%         axis tight; ax = axis;
%         plot([t0, t0], [1.2*ax(3), 1.2*ax(4)], ':r');
%         title(sprintf('Shift = %.0fms', -t0*1000), 'FontSize', 6);
%         for idir = 1 : size(MuseStruct_template_zeroed{ipatient}{ipart}, 2)
%             MuseStruct_template_zeroed{ipatient}{ipart}{idir}.markers.(['template', num2str(imarker)]).t0_shift = t0;
%             MuseStruct_template_zeroed{ipatient}{ipart}{idir}.markers.(['template', num2str(imarker)]).synctime = MuseStruct_template_zeroed{ipatient}{ipart}{idir}.markers.(['template', num2str(imarker)]).synctime + t0;
%             MuseStruct_template_zeroed{ipatient}{ipart}{idir}.markers.(['template', num2str(imarker)]).clock = MuseStruct_template_zeroed{ipatient}{ipart}{idir}.markers.(['template', num2str(imarker)]).clock + seconds(t0);
%         end
%     end
% end
% fname = fullfile(config{ipatient}.imagesavedir, 'alignment', [config{ipatient}.prefix, 'tzeroed']);
% disp(['Exporting figure ', fname])
% exportgraphics(fig, [fname, '.pdf']);
% exportgraphics(fig, [fname, '.tiff'], 'resolution', 600);
% 
% % update - add any new artefacts
% MuseStruct_template{ipatient}                                          = updateMarkers(config{ipatient}, MuseStruct_template_zeroed{ipatient}, labels);
% 
% % rename markers to combined markers (see config)
% MuseStruct_combined{ipatient}                                          = editMuseMarkers(config{ipatient}, MuseStruct_template{ipatient});
% 

% % from now on work on manual and combined templates
% itemp = 1;
% config{ipatient}.name = [];
% for markername = string(unique(config{ipatient}.editmarkerfile.torename(:, 2))')
%     config{ipatient}.name{itemp}                      = markername;
%     config{ipatient}.muse.startmarker.(markername)    = markername;
%     config{ipatient}.muse.endmarker.(markername)      = markername;
%     config{ipatient}.epoch.toi.(markername)           = [-0.5  1];
%     config{ipatient}.epoch.pad.(markername)           = 0.5;
%     config{ipatient}.LFP.baselinewindow.(markername)  = [-0.5  1];
%     config{ipatient}.LFP.baselinewindow.(markername)  = [-0.5  1];
%     config{ipatient}.LFP.name{itemp}                  = markername;
%     config{ipatient}.hyp.markers{itemp}               = markername;
%     itemp = itemp + 1;
% end


% [marker{ipatient}, hypnogram{ipatient}, hypmusestat{ipatient}] = hypnogramMuseStats(config{ipatient}, MuseStruct_combined{ipatient}, false);
% [LFP{ipatient}]                                                = readLFP(config{ipatient}, MuseStruct_combined{ipatient}, false);
% [TFR{ipatient}]                                                = TFRtrials(config{ipatient}, LFP{ipatient}, false);

% trim files to only those within a hypnogram
% MuseStruct_trimmed{ipatient}  = MuseStruct_combined{ipatient};
% config_trimmed{ipatient}      = config{ipatient};
% for ipart = 1 : 3
%     sel     = hypnogram{ipatient}.directory(hypnogram{ipatient}.part == ipart);
%     first   = find(strcmp(config{ipatient}.directorylist{ipart}, sel(1)));
%     last    = find(strcmp(config{ipatient}.directorylist{ipart}, sel(end)));
%     config_trimmed{ipatient}.directorylist{ipart}   = config{ipatient}.directorylist{ipart}(first:last);
%     MuseStruct_trimmed{ipatient}{ipart}             = MuseStruct_combined{ipatient}{ipart}(first:last);
% end
% 
% writeSpykingCircusDeadfiles(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, true);
% writeSpykingCircusFileList(config_trimmed{ipatient}, true);
% writeSpykingCircusParameters_new(config_trimmed{ipatient});

% 
% % read spike data from Phy as one continuous trial
% SpikeRaw{ipatient}                    = readSpikeRaw_Phy_new(config_trimmed{ipatient}, true);
% %
% % if ipatient == 3
% %     SpikeRaw{ipatient}{1}.template_maxchan(1) = 1;
% % end
% 
% % SpikeWaveforms_fromRaw{ipatient}              = readSpikeWaveforms_fromRaw(config_trimmed{ipatient}, SpikeRaw{ipatient}, true);
% 
% % segment into trials based on IED markers
% SpikeTrials_timelocked{ipatient}      = readSpikeTrials_MuseMarkers(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, true);
% SpikeDensity_timelocked{ipatient}     = spikeTrialDensity(config_trimmed{ipatient}, SpikeTrials_timelocked{ipatient}, true);
% 
% % segment into equal periods
% SpikeTrials_windowed{ipatient}        = readSpikeTrials_windowed(config_trimmed{ipatient}, MuseStruct_trimmed{ipatient}, SpikeRaw{ipatient}, true);
% SpikeStats_windowed{ipatient}         = spikeTrialStats(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);
% 
% ipart = 1;
% nrtrials = 10;
% for markername = string(config{ipatient}.LFP.name)
%     nrdirs = size(config{ipatient}.directorylist{ipart}, 2);
%     config{ipatient}.plot.trial.(markername) = repmat(1:nrtrials, 1, nrdirs);
%     config{ipatient}.plot.dir.(markername) = [];
%     for idir = 1 : nrdirs
%         config{ipatient}.plot.dir.(markername) = [config{ipatient}.plot.dir.(markername), ones(1, nrtrials) * idir];
%     end
% end
% 
% plot_patterns_multilevel_examples(config{ipatient});
% 
% % plot_patterns_multilevel(config{ipatient});
% 
% % SpikeWaveforms{ipatient}              = readSpikeWaveforms(config_trimmed{ipatient}, SpikeTrials_windowed{ipatient}, true);
