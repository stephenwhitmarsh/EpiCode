function preictal_project_SW

restoredefaultpath
if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\dtx\to_share'))
    addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/EpiCode/projects/dtx/to_share'))
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
end
ft_defaults

config = preictal_setparams;

for ielec = 1 : 9
    
    % read muse markers
    MuseStruct = readMuseMarkers(config{ielec}, true);
    
     
    % write artefacts to file
    writeSpykingCircusDeadfiles(config{ielec}, MuseStruct, true);
    
    %%
    % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresNotRemoved
    %%
    
    % add arteafct marker from seizure to end of file
    MuseStruct = addMuseBAD(config{ielec}, MuseStruct);
   
    %%
    % NOW RENAME SpykingCircus_artefacts_sample -> SpykingCircus_artefacts_samples_SeizuresRemoved
    %%
    
    % write parameters file for spyking circus
    writeSpykingCircusParameters(config{ielec});
    
    % write file list for spyking circus
    [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(config{ielec}, true);

    % write a new joblist for the cluster
    preictal_spikes_slurm_joblist
    
    %%
    % Now do your spike sorting
    %%
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ielec}, true);
    
    %read spike waveforms
    SpikeWaveforms = readSpikeWaveforms(config{ipatient}, SpikeRaw, true);

    
    %adapt spike structure to be consistent with spike trial structure in other projects
%     SpikeRaw{1}.(config{ielec}.name{1}) = SpikeRaw{1};
%     SpikeRaw{1} = keepfields(SpikeRaw{1}, (config{ielec}.name{1}));
    
    
    %compute stats over time, for each unit
    stats = spikestatsOverTime(config{ielec}, SpikeRaw, true);
    
    %compute stats after removing bursts (for regularity)
    cfgtemp = config{ielec};
    cfgtemp.statstime.removebursts        = 'yes';
    cfgtemp.statstime.suffix              = '_withoutbursts';
    stats_without_bursts                  = spikestatsOverTime(cfgtemp, SpikeRaw, true);
    
    %find stats windows which intersect BAD Muse markers
    % it would be normal to have lot of windows removed if the seizure occurs
    % early and you added BAD with addMuseBAD in all data after this seizure
    for ipart = 1:size(stats, 2)
        for markername = string(fieldnames(stats{ipart})')
            bad_start = concatenateMuseMarker(config{ielec},MuseStruct, ipart, 'BAD__START__');
            bad_end   = concatenateMuseMarker(config{ielec},MuseStruct, ipart, 'BAD__END__');
            if length(bad_start.synctime) ~= length(bad_end.synctime)
                error('Not the same amount of BAD start and end markers');
            end
            
            for i_unit = 1:size(SpikeRaw{ipart}.(markername).label, 2)
                for itrial = 1:size(SpikeRaw{ipart}.(markername).trialinfo, 1)
                    hasartefact = false(1,size(stats{ipart}.(markername).time{i_unit}{itrial},2));
                    for i_window = 1:size(stats{ipart}.(markername).time{i_unit}{itrial},2)
                        %go through each BAD interval
                        for iBAD = 1:size(bad_start.synctime,2)
                            %ignore bad markers if period is less than some length
                            if bad_end.synctime(iBAD) - bad_start.synctime(iBAD) < config{ielec}.statstime.minbadtime
                                continue
                            end
                            if bad_start.synctime(iBAD) > stats{ipart}.(markername).window_times{i_unit}{itrial}(i_window,1) && bad_start.synctime(iBAD) > stats{ipart}.(markername).window_times{i_unit}{itrial}(i_window,2)
                                continue
                            elseif bad_end.synctime(iBAD) < stats{ipart}.(markername).window_times{i_unit}{itrial}(i_window,1) && bad_end.synctime(iBAD) < stats{ipart}.(markername).window_times{i_unit}{itrial}(i_window,2)
                                continue
                            else
                                hasartefact(i_window) = true;
                            end
                        end
                    end
                    
                    % remove stats windows which intersect bad markers
                    stats{ipart}.(markername).window_times{i_unit}{itrial}           = stats{ipart}.(markername).window_times{i_unit}{itrial}(~hasartefact,:);
                    stats_without_bursts{ipart}.(markername).window_times{i_unit}{itrial}    = stats_without_bursts{ipart}.(markername).window_times{i_unit}{itrial}(~hasartefact,:);
                    for ifield = ["time","isi_smooth","freq","amplitude","cv","fanofactor","burstindex","cv2"]
                        stats{ipart}.(markername).(ifield){i_unit}{itrial}                   = stats{ipart}.(markername).(ifield){i_unit}{itrial}(~hasartefact);
                        stats_without_bursts{ipart}.(markername).(ifield){i_unit}{itrial} 	= stats_without_bursts{ipart}.(markername).(ifield){i_unit}{itrial}(~hasartefact);
                    end
                    fprintf('Removed %d artefacted windows from %d\n', sum(hasartefact), size(hasartefact,2));
                end %itrial
            end %i_unit
        end %ilabel
    end %ipart
    
    %get timing of seizure, and convert from cell to mat
    seizure_start  = concatenateMuseMarker(config{ielec},MuseStruct,ipart,'CriseStart'); %cfg.preictal.crisestart
    seizure_end    = concatenateMuseMarker(config{ielec},MuseStruct,ipart,'CriseEnd'); %cfg.preictal.crisestart
    
    %plot results and compute stats on spike morpho ; and compute avg discharge values on all the data
    
    cfgtemp                             = config{ielec};
    cfgtemp.spikewaveforms              = SpikeWaveforms;
    cfgtemp.statstime.plot.toi          = [0 seizure_end.synctime(end) + config{ielec}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.marker1      = seizure_start.synctime; %red
    cfgtemp.statstime.plot.marker2      = seizure_end.synctime; %blue
    cfgtemp.statstime.plot.bad_start    = bad_start.synctime;
    cfgtemp.statstime.plot.bad_end      = bad_end.synctime;
    cfgtemp.statstime.plot.suffix       = 'alldata';
    stats = preictal_plot_spikestats(cfgtemp,stats,SpikeRaw);
    
    cfgtemp.statstime.plot.toi = [seizure_start.synctime(end) + config{ielec}.statstime.plot.toi_seizure(1), seizure_end.synctime(end) + config{ielec}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix       = 'seizurezoom';
    preictal_plot_spikestats(cfgtemp,stats,SpikeRaw);
    
    cfgtemp.statstime.plot.toi      = [0 seizure_end.synctime(end) + config{ielec}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix   = 'alldata_withoutbursts';
    stats_without_bursts = preictal_plot_spikestats(cfgtemp,stats_without_bursts,SpikeRaw);
    
    cfgtemp.statstime.plot.toi = [seizure_start.synctime(end) + config{ielec}.statstime.plot.toi_seizure(1), seizure_start.synctime(end) + config{ielec}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix       = 'seizurezoom_withoutbursts';
    preictal_plot_spikestats(cfgtemp,stats_without_bursts,SpikeRaw);
    
    %save stats
    save(fullfile(config{ielec}.datasavedir,[config{ielec}.prefix, 'spikestatsOverTime.mat']), 'stats', '-v7.3');
    save(fullfile(config{ielec}.datasavedir,[config{ielec}.prefix, 'spikestatsOverTime_withoutbursts.mat']), 'stats_without_bursts', '-v7.3');
    
end