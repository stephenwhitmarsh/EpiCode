function preictal_project(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\projects\dtx\to_share'))
    addpath \\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/dtx/to_share'))
    addpath /network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/fieldtrip-20200607
end


ft_defaults

config = preictal_setparams;

%% create Spyking circus files

for ipatient = slurm_task_id%1 %[2 4] %1:6
    
    ipart = 1;
    
    %read muse markers
    MuseStruct = readMuseMarkers(config{ipatient}, true);
    
    %remove the end of the file, after the seizure to analyse (because we do
    %not use the post-ictal data, so better to remove it of the spike sorting)
    MuseStruct                  = addMuseBAD(config{ipatient},MuseStruct);
    
    %read header of the neuralynx concatenated file (assuming all channels have the same sampling freq and length)
    header = ft_read_header(fullfile(config{ipatient}.datasavedir,config{ipatient}.prefix(1:end-1),['p',num2str(ipart)],[config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs']));
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ipatient},true);
    
    %convert raw spikedata{ipart} into 1-trial spikedata{ipart}, so it can
    %be used with the same functions as SpikeTrials
    %add 2 '0' colomn to be consistent with startsample and endsample in trialinfos of other analysis
    %add '{1}' to be consistent with the other spiketrials structures in
    %other projects
    for ipart = 1:size(SpikeRaw, 1)
        filebegin                   = 0-header.nSamplesPre;
        fileend                     = header.nSamples-header.nSamplesPre;
        cfgtemp                     = [];
        cfgtemp.trl                 = [filebegin, fileend, 0, 0, 0, filebegin, fileend];
        cfgtemp.trlunit             = 'samples';
        cfgtemp.hdr                 = header;
        cfgtemp.timestampspersecond = header.TimeStampPerSample * header.Fs;
        SpikeTrials{ipart}{1}       = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
        SpikeTrials{ipart}{1}.hdr           = header;
        SpikeTrials{ipart}{1}.analysis_name = config{ipatient}.name{1};
    end
    
    %read spike waveforms
    %if you do not want spike waveform, switch comments of the 2 following lines
    if strcmp(config{irat}.compute_spikewaveforms, 'yes')
        SpikeWaveforms = readSpikeWaveforms(config{ipatient},SpikeTrials, false);
    else
        SpikeWaveforms = [];
    end
    
    
    %compute stats over time, for each unit
    stats = spikestatsOverTime(config{ipatient}, SpikeTrials, false);
    
    %compute stats after removing bursts (for regularity)
    cfgtemp = config{ipatient};
    cfgtemp.statstime.removebursts        = 'yes';
    cfgtemp.statstime.suffix              = '_withoutbursts';
    stats_without_bursts                  = spikestatsOverTime(cfgtemp, SpikeTrials, false);
    
    %find stats windows which intersect BAD Muse markers
    for ipart = 1:size(stats, 2)
        for ilabel = 1:size(stats{ipart},2)
            [~, bad_start_synctime] = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__START__');
            [~, bad_end_synctime]   = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__END__');
            if length(bad_start_synctime) ~= length(bad_end_synctime)
                error('Not the same amount of BAD start and end markers');
            end
            
            for i_unit = 1:size(SpikeTrials{ipart}{ilabel}.label, 2)
                for itrial = 1:size(SpikeTrials{ipart}{ilabel}.trialinfo, 1)
                    hasartefact = false(1,size(stats{ipart}{ilabel}.time{i_unit}{itrial},2));
                    for i_window = 1:size(stats{ipart}{ilabel}.time{i_unit}{itrial},2)
                        %go through each BAD interval
                        for iBAD = 1:size(bad_start_synctime,2)
                            %ignore bad markers if period is less than some length
                            if bad_end_synctime{1,iBAD} - bad_start_synctime{1,iBAD} < config{ipatient}.statstime.minbadtime
                                continue
                            end
                            if bad_start_synctime{1,iBAD} > stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,1) && bad_start_synctime{1,iBAD} > stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,2)
                                continue
                            elseif bad_end_synctime{1,iBAD} < stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,1) && bad_end_synctime{1,iBAD} < stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,2)
                                continue
                            else
                                hasartefact(i_window) = true;
                            end
                        end
                    end
                    
                    % remove stats windows which intersect bad markers
                    stats{ipart}{ilabel}.window_times{i_unit}{itrial}           = stats{ipart}{ilabel}.window_times{i_unit}{itrial}(~hasartefact,:);
                    stats_without_bursts{ipart}{ilabel}.window_times{i_unit}{itrial}    = stats_without_bursts{ipart}{ilabel}.window_times{i_unit}{itrial}(~hasartefact,:);
                    for ifield = ["time","isi_smooth","freq","amplitude","cv","fanofactor","burstindex","cv2"]
                        stats{ipart}{ilabel}.(ifield){i_unit}{itrial}                   = stats{ipart}{ilabel}.(ifield){i_unit}{itrial}(~hasartefact);
                        stats_without_bursts{ipart}{ilabel}.(ifield){i_unit}{itrial} 	= stats_without_bursts{ipart}{ilabel}.(ifield){i_unit}{itrial}(~hasartefact);
                    end
                    fprintf('Removed %d artefacted windows from %d\n', sum(hasartefact), size(hasartefact,2));
                end %itrial
            end %i_unit
        end %ilabel
    end %ipart
    
    %get timing of seizure, and convert from cell to mat
    [~, seizure_start]  = concatenateMuseMarkers(MuseStruct,ipart,'CriseStart'); %cfg.preictal.crisestart
    [~, seizure_end]    = concatenateMuseMarkers(MuseStruct,ipart,'CriseEnd'); %cfg.preictal.crisestart
    seizure_start       = cell2mat(seizure_start);
    seizure_end         = cell2mat(seizure_end);
    bad_start_synctime  = cell2mat(bad_start_synctime);
    bad_end_synctime    = cell2mat(bad_end_synctime);
    
    %plot results and compute stats on spike morpho ; and compute avg discharge values on all the data
        
    cfgtemp                             = config{ipatient};
    cfgtemp.spikewaveforms              = SpikeWaveforms;
    cfgtemp.statstime.plot.toi          = [0 seizure_end(1,end) + config{ipatient}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.marker1      = seizure_start(1,:);
    cfgtemp.statstime.plot.marker2      = seizure_end(1,:);
    cfgtemp.statstime.plot.bad_start    = bad_start_synctime(1,:);
    cfgtemp.statstime.plot.bad_end      = bad_end_synctime(1,:);
    cfgtemp.statstime.plot.suffix       = 'alldata';
    stats = plot_spikestats_preictal(cfgtemp,stats,SpikeTrials);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end) + config{ipatient}.statstime.plot.toi_seizure(1), seizure_end(1,end) + config{ipatient}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix       = 'seizurezoom';
    plot_spikestats_preictal(cfgtemp,stats,SpikeTrials);
    
    cfgtemp.statstime.plot.toi      = [0 seizure_end(1,end) + config{ipatient}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix   = 'alldata_withoutbursts';
    stats_without_bursts = plot_spikestats_preictal(cfgtemp,stats_without_bursts,SpikeTrials);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end) + config{ipatient}.statstime.plot.toi_seizure(1), seizure_start(1,end) + config{ipatient}.statstime.plot.toi_seizure(2)];
    cfgtemp.statstime.plot.suffix       = 'seizurezoom_withoutbursts';
    plot_spikestats_preictal(cfgtemp,stats_without_bursts,SpikeTrials);
        
    %save stats
    save(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime.mat']), 'stats', '-v7.3');
    save(fullfile(config{ipatient}.datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime_withoutbursts.mat']), 'stats_without_bursts', '-v7.3');
    
end %ipatient


end