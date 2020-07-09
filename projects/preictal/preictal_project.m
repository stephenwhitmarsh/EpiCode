function preictal_project(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\projects\dtx'))
    addpath \\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/dtx'))
    addpath /network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/fieldtrip-20200607
end


ft_defaults

config = preictal_setparams;

%% create Spyking circus files

for ipatient = slurm_task_id%1 %[2 4] %1:6
    
    %read muse markers
    MuseStruct = readMuseMarkers(config{ipatient}, true);
    
    %remove the end of the file, after the seizure to analyse (because we do
    %not use the post-ictal data, so better to remove it of the spike sorting)
    cfgtemp                     = [];
    cfgtemp.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
    cfgtemp.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
    cfgtemp.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers
    cfgtemp.bad.sample_list     = ft_getopt(config{ipatient},'seizure_index', 'last'); %index of the seizure to analyze, on the LAST dir. can be 'last' (default)
    cfgtemp.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
    MuseStruct                  = addMuseBAD(cfgtemp,MuseStruct);
    
    %read header of the neuralynx concatenated file (assuming all channels have the same sampling freq and length)
    header = ft_read_header(fullfile(config{ipatient}.datasavedir,config{ipatient}.prefix(1:end-1),['p',num2str(ipart)],[config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs']));
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ipatient},true);
    
    %convert raw spikedata{ipart} into 1-trial spikedata{ipart}, so it can
    %be used with the same functions as SpikeTrials
    %add 2 '0' colomn to be consistent with startsample and endsample in trialinfos of other analysis
    %add '{1}' to be consistent with the other spiketrials structures in
    %other projects
    if strcmp(config{ipatient}.statstime.latency, 'all')
        config{ipatient}.statstime.latency = [0-header.nSamplesPre header.nSamples] ./ header.Fs;
    end
    for ipart = 1:size(SpikeRaw, 1)
        cfgtemp                     = [];
        cfgtemp.trl                 = [config{ipatient}.statstime.latency, 0, 0, 0, config{ipatient}.statstime.latency] .* header.Fs;
        cfgtemp.trlunit             = 'samples';
        cfgtemp.hdr                 = header;
        cfgtemp.timestampspersecond = header.TimeStampPerSample * header.Fs;
        SpikeTrials{ipart}{1}       = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
        SpikeTrials{ipart}{1}.hdr           = header;
        SpikeTrials{ipart}{1}.analysis_name = config{ipatient}.name{1};
    end
    
    %read spike waveforms
    SpikeWaveforms = readSpikeWaveforms(config{ipatient},SpikeTrials, false);
    
    %compute stats over time, for each unit
    stats = spikestatsOverTime(config{ipatient}, SpikeTrials, header, true);
    
    %compute stats after removing bursts (for regularity)
    cfgtemp = config{ipatient};
    cfgtemp.statstime.removebursts        = 'yes';
    stats_without_bursts = spikestatsOverTime(config{ipatient}, SpikeTrials, header);
    
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
                    for i_window = 1:size(stats{ipart}{ilabel}.time{i_unit}{itrial}{itrial},2)
                        %go through each BAD interval
                        for iBAD = 1:size(bad_start_synctime,2)
                            %ignore bad markers if period is less than some length
                            if bad_end_synctime{1,iBAD} - bad_start_synctime{1,iBAD} < config{ipatient}.statstime.minbadtime
                                continue
                            end
                            if bad_start_synctime{1,iBAD} > stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}(i_window,1) && bad_start_synctime{1,iBAD} > stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}(i_window,2)
                                continue
                            elseif bad_end_synctime{1,iBAD} < stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}(i_window,1) && bad_end_synctime{1,iBAD} < stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}(i_window,2)
                                continue
                            else
                                hasartefact(i_window) = true;
                            end
                        end
                    end
                    
                    % remove stats windows which intersect bad markers
                    stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}           = stats{ipart}{ilabel}.window_times{i_unit}{itrial}{itrial}(~hasartefact,:);
                    stats_without_bursts{ipart}.window_times{i_unit}{itrial}{itrial}    = stats_without_bursts{ipart}.window_times{i_unit}{itrial}{itrial}(~hasartefact,:);
                    for ifield = ["time","isi_smooth","freq","amplitude","cv","fanofactor","burstindex","cv2"]
                        stats{ipart}{ilabel}.(ifield){i_unit}{itrial}{itrial}                   = stats{ipart}{ilabel}.(ifield){i_unit}{itrial}{itrial}(~hasartefact);
                        stats_without_bursts{ipart}.(ifield){i_unit}{itrial}{itrial} 	= stats_without_bursts{ipart}.(ifield){i_unit}{itrial}{itrial}(~hasartefact);
                    end
                    fprintf('Removed %d artefacted windows from %d\n', sum(hasartefact), size(hasartefact,2));
                end %itrial
            end %i_unit
        end %ilabel
    end %ipart
    
    %get timing of seizure, and convert from cell to mat
    [~, seizure_start] = concatenateMuseMarkers(MuseStruct,ipart,'CriseStart'); %cfg.preictal.crisestart
    seizure_start = cell2mat(seizure_start);
    bad_start_synctime = cell2mat(bad_start_synctime);
    bad_end_synctime = cell2mat(bad_end_synctime);
    
    %plot results and compute stats on spike morpho ; and compute avg discharge values on all the data
    cfgtemp                             = config{ipatient};
    cfgtemp.spikewaveforms              = SpikeWaveforms;
    cfgtemp.statstime.plot.toi          = [0 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.marker       = seizure_start(1,:);
    cfgtemp.statstime.plot.bad_start    = bad_start_synctime(1,:);
    cfgtemp.statstime.plot.bad_end      = bad_end_synctime(1,:);
    cfgtemp.statstime.plot.suffix       = '_alldata';
    stats = plot_spikestats_preictal(cfgtemp,stats,SpikeTrials);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end)-100 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.suffix       = '_seizurezoom';
    plot_spikestats_preictal(cfgtemp,stats,SpikeTrials);
    
    cfgtemp.statstime.plot.toi      = [0 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.suffix   = '_alldata_withoutbursts';
    stats_without_bursts = plot_spikestats_preictal(cfgtemp,stats_without_bursts,SpikeTrials);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end)-100 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.suffix       = '_seizurezoom_withoutbursts';
    plot_spikestats_preictal(cfgtemp,stats_without_bursts,SpikeTrials);
    
    %save stats
    save(fullfile(datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime.mat']), stats, '-v7.3');
    save(fullfile(datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime_withoutbursts.mat']), stats_without_bursts, '-v7.3');
    
end %ipatient


end