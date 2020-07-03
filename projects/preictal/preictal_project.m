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
    
    %read spike data
    SpikeRaw = readSpikeRaw_Phy(config{ipatient},true);
    
    %read spike waveforms
    SpikeWaveforms = readSpikeWaveforms_raw(config{ipatient},SpikeRaw, false);
    
    %read header of the neuralynx concatenated file (assuming all channels have the same sampling freq and length)
    header = ft_read_header(fullfile(config{ipatient}.datasavedir,config{ipatient}.prefix(1:end-1),['p',num2str(ipart)],[config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs']));
    
    %compute stats over time, for each unit
    stats = spikestatsOverTime_raw(config{ipatient}, SpikeRaw, header, true);
    
    %compute stats after removing bursts (for regularity)
    cfgtemp = config{ipatient};
    cfgtemp.statstime.removebursts        = 'yes';
    stats_without_bursts = spikestatsOverTime_raw(config{ipatient}, SpikeRaw, header);
    
    %find stats windows which intersect BAD Muse markers
    for ipart = 1:size(stats, 1)
        
        [~, bad_start_synctime] = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__START__');
        [~, bad_end_synctime]   = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__END__');
        if length(bad_start_synctime) ~= length(bad_end_synctime)
            error('Not the same amount of BAD start and end markers');
        end
        
        for i_unit = 1:size(SpikeRaw{ipart}.label, 2)
            hasartefact = false(1,size(stats{ipart}.time{i_unit},2));
            ft_progress('init', 'text',     'Removing stats windows with artefacts');
            for i_window = 1:size(stats{ipart}.time{i_unit},2)
                ft_progress(i_window/size(stats{ipart}.time{i_unit},2), 'Checking artefacts on window %d from %d', i_window, size(stats{ipart}.time{i_unit},2));
                %go through each BAD interval
                for iBAD = 1:size(bad_start_synctime,2)
                    %ignore bad markers if periode is less than 500ms
                    if bad_end_synctime{1,iBAD} - bad_start_synctime{1,iBAD} < config{ipatient}.statstime.minbadtime
                        continue
                    end
                    if bad_start_synctime{1,iBAD} > stats{ipart}.window_times{i_unit}(i_window,1) && bad_start_synctime{1,iBAD} > stats{ipart}.window_times{i_unit}(i_window,2)
                        continue
                    elseif bad_end_synctime{1,iBAD} < stats{ipart}.window_times{i_unit}(i_window,1) && bad_end_synctime{1,iBAD} < stats{ipart}.window_times{i_unit}(i_window,2)
                        continue
                    else
                        hasartefact(i_window) = true;
                    end
                end
            end
            ft_progress('close');
            
            % remove stats windows which intersect bad markers
            stats{ipart}.window_times{i_unit}                   = stats{ipart}.window_times{i_unit}(~hasartefact,:);
            stats_without_bursts{ipart}.window_times{i_unit}    = stats_without_bursts{ipart}.window_times{i_unit}(~hasartefact,:);
            for ifield = ["time","isi_smooth","freq","amplitude","cv","fanofactor","burstindex","cv2"]
                stats{ipart}.(ifield){i_unit}                   = stats{ipart}.(ifield){i_unit}(~hasartefact);
                stats_without_bursts{ipart}.(ifield){i_unit} 	= stats_without_bursts{ipart}.(ifield){i_unit}(~hasartefact);
            end
            fprintf('Removed %d artefacted windows from %d\n', sum(hasartefact), size(hasartefact,2));
        end %i_unit
    end %ipart
        
    %get timing of seizure, and convertfrom cell to mat
    [~, seizure_start] = concatenateMuseMarkers(MuseStruct,ipart,'CriseStart'); %cfg.preictal.crisestart
    seizure_start = cell2mat(seizure_start);
    bad_start_synctime = cell2mat(bad_start_synctime);
    bad_end_synctime = cell2mat(bad_end_synctime);
    
    %plot results
    cfgtemp                             = config{ipatient};
    cfgtemp.spikewaveforms              = SpikeWaveforms;
    cfgtemp.statstime.plot.toi          = [0 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.marker       = seizure_start(1,:);
    cfgtemp.statstime.plot.bad_start    = bad_start_synctime(1,:);
    cfgtemp.statstime.plot.bad_end      = bad_end_synctime(1,:);
    stats = plot_spikestatsOverTime_raw(cfgtemp,stats,SpikeRaw);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end)-100 seizure_start(1,end)+100];
    plot_spikestatsOverTime_raw(cfgtemp,stats,SpikeRaw);
    
    cfgtemp.statstime.plot.toi      = [0 seizure_start(1,end)+100];
    cfgtemp.statstime.plot.suffix   = 'withoutbursts';
    stats_without_bursts = plot_spikestatsOverTime_raw(cfgtemp,stats_without_bursts,SpikeRaw);
    
    cfgtemp.statstime.plot.toi = [seizure_start(1,end)-100 seizure_start(1,end)+100];
    plot_spikestatsOverTime_raw(cfgtemp,stats_without_bursts,SpikeRaw);
    
    %save stats
    save(fullfile(datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime.mat']), stats, '-v7.3');
    save(fullfile(datasavedir,[config{ipatient}.prefix, 'spikestatsOverTime_withoutbursts.mat']), stats_without_bursts, '-v7.3');
    
end


end