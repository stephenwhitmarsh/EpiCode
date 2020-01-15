%% Analysis script
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/hspike/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/shared/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip/
ft_defaults

% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epilepsy\hspike\
% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epilepsy\shared\
% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip\
% ft_defaults

% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/mlib6/
% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/subaxis/

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
% maxNumCompThreads(4)


%% General analyses


for ipatient = 1
    
    config = hspike_setparams([]);
    
    % export hypnogram to muse
%     export_hypnogram(config{ipatient});

    % read muse markers
    [MuseStruct_micro, MuseStruct_macro] = readMuseMarkers_parts(config{ipatient}, false);
    
    % update paths for different OS
    [MuseStruct_micro, MuseStruct_macro] = MuseMarkers_update_filepath_parts(config{ipatient},MuseStruct_micro, MuseStruct_macro);

    % plot hypnogram
%     plotHypnogram(config{ipatient},MuseStruct_micro)

    
    % align data
    
    % read LFP data
%     [dat_micro, dat_macro] = readLFP_parts(config{ipatient}, MuseStruct_micro, MuseStruct_macro, false, false);
    
    % write data concatinated for SC, and update config with sampleinfo
    config{ipatient} = writeSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true, false); 
    
    % write parameters for spyking circus   
%     writeSpykingCircusParameters_parts(config{ipatient})
    
    % read spyking circus results
    [SpikeRaw, SpikeTrials] = readSpykingCircusSleepStage(config{ipatient}, MuseStruct_micro, true, 'all');

    % read and plot spikerate overview, and get the stats
%     [SpikeRateStats{ipatient}, stats_bar{ipatient}, sdf_orig_out{ipatient}, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, true, 1);
    
    % read and plot LFP of spike events
%     [spike_LFP]  = spikeLFP(config{ipatient},SpikeRaw, false);
    
    
end

    