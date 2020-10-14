%% Analysis script DT PROJECT
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/shared/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/dtx/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip/
addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/spikes-master'));
addpath(genpath('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/npy-matlab-master'));
ft_defaults

% 
% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epilepsy\shared
% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\epilepsy\dtx
% addpath \\lexport\iss01.charpier\analyses\stephen.whitmarsh\fieldtrip
% addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\spikes-master'));
% addpath(genpath('\\lexport\iss01.charpier\analyses\stephen.whitmarsh\scripts\npy-matlab-master'));

% ft_defaults

% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/mlib6/
% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/subaxis/

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
% maxNumCompThreads(4)


%% General analyses

for ipatient = 1:1
    
    config = dtx_setparams([]);    
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
    
    % align Muse markers according to peaks and detect whether they contain artefacts
    [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);   
    [MuseStruct_micro, MuseStruct_macro]    = MuseMarkers_update_filepath(config{ipatient},MuseStruct_micro, MuseStruct_macro);
    
    
    % read LFP data
%    [dat_micro, dat_macro] = readLFP(config{ipatient}, MuseStruct_micro, MuseStruct_macro, false, false);
    
    % create 'artefacts' to remove seizures and time from seizures
%     for ipart = 1 : size(MuseStruct_micro,2)
%         
%         % only when seizures are present
%         if isfield(MuseStruct_micro{ipart}.markers,'Crise_Start')
%             if isfield(MuseStruct_micro{ipart}.markers.Crise_Start,'offset')
%                 
%                 % if no artefact fields exist, create empty ones
%                 if ~isfield(MuseStruct_micro{ipart}.markers,'BAD__START__')
%                     MuseStruct_micro{ipart}.markers.BAD__START__.offset      = [];
%                     MuseStruct_micro{ipart}.markers.BAD__START__.synctime    = [];
%                     MuseStruct_micro{ipart}.markers.BAD__START__.clock       = [];
%                     
%                     MuseStruct_micro{ipart}.markers.BAD__END__.offset        = [];
%                     MuseStruct_micro{ipart}.markers.BAD__END__.synctime      = [];
%                     MuseStruct_micro{ipart}.markers.BAD__END__.clock         = [];
%                 end
%                 
%                 MuseStruct_micro{ipart}.markers.BAD__START__.offset      = [MuseStruct_micro{ipart}.markers.BAD__START__.offset,   MuseStruct_micro{ipart}.markers.Crise_Start.offset];
%                 MuseStruct_micro{ipart}.markers.BAD__START__.synctime    = [MuseStruct_micro{ipart}.markers.BAD__START__.synctime, MuseStruct_micro{ipart}.markers.Crise_Start.synctime];
%                 MuseStruct_micro{ipart}.markers.BAD__START__.clock       = [MuseStruct_micro{ipart}.markers.BAD__START__.clock,    MuseStruct_micro{ipart}.markers.Crise_Start.clock];
%                 
%                 MuseStruct_micro{ipart}.markers.BAD__END__.offset        = [MuseStruct_micro{ipart}.markers.BAD__END__.offset,     MuseStruct_micro{ipart}.markers.Crise_End.offset];
%                 MuseStruct_micro{ipart}.markers.BAD__END__.synctime      = [MuseStruct_micro{ipart}.markers.BAD__END__.synctime,   MuseStruct_micro{ipart}.markers.Crise_End.synctime];
%                 MuseStruct_micro{ipart}.markers.BAD__END__.clock         = [MuseStruct_micro{ipart}.markers.BAD__END__.clock,      MuseStruct_micro{ipart}.markers.Crise_End.clock];
%             end
%         end       
%     end
%     
    % write data concatinated for SC, and update config with sampleinfo
    config{ipatient} = writeSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true, true);
    
    % create parameter and probe file for spyking circus
    writeSpykingCircusParameters(config{ipatient});
        
    % read raw spike data from SC, and segment into trials, requires 
%     [SpikeRaw, SpikeTrials]                 = readSpykingCircus_phy(config{ipatient}, MuseStruct_micro, true);
%     [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, true);

    % read and plot spikerate overview, and get the stats
%     [SpikeRateStats{ipatient}, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, true);
 

    
end
