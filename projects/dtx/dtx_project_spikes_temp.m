function dtx_project_spikes_temp(slurm_task_id)

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
end

ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx




%% cfg modifs for batching
sc_preparation = false;

config = dtx_setparams_probe_spikes([]);
irat = 1;

config{irat}.circus.postfix = '-thresh12-merge_after_tempmatch';
config{irat}.imagesavedir = fullfile(config{1}.imagesavedir, sprintf('spikerate%s',config{irat}.circus.postfix));

%fprintf('for postfix %s \n', circus_names{slurm_task_id});

%% 
pause(60*(slurm_task_id-1)); %to avoid loading the same data at the same time

[MuseStruct]                     = readMuseMarkers(config{irat}, false);

% align Muse markers according to peaks and detect whether they contain artefacts
[MuseStruct]                     = alignMuseMarkers(config{irat},MuseStruct, false);

[MuseStruct]                     = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);


%% Create files for spyking circus

%     %remove seizure from 5s after SlowWave to Crise_End
%     [MuseStruct_without_seizures]    = addMuseBAD(MuseStruct, 'all', 'all', 'SlowWave', 'Crise_End', 'all', -2, 2);
%     [sampleinfo] = writeSpykingCircus(config{irat}, MuseStruct_without_seizures, true, true);
%     writeSpykingCircusParameters(config{irat}) %Corriger les noms des fichiers lus par SC

    
%% read spike-clustering results, and epoch around events

[SpikeRaw] = readSpykingCircus_SpikeRaw(config{irat},false,'all');
[SpikeTrials] = readSpykingCircus(config{irat}, MuseStruct,SpikeRaw, false, 1);
statsspiketrials_norm(config{irat},SpikeTrials,1,3, slurm_task_id);

if slurm_task_id == 1
    config{1}.spikestage = [];
    
    config{1}.spikestage{1}.markerstart          = 'Baseline_Start';
    config{1}.spikestage{1}.markerstartoffset    = 0; %seconds
    config{1}.spikestage{1}.markerend            = 'Injection';
    config{1}.spikestage{1}.indexEnd             = 0;
    config{1}.spikestage{1}.periodmin            = 600; %seconds
    config{1}.spikestage{1}.triallength          = 60; %seconds
    config{1}.spikestage{1}.subdivision          = 'all'; %minimum 2
    config{1}.spikestage{1}.includeend           = false;
    config{1}.spikestage{1}.markerendoffset      = 0; %seconds
    config{1}.spikestage{1}.removeartefact       = 'trial';
    config{1}.spikestage{1}.stagenumber          = 0; %array of integers. to give the same number to all the trials, write one single value
    
    config{1}.spikestage{2}.markerstart          = 'Crise_End';
    config{1}.spikestage{2}.markerstartoffset    = 2; %seconds
    config{1}.spikestage{2}.markerend            = 'SlowWave';
    config{1}.spikestage{2}.indexEnd             = 1; %if 1 : take the i start marker and the i+1 end marker
    config{1}.spikestage{2}.periodmin            = 20; %seconds
    config{1}.spikestage{2}.triallength          = 10; %seconds
    config{1}.spikestage{2}.subdivision          = 'all'; %minimum 2
    config{1}.spikestage{2}.includeend           = false; %add the end of the period to the previous subdivisions
    config{1}.spikestage{2}.markerendoffset      = -2; %seconds
    config{1}.spikestage{2}.removeartefact       = 'period'; %period, trial, none
    config{1}.spikestage{2}.stagenumber          = 1; %array of integers. if single value : gives the same number to all the trials


config{1}.imagesavedir = fullfile(config{1}.imagesavedir,'Compare_Baseline');
% spikeratestatsEvents(config{irat}, SpikeRaw, SpikeTrials, true);

SpikeTrials = [];
[SpikeTrials] = readSpykingCircus_spikestage(config{irat}, SpikeRaw, MuseStruct, 1, 'all', true);
[stats] = spikeratestatsPSG(config{irat},SpikeRaw,SpikeTrials,true,true);

end
%
% %% General analyses
% 
% for ipatient = 1:1
%     
%     config = dtx_setparams([]);    
%     
%     % read muse markers
%     [MuseStruct_micro, MuseStruct_macro]    = readMuseMarkers(config{ipatient}, false);
%     
%     % align Muse markers according to peaks and detect whether they contain artefacts
%     [MuseStruct_micro, MuseStruct_macro]    = alignMuseMarkers(config{ipatient},MuseStruct_micro, MuseStruct_macro, false);   
%     [MuseStruct_micro, MuseStruct_macro]    = MuseMarkers_update_filepath(config{ipatient},MuseStruct_micro, MuseStruct_macro);
%     
%     
%     % read LFP data
% %    [dat_micro, dat_macro] = readLFP(config{ipatient}, MuseStruct_micro, MuseStruct_macro, false, false);
%     
%     % create 'artefacts' to remove seizures and time from seizures
% %     for ipart = 1 : size(MuseStruct_micro,2)
% %         
% %         % only when seizures are present
% %         if isfield(MuseStruct_micro{ipart}.markers,'Crise_Start')
% %             if isfield(MuseStruct_micro{ipart}.markers.Crise_Start,'offset')
% %                 
% %                 % if no artefact fields exist, create empty ones
% %                 if ~isfield(MuseStruct_micro{ipart}.markers,'BAD__START__')
% %                     MuseStruct_micro{ipart}.markers.BAD__START__.offset      = [];
% %                     MuseStruct_micro{ipart}.markers.BAD__START__.synctime    = [];
% %                     MuseStruct_micro{ipart}.markers.BAD__START__.clock       = [];
% %                     
% %                     MuseStruct_micro{ipart}.markers.BAD__END__.offset        = [];
% %                     MuseStruct_micro{ipart}.markers.BAD__END__.synctime      = [];
% %                     MuseStruct_micro{ipart}.markers.BAD__END__.clock         = [];
% %                 end
% %                 
% %                 MuseStruct_micro{ipart}.markers.BAD__START__.offset      = [MuseStruct_micro{ipart}.markers.BAD__START__.offset,   MuseStruct_micro{ipart}.markers.Crise_Start.offset];
% %                 MuseStruct_micro{ipart}.markers.BAD__START__.synctime    = [MuseStruct_micro{ipart}.markers.BAD__START__.synctime, MuseStruct_micro{ipart}.markers.Crise_Start.synctime];
% %                 MuseStruct_micro{ipart}.markers.BAD__START__.clock       = [MuseStruct_micro{ipart}.markers.BAD__START__.clock,    MuseStruct_micro{ipart}.markers.Crise_Start.clock];
% %                 
% %                 MuseStruct_micro{ipart}.markers.BAD__END__.offset        = [MuseStruct_micro{ipart}.markers.BAD__END__.offset,     MuseStruct_micro{ipart}.markers.Crise_End.offset];
% %                 MuseStruct_micro{ipart}.markers.BAD__END__.synctime      = [MuseStruct_micro{ipart}.markers.BAD__END__.synctime,   MuseStruct_micro{ipart}.markers.Crise_End.synctime];
% %                 MuseStruct_micro{ipart}.markers.BAD__END__.clock         = [MuseStruct_micro{ipart}.markers.BAD__END__.clock,      MuseStruct_micro{ipart}.markers.Crise_End.clock];
% %             end
% %         end       
% %     end
% %     
%     % write data concatinated for SC, and update config with sampleinfo
%     config{ipatient} = writeSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true, true);
%     
%     % create parameter and probe file for spyking circus
%     writeSpykingCircusParameters(config{ipatient});
%         
%     % read raw spike data from SC, and segment into trials, requires 
% %     [SpikeRaw, SpikeTrials]                 = readSpykingCircus_phy(config{ipatient}, MuseStruct_micro, true);
% %     [SpikeRaw, SpikeTrials]                 = readSpykingCircus(config{ipatient}, MuseStruct_micro, true);
% 
%     % read and plot spikerate overview, and get the stats
% %     [SpikeRateStats{ipatient}, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, true);
%  
% 
%     
% end
