function dtx_project_spikes(slurm_task_id)

if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    
end

ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx




%% prepare data for spyking circus analysis
%     for irat = 3:5
%         [MuseStruct]                     = readMuseMarkers(config{irat}, false);
%         [MuseStruct]                     = dtx_remove_wrong_seizure(config{irat}, MuseStruct,true, true);
%         %remove all seizures
%         [MuseStruct_without_seizures]    = addMuseBAD(MuseStruct, 'all', 'all', 'SlowWave', 'Crise_End', 'all', -2, 2);
%         %multifile and dead file
%         writeSpykingCircus(config{irat}, MuseStruct_without_seizures, true, true);
%         %write deadfile with the seizures
%         writeSpykingCircusDeadFile(config{irat},MuseStruct, 'withSeizures', 'all');
%         %param file and probe file
%         writeSpykingCircusParameters(config{irat})
%         %FIXME : vérifier que l'ajout des BAD a bien fonctionné
%     end

config = dtx_setparams_probe_spikes([]);


%% analyse spyking circus output

for irat = slurm_task_id
    

%     config{irat}.imagesavedir = fullfile(config{irat}.imagesavedir,'Temporary_check_remove_seizures'); %REMOVEME FIXME
    
    %align markers and remove wrong seizures
    [MuseStruct]                    = readMuseMarkers(config{irat}, false);
    [MuseStruct]                    = alignMuseMarkers(config{irat},MuseStruct, false);
    [MuseStruct]                    = dtx_remove_wrong_seizure(config{irat}, MuseStruct,false, false);
    
    %read LFP data
    dat_LFP                         = readLFP(config{irat}, MuseStruct, false, false);
    [config{irat},dat_LFP]          = dtx_correctDTX2name(config{irat},dat_LFP);
    
    %read spike data
    [SpikeRaw]                      = readSpikeRaw_Phy(config{irat},false,'all');
    [SpikeTrials]                   = readSpykingCircus(config{irat}, MuseStruct,SpikeRaw, false, 'all'); %FIXME rename readSpikeTrials_Muse
    
    %remove BAD LFP/Spike trials
    cfgtemp                         = [];
    cfgtemp                         = config{irat};
    cfgtemp.LFP.electrodetoplot     = {config{irat}.align.channel{1},config{irat}.align.channel{1},config{irat}.align.channel{1}};
    cfgtemp.method                  = 'remove';
    cfgtemp.markerstart             = 'BAD__START__';
    cfgtemp.markerend               = 'BAD__END__';
    cfgtemp.indexstart              = 0;
    cfgtemp.indexend                = 0;
    cfgtemp.timefrombegin           = 0;
    cfgtemp.timefromend             = 0;
    cfgtemp.plotdata                = 'yes';
    [dat_LFP, ~]                    = removetrials_MuseMarkers(cfgtemp, dat_LFP, MuseStruct, 'all', 'all');
    [SpikeTrials, ~]                = removetrials_MuseMarkers(cfgtemp, SpikeTrials, MuseStruct, 'all', 'all');
    
    %read spike waveforms
    [SpikeWaveforms]                = readSpikeWaveforms(config{irat}, SpikeTrials, false, 'all');
    %FIXME remove _1000 for loading precomputed data
    
    %stats per unit
    stats                           = spikeratestats_Events_Interictal(config{irat},SpikeTrials,SpikeWaveforms,dat_LFP,'all',true);
    %FIXME : remove try catch x2
    %FIXME : ajouter doublets de PA avec spikestatsOverTime
    
    
    clear all
end %irat

end

