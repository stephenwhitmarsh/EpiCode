%% Analysis script
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/shared/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/seg
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip/
ft_defaults

% addpath Z:/scripts/epilepsy/seg/
% addpath Z:/scripts/epilepsy/shared/
% addpath Z:/fieldtrip/
% ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% General analyses

for ipatient = 1
    
    % load settings
    config = seg_setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro] = readMuseMarkers_parts(config{ipatient}, false); % false = reload from disk, true = redo the analysis if anything has changed
    
    [MuseStruct_micro, MuseStruct_macro] = MuseMarkers_update_filepath_parts(config{ipatient},MuseStruct_micro, MuseStruct_macro);
    
    %%
    for ipart = 1 : size(config{ipatient}.directorylist,2)
        for ifile = 1 : 2
            % move all seizures in first file to artefact.
            % only when seizures are present
            if isfield(MuseStruct_micro{ipart}{ifile}.markers,'Crise_Start')
                if isfield(MuseStruct_micro{ipart}{ifile}.markers.Crise_Start,'offset')
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__START__.offset  (end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseStart.offset))    = MuseStruct_micro{ipart}{ifile}.markers.CriseStart.offset;
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__START__.synctime(end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseStart.synctime))  = MuseStruct_micro{ipart}{ifile}.markers.CriseStart.synctime;
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__START__.clock   (end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseStart.clock))     = MuseStruct_micro{ipart}{ifile}.markers.CriseStart.clock;
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__END__.offset    (end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.offset))      = MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.offset;
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__END__.synctime  (end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.synctime))    = MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.synctime;
                    MuseStruct_micro{ipart}{ifile}.markers.BAD__END__.clock     (end+1:end+length(MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.clock))       = MuseStruct_micro{ipart}{ifile}.markers.CriseEnd.clock;
                end
            end
        end
    end
    
            %%
    
    % write data concatinated for SC, and update config with sampleinfo
    %     config{ipatient} = removefields( config{ipatient},'fnames_ncs');
    config{ipatient} = writeSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true, true);
    
%     writeSpykingCircusParameters_parts(config{ipatient});
    
    %     [SpikeRaw, SpikeTrials] = readSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true);
    %     [SpikeRateStats{ipatient}, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, true, 1);
    
   
end

