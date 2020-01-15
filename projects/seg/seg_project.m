%% Analysis script 
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/shared/
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/scripts/epilepsy/seg
addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip/
% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/mlib6/
% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/subaxis/
ft_defaults

addpath Z:/scripts/epilepsy/seg/
addpath Z:/scripts/epilepsy/shared/
addpath Z:/fieldtrip/
ft_defaults

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% to not overload server
% maxNumCompThreads(4)


%% General analyses


for ipatient = 1
   
    % load settings
    config = seg_setparams([]);
    
    % read muse markers
    [MuseStruct_micro, MuseStruct_macro] = readMuseMarkers_parts(config{ipatient}, false); % false = reload from disk, true = redo the analysis if anything has changed
    
    [MuseStruct_micro, MuseStruct_macro] = MuseMarkers_update_filepath_parts(config{ipatient},MuseStruct_micro, MuseStruct_macro);

    switch ipatient  
        case 1   
            % move all seizures in first file to artefact.
            MuseStruct_micro{1}{1}.markers.BAD__START__.offset  (end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseStart.offset))    = MuseStruct_micro{1}{1}.markers.CriseStart.offset;
            MuseStruct_micro{1}{1}.markers.BAD__START__.synctime(end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseStart.synctime))  = MuseStruct_micro{1}{1}.markers.CriseStart.synctime;
            MuseStruct_micro{1}{1}.markers.BAD__START__.clock   (end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseStart.clock))     = MuseStruct_micro{1}{1}.markers.CriseStart.clock;
            MuseStruct_micro{1}{1}.markers.BAD__END__.offset    (end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseEnd.offset))      = MuseStruct_micro{1}{1}.markers.CriseEnd.offset;
            MuseStruct_micro{1}{1}.markers.BAD__END__.synctime  (end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseEnd.synctime))    = MuseStruct_micro{1}{1}.markers.CriseEnd.synctime;
            MuseStruct_micro{1}{1}.markers.BAD__END__.clock     (end+1:end+length(MuseStruct_micro{1}{1}.markers.CriseEnd.clock))       = MuseStruct_micro{1}{1}.markers.CriseEnd.clock;
            
            % move all seizures in second file to artefact.        
            MuseStruct_micro{1}{2}.markers.BAD__START__.offset  (end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseStart.offset))    = MuseStruct_micro{1}{2}.markers.CriseStart.offset;
            MuseStruct_micro{1}{2}.markers.BAD__START__.synctime(end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseStart.synctime))  = MuseStruct_micro{1}{2}.markers.CriseStart.synctime;
            MuseStruct_micro{1}{2}.markers.BAD__START__.clock   (end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseStart.clock))     = MuseStruct_micro{1}{2}.markers.CriseStart.clock;
            MuseStruct_micro{1}{2}.markers.BAD__END__.offset    (end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseEnd.offset))      = MuseStruct_micro{1}{2}.markers.CriseEnd.offset;
            MuseStruct_micro{1}{2}.markers.BAD__END__.synctime  (end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseEnd.synctime))    = MuseStruct_micro{1}{2}.markers.CriseEnd.synctime;
            MuseStruct_micro{1}{2}.markers.BAD__END__.clock     (end+1:end+length(MuseStruct_micro{1}{2}.markers.CriseEnd.clock))       = MuseStruct_micro{1}{2}.markers.CriseEnd.clock;
           
            % extend last seizure to end of file
            MuseStruct_micro{1}{2}.markers.BAD__START__.offset(end+1)   = MuseStruct_micro{1}{2}.markers.CriseStart.offset(4);
            MuseStruct_micro{1}{2}.markers.BAD__START__.synctime(end+1) = MuseStruct_micro{1}{2}.markers.CriseStart.synctime(4);
            MuseStruct_micro{1}{2}.markers.BAD__START__.clock(end+1)    = MuseStruct_micro{1}{2}.markers.CriseStart.clock(4);
            MuseStruct_micro{1}{2}.markers.BAD__END__.offset(end+1)     = MuseStruct_micro{1}{2}.markers.CriseStart.endsample;
            MuseStruct_micro{1}{2}.markers.BAD__END__.synctime(end+1)   = MuseStruct_micro{1}{2}.markers.CriseStart.endsample/MuseStruct_micro{1}{2}.Fs;
            MuseStruct_micro{1}{2}.markers.BAD__END__.clock(end+1)      = MuseStruct_micro{1}{2}.markers.BAD__END__.clock(end) + seconds(MuseStruct_micro{1}{2}.markers.CriseStart.endsample);      
    end
    
    % write data concatinated for SC, and update config with sampleinfo
    config{ipatient} = removefields( config{ipatient},'fnames_ncs');
    config{ipatient} = writeSpykingCircus_parts(config{ipatient}, MuseStruct_micro, false, false);
    
    writeSpykingCircusParameters_parts(config{ipatient});

    [SpikeRaw, SpikeTrials] = readSpykingCircus_parts(config{ipatient}, MuseStruct_micro, true);
    [SpikeRateStats{ipatient}, stats_bar, sdf_orig_out, sdf_bar_out] = spikeratestats(config{ipatient}, SpikeRaw, SpikeTrials, true, 1);


end
    
    % plot seizures
%     [MuseStruct_micro, MuseStruct] = readMuseMarkers_parts(config{ipatient}, true); % false = reload from disk, true = redo the analysis if anything has changed    
%     plotSeizureSegmentation(config{ipatient},MuseStruct_micro);
    
    % pre-ictal period, second file always has the seizure
%     [MuseStruct_micro, MuseStruct_macro] = readMuseMarkers_parts(config{ipatient}, true); % false = reload from disk, true = redo the analysis if anything has changed
    
    
