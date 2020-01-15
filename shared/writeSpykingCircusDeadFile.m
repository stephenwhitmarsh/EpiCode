function writeSpykingCircusDeadFile(MuseStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function writeSpykingCircusDeadFile(MuseStruct)
%
% Crawls and searches through patient- and recording-directories and writes
% artifact definitions for SpykingCircus. The input structure can be read
% by readMuseMarkers.m, and edited before writing it to artefact files for
% SpikingCircus
%
% Example:
% 
% writeSpykingCircusDeadFile(MuseStruct);
%
% Code by Stephen Whitmarsh (stephen.whitmarsh@icm-institute)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idir = 1 : size(MuseStruct,2)
    
    if isfield(MuseStruct{idir},'markers')
        if isfield(MuseStruct{idir}.markers,'BAD__START__')
            if isfield(MuseStruct{idir}.markers.BAD__START__,'synctime')
                
                filename = fullfile(MuseStruct{idir}.directory,[MuseStruct{idir}.filenames(1).name(1:end-4),'.dead']);
                fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
  
                % check if there is an equal amount of start and end markers
                if size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) == 0
                    fprintf('Great, recovered same number of starts and end markers \n')
                    dlmwrite(filename,[MuseStruct{idir}.markers.BAD__START__.synctime'*1000,MuseStruct{idir}.markers.BAD__END__.synctime'*1000],'delimiter','	','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense

                elseif size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) > 0
                    fprintf('ERROR! more starts than ends found in %s \n',MuseStruct{idir}.directory);
                elseif size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) < 0
                    fprintf('ERROR! more ends than starts found in %s - CORRECTING \n',MuseStruct{idir}.directory)
                    for i = 1 : 10
                        for itrial = 1 : length(MuseStruct{idir}.markers.BAD__START__.synctime)
                            start(itrial) = MuseStruct{idir}.markers.BAD__START__.synctime(itrial);
                            stop(itrial)  = MuseStruct{idir}.markers.BAD__END__.synctime(itrial);
                        end
                        
                        x = find(start > stop,1,'first');
                        if ~isempty(x)
                            MuseStruct{idir}.markers.BAD__END__.synctime(x) = [];
                        end
                    end
                    dlmwrite(filename,[MuseStruct{idir}.markers.BAD__START__.synctime'*1000,MuseStruct{idir}.markers.BAD__END__.synctime'*1000],'delimiter','	','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense

                end
                
            end
        end
    end
end
