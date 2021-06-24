function writeSpykingCircusDeadFile_concatinated(MuseStruct)

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

deadfile_ms         = [];
deadfile_samples    = [];
last_ms             = 0;
last_samples        = 0;
dirlist             = [];

for idir = 1 : size(MuseStruct,2)
    
    if isfield(MuseStruct{idir},'markers')
        if isfield(MuseStruct{idir}.markers,'BAD__START__')
            if isfield(MuseStruct{idir}.markers.BAD__START__,'synctime')
                
                %                 filename = fullfile(MuseStruct{idir}.directory,[MuseStruct{idir}.filenames(1).name(1:end-4),'.dead']);
                %                 fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
                
                % check if there is an equal amount of start and end markers
                if size(MuseStruct{idir}.markers.BAD__START__.events,2)-size(MuseStruct{idir}.markers.BAD__END__.events,2) == 0
                    fprintf('Great, recovered same number of starts and end markers \n')
                    %                     dlmwrite(filename,[MuseStruct{idir}.markers.BAD__START__.synctime'*1000,MuseStruct{idir}.markers.BAD__END__.synctime'*1000],'delimiter','	','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense
                    
                    deadfile_ms         = [deadfile_ms;         MuseStruct{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                    deadfile_samples    = [deadfile_samples;    MuseStruct{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{idir}.markers.BAD__END__.offset'+last_samples];
                    
                    hdr                 = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(1).name));
                    last_samples        = last_samples + hdr.nSamples;
                    last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                    
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
                            MuseStruct{idir}.markers.BAD__END__.offset(x) = [];
                            
                        end
                    end
                    %                     dlmwrite(filename,[MuseStruct{idir}.markers.BAD__START__.synctime'*1000,MuseStruct{idir}.markers.BAD__END__.synctime'*1000],'delimiter','	','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense
                    deadfile_ms         = [deadfile_ms;         MuseStruct{idir}.markers.BAD__START__.synctime'*1000+last_ms,      MuseStruct{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                    deadfile_samples    = [deadfile_samples;    MuseStruct{idir}.markers.BAD__START__.offset'+last_samples,        MuseStruct{idir}.markers.BAD__END__.offset'+last_samples];
                    
                    hdr                 = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(1).name));
                    last_samples        = last_samples + hdr.nSamples;
                    last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                    
                end
                
            else
                fprintf('Found no artefacts for file: %s\n',MuseStruct{idir}.directory);
            end
        else
            fprintf('Found no artefacts for file: %s\n',MuseStruct{idir}.directory);
        end
    end
    dirlist = [dirlist; MuseStruct{idir}.directory];
    fprintf('%d\n',idir);
end

filename = fullfile(fileparts(MuseStruct{idir}.directory),'artefacts_ms.dead');
fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
dlmwrite(filename,deadfile_ms,'delimiter','	','precision','%.4f');


filename = fullfile(fileparts(MuseStruct{idir}.directory),'artefacts_samples.dead');
fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
dlmwrite(filename,deadfile_samples,'delimiter','	','precision','%.4f');

filename = fullfile(fileparts(MuseStruct{idir}.directory),'dirlist.txt');
fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
dlmwrite(filename,dirlist);

fid = fopen(filename,'w');
for r=1:size(dirlist,1)
    fprintf(fid,'%s\n',dirlist(r,:));
end
fclose(fid);
