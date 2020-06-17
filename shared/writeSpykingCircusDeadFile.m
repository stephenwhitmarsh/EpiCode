function writeSpykingCircusDeadFile(cfg, MuseStruct, suffix, part_list)

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
% Paul Baudin : add suffix, to more easily create several different dead
% files for the same analysis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(part_list, 'all')
    part_list = 1:size(MuseStruct,2);
end

for ipart = part_list
    
    deadfile_ms         = [];
    deadfile_samples    = [];
    last_ms             = 0;
    last_samples        = 0;
    dirlist{ipart}      = [];
    
    for idir = 1 : size(MuseStruct{ipart},2)
        if isfield(MuseStruct{ipart}{idir},'markers')
            if isfield(MuseStruct{ipart}{idir}.markers,'BAD__START__')
                if isfield(MuseStruct{ipart}{idir}.markers.BAD__START__,'synctime')
                    
                    % check if there is an equal amount of start and end markers
                    if size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) == 0
                        
                        temp = [];
                        hdr = [];
                        temp    = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.channel{1},'.ncs'])); %firts chan. all chans have the same nr of samples
                        hdr     = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, temp.name));
                        
                        fprintf('Great, recovered same number of start and end markers \n')
                        deadfile_ms         = [deadfile_ms;         MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                        deadfile_samples    = [deadfile_samples;    MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples, MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
                        last_samples        = last_samples  + hdr.nSamples;
                        last_ms             = last_ms       + hdr.nSamples/hdr.Fs * 1000;
                        
                    elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) > 0
                        
                        fprintf('ERROR! more start than end found in %s \n',MuseStruct{ipart}{idir}.directory);
                        
                    elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) < 0
                        
                        fprintf('ERROR! more end than start found in %s - CORRECTING \n',MuseStruct{ipart}{idir}.directory)
                        for i = 1 : 10
                            for itrial = 1 : length(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime)
                                start(itrial) = MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(itrial);
                                stop(itrial)  = MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(itrial);
                            end
                            
                            x = find(start > stop,1,'first');
                            if ~isempty(x)
                                MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(x) = [];
                            end
                        end
                        
                        deadfile_ms         = [deadfile_ms;         MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                        deadfile_samples    = [deadfile_samples;    MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples, MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
                        last_samples        = last_samples + hdr.nSamples;
                        last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                    end
                else
                    fprintf('Found no artefacts for file: %s\n',MuseStruct{ipart}{idir}.directory);
                end
            else
                fprintf('Found no artefacts for file: %s\n',MuseStruct{ipart}{idir}.directory);
            end
        end
        dirlist{ipart} = [dirlist{ipart}; MuseStruct{ipart}{idir}.directory];
        fprintf('%d\n',idir);
    end % idir
    
    % write artefacts to txt file for spyking circus
    subjdir     = cfg.prefix(1:end-1);
    partdir     = ['p',num2str(ipart)];
    
    filename    = sprintf('SpykingCircus_artefacts_ms%s.dead', suffix);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(fullfile(cfg.datasavedir,subjdir,partdir,filename),deadfile_ms,'delimiter','	','precision','%.4f');
    
    filename = sprintf('SpykingCircus_artefacts_samples%s.dead', suffix);
    fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
    dlmwrite(fullfile(cfg.datasavedir,subjdir,partdir,filename),deadfile_samples,'delimiter','	','precision','%.4f');
    
end %ipart