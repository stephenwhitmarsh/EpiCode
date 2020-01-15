function writeMuseMarkers(MuseStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function writeMuseMarkers(MuseStruct)
%
% Crawls and searches through patient- and recording-directories and
% writes markerfiles for Muse. The input structure can be read by
% readMuseMarkers.m, and edited before writing it back into the Muse file.
%
% Example:
% 
% writeMuseMarkers(MuseStruct);
%
% Code by Stephen Whitmarsh (stephen.whitmarsh@icm-institute)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idir = 1 : size(MuseStruct,1)
    
    names           = fieldnames(MuseStruct{idir}.markers);
    mrk_number      = length(fieldnames(MuseStruct{idir}.markers));
    fname           = fullfile(MuseStruct{idir}.directory,'Events_test.mrk');
    [fid, message]  = fopen(fname,'w');
    
    if message ~= 0
        disp('ERROR, something went wrong with reading the markerfile %s',fname);
        disp(message);
        return
    else
        fprintf('Succesfully opened %s\n',fname);
    end
    
    % add header
    fprintf(fid,'%s\n', 'PATH OF DATASET:');
    fprintf(fid,'\n\n\n');
    fprintf(fid,'%s\n', 'NUMBER OF MARKERS:');
    fprintf(fid,'%s\n',num2str(mrk_number));
    fprintf(fid,'\n\n');
    
    % add marker
    for imarker = 1 : mrk_number
        
        fprintf(fid,'%s\n', 'CLASSGROUPID:');
        fprintf(fid,'%s\n', MuseStruct{idir}.markers.(names{imarker}).classgroupid);
        fprintf(fid,'%s\n', 'NAME:');
        fprintf(fid,'%s\n', names{imarker});
        fprintf(fid,'%s\n', 'COMMENT:');
        fprintf(fid,'%s\n', MuseStruct{idir}.markers.(names{imarker}).comment);
        fprintf(fid,'%s\n', 'COLOR:');
        fprintf(fid,'%s\n', MuseStruct{idir}.markers.(names{imarker}).color);
        fprintf(fid,'%s\n', 'EDITABLE:');
        fprintf(fid,'%s\n', MuseStruct{idir}.markers.(names{imarker}).editable);
        fprintf(fid,'%s\n', 'CLASSID:');
        fprintf(fid,'%s\n', MuseStruct{idir}.markers.(names{imarker}).classid);
        
        if isfield(MuseStruct{idir}.markers.(names{imarker}),'synctime')
            nr_samples = size(MuseStruct{idir}.markers.(names{imarker}).synctime,1);
        else
            nr_samples = 0;
        end
        
        fprintf(fid,'%s\n', 'NUMBER OF SAMPLES:');
        fprintf(fid,'%s\n', num2str(nr_samples));
        fprintf(fid,'%s\n', 'LIST OF SAMPLES:');
        fprintf(fid,'%s\n', 'TRIAL NUMBER		TIME FROM SYNC POINT (in seconds)');
        
        if nr_samples > 0
            % add samples
            for itrial = 1 : size(MuseStruct{idir}.markers.(names{imarker}).synctime,1)
                fprintf(fid,'                  +%d\t\t\t\t+%.10f\n',MuseStruct{idir}.markers.(names{imarker}).trialnum(itrial),MuseStruct{idir}.markers.(names{imarker}).synctime(itrial));
            end
        end
        fprintf(fid,'\n\n'); 
   
    end
    
    exitcode = fclose(fid);
    if exitcode ~= 0
        disp('ERROR, something went wrong with writing the markerfile. Is it open?');
    else
        fprintf('Succesfully written markers to: %s\n',fname);
    end
    
end
