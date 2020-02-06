function writeMuseMarkerfile(MuseStruct,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function writeMuseMarkers(MuseStruct)
%
% Writes markerfiles for Muse. The input structure can be read by
% readMuseMarkers.m, and edited before writing it back into the Muse file.
%
% Example:
%
% writeMuseMarkers(MuseStruct);
%
% Code by Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names           = fieldnames(MuseStruct.markers);
mrk_number      = length(fieldnames(MuseStruct.markers));
[fid, message]  = fopen(fname,'w+'); % open/create and discard content

if message ~= 0
    fprintf('ERROR, something went wrong with reading the markerfile %s\n',fname);
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
    fprintf(fid,'%s\n', MuseStruct.markers.(names{imarker}).classgroupid);
    fprintf(fid,'%s\n', 'NAME:');
    fprintf(fid,'%s\n', names{imarker});
    fprintf(fid,'%s\n', 'COMMENT:');
    fprintf(fid,'%s\n', MuseStruct.markers.(names{imarker}).comment);
    fprintf(fid,'%s\n', 'COLOR:');
    fprintf(fid,'%s\n', MuseStruct.markers.(names{imarker}).color);
    fprintf(fid,'%s\n', 'EDITABLE:');
    fprintf(fid,'%s\n', MuseStruct.markers.(names{imarker}).editable);
    fprintf(fid,'%s\n', 'CLASSID:');
    fprintf(fid,sprintf('+%d\n',imarker)); 
%     fprintf(fid,'%s\n', MuseStruct.markers.(names{imarker}).classid);
    
    if isfield(MuseStruct.markers.(names{imarker}),'synctime')
        nr_samples = length(MuseStruct.markers.(names{imarker}).synctime);
    else
        nr_samples = 0;
    end
    
    fprintf(fid,'%s\n', 'NUMBER OF SAMPLES:');
    fprintf(fid,'%s\n', num2str(nr_samples));
    fprintf(fid,'%s\n', 'LIST OF SAMPLES:');
    fprintf(fid,'%s\n', 'TRIAL NUMBER		TIME FROM SYNC POINT (in seconds)');
    
    if nr_samples > 0
        % add samples
        fprintf('Wrote %d events of: %s\n',length(MuseStruct.markers.(names{imarker}).synctime),names{imarker});
        for itrial = 1 : length(MuseStruct.markers.(names{imarker}).synctime)
            %             if ~isfield(MuseStruct.markers.(names{imarker}),'trialnum')
            fprintf(fid,'                  +0\t\t\t\t+%.10f\n',MuseStruct.markers.(names{imarker}).synctime(itrial));
            %             else
            %                 fprintf(fid,'                  +%d\t\t\t\t+%.10f\n',MuseStruct.markers.(names{imarker}).trialnum(itrial),MuseStruct.markers.(names{imarker}).synctime(itrial))
            %             end
        end
    end
    fprintf(fid,'\n\n');
    
end

exitcode = fclose(fid);
if exitcode ~= 0
    disp('ERROR, something went wrong with writing the markerfile. Is it open?\n');
else
    fprintf('Succesfully written markers to: %s\n',fname);
end

