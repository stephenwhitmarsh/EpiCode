function writeMuseMarkerfile(MuseStruct, fname)

% WRITEMUSEMARKERFILE writes markerfiles for Muse.
% use as
%   function writeMuseMarkers(MuseStructt, fname)
%
% The input structure can be read by readMuseMarkers.m, and edited before
% writing it back into the Muse file.

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

names           = fieldnames(MuseStruct.markers);
mrk_number      = length(fieldnames(MuseStruct.markers));
[fid, message]  = fopen(fname, 'w+'); % open/create and discard content

if message ~= 0
    fprintf('ERROR, something went wrong with reading the markerfile %s\n',fname);
    return
else
    fprintf('Successfully opened %s\n',fname);
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

    if isfield(MuseStruct.markers.(names{imarker}),'synctime')
        nr_samples = length(MuseStruct.markers.(names{imarker}).synctime);
    else
        nr_samples = 0;
    end

    fprintf(fid,'%s\n', 'NUMBER OF SAMPLES:');
    fprintf(fid,'%s\n', num2str(nr_samples));
    fprintf(fid,'%s\n', 'LIST OF SAMPLES:');
    fprintf(fid,'%s\n', 'TRIAL NUMBER		TIME FROM SYNC POINT (in seconds)');
    
    if nr_samples == 0
        continue
    end
    
    fprintf('Writing %d events of: %s\n',length(MuseStruct.markers.(names{imarker}).synctime),names{imarker});
    for itrial = 1 : length(MuseStruct.markers.(names{imarker}).synctime)
        fprintf(fid,'                  +0\t\t\t\t+%.10f\n',MuseStruct.markers.(names{imarker}).synctime(itrial));
    end
    fprintf(fid,'\n\n');
end

exitcode = fclose(fid);
if exitcode ~= 0
    disp('ERROR, something went wrong with writing the markerfile. Is it open?\n');
else
    fprintf('Succesfully written markers to: %s\n', fname);
end
