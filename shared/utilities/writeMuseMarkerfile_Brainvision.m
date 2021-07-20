function writeMuseMarkerfile_Brainvision(MuseStruct, fname)

% WRITEMUSEMARKERFILE writes markerfiles for Muse from Brainvision data format.
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

[folder, filename] = fileparts(fname);
hdr = ft_read_header(fullfile(folder, [filename, '.eeg']));

if message ~= 0
    fprintf('ERROR, something went wrong with reading the markerfile %s\n',fname);
    return
else
    ft_info('Successfully opened %s\n',fname);
end

% add header
fprintf(fid,'#PTX_V2.0\n');
fprintf(fid,'#ASSOCIATED HEADER FILE\n');
fprintf(fid,'#FORMAT: TRIAL	TIME	NAME	DURATION	OFFSET\n');
fprintf(fid,'#TRIGGER COMMENTS\n');

% add marker names
for imarker = 1 : mrk_number
    fprintf(fid,'#	%s	mrk	#color=#%s\n', names{imarker}, MuseStruct.markers.(names{imarker}).color);
end

% add the end of the header
fprintf(fid,'Brain Vision Data Exchange Marker File, Version 1.0\n');
fprintf(fid,'[Common Infos]\n');
fprintf(fid,'Codepage=UTF-8\n');
fprintf(fid,'DataFile=\n');
fprintf(fid,'[Marker Infos]\n');

% sort and write events
ievent = 0;
for imarker = 1 : mrk_number
    if ~isfield(MuseStruct.markers.(names{imarker}), 'synctime')
        continue
    end
    if isempty(MuseStruct.markers.(names{imarker}).synctime)
        continue
    end
    for i = 1:size(MuseStruct.markers.(names{imarker}).synctime, 2)
        ievent              = ievent +1;
        mrk_name{ievent}    = names{imarker};
        event_time(ievent)  = MuseStruct.markers.(names{imarker}).synctime(i);
    end
end
if ievent > 0
    [~, sort_idx]   = sort(event_time);
    event_time      = event_time(sort_idx);
    mrk_name        = mrk_name(sort_idx);
    event_sample    = round(event_time .* hdr.Fs);
    
    %write events
    for ievent = 1:size(event_sample, 2)
        fprintf(fid,'Mk%d=mrk,%s,%d,1,0\n', ievent, mrk_name{ievent}, event_sample(ievent));
    end
end

exitcode = fclose(fid);
if exitcode ~= 0
    disp('ERROR, something went wrong with writing the markerfile. Is it open?\n');
else
    ft_info('Succesfully written markers to: %s\n', fname);
end
