function [MuseStruct]  = readMuseMarker(name_mrk)

% READMUSEMARKER extracts marker timings (time and samples) created by Muse.
% use as
%    [MuseStruct]  = readMuseMarker(name_mrk)
%
% with help (code) from Jean-Didier Lemarechal & Craig Richter

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

if ~exist(name_mrk)
    error('%s not found', name_mrk);
end

f = fopen(name_mrk, 'rt');
markfile = {};
while true
    l = fgetl(f);
    if ~ischar(l)
        break
    end
    markfile{end + 1} = l;
end
fclose(f);

if ~isempty(markfile)
    nmarkers        = str2num(markfile{strmatch('NUMBER OF MARKERS:', markfile, 'exact') + 1});
    classgroupid    = markfile(strmatch('CLASSGROUPID:', markfile, 'exact') + 1);
    name            = markfile(strmatch('NAME:', markfile, 'exact') + 1);
    comment         = markfile(strmatch('COMMENT:', markfile, 'exact') + 1);
    color           = markfile(strmatch('COLOR:', markfile, 'exact') + 1);
    editable        = markfile(strmatch('EDITABLE:', markfile, 'exact') + 1);
    classid         = markfile(strmatch('CLASSID:', markfile, 'exact') + 1);
    nrEvents        = str2num(char(markfile(strmatch('NUMBER OF SAMPLES:', markfile, 'exact') + 1)));

    % Get the events, time in seconds from onset of file
    j = strmatch('LIST OF SAMPLES:', markfile, 'exact') + 2;
    for i = 1 : nmarkers
        marks{i} = str2num(char(markfile(j(i):j(i) + nrEvents(i) - 1)));
        fprintf('Found %d occurances of %s \n', size(marks{i},1), name{i});

        % Convert from index origin 0 to 1
        if nrEvents(i) ~= 0
            marks{i}(:, 1) = marks{i}(:, 1) + 1;
        end
    end

    for imarker = 1 : nmarkers
        name{imarker} = strrep(name{imarker},'-','_'); % cant make fieldnames with minusses
        MuseStruct.markers.(name{imarker}).events         = [];
        MuseStruct.markers.(name{imarker}).comment        = comment{imarker};
        MuseStruct.markers.(name{imarker}).color          = color{imarker};
        MuseStruct.markers.(name{imarker}).editable       = editable{imarker};
        MuseStruct.markers.(name{imarker}).classid        = classid{imarker};
        MuseStruct.markers.(name{imarker}).classgroupid   = classgroupid{imarker};

        for ievent = 1 : size(marks{imarker},1)
            MuseStruct.markers.(name{imarker}).trialnum                 = marks{imarker}(ievent,1);
            MuseStruct.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent,2);
        end

    end
else
    fprintf('\n\n %s is empty!!! \n\n',name_mrk);
end
