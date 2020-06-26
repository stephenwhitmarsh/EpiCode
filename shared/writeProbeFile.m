function writeProbeFile(nb_channels, fname)

% WRITEPROBEFILE writes a SpykingCircus probe file
% use as
%   writeProbeFile(nb_channels, fname)

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

fid = fopen(fname,'w+');

fprintf(fid,'# Probe file for Spyking Circus \n');
fprintf(fid,'# Created automatically by writeSpykingCircus.m \n \n');

fprintf(fid,'total_nb_channels = %d;\n',nb_channels);
fprintf(fid,'radius            = 10;\n\n');
fprintf(fid,'channel_groups = {\n');
fprintf(fid,'\t1: {\n');
fprintf(fid,"\t\t'channels':range(0,%d),\n",nb_channels);
fprintf(fid,"\t\t'geometry': {\n");
for i = 1 : nb_channels
    fprintf(fid,'\t\t\t\t\t\t%d: [0, %d],\n',i-1,i*50);
end
fprintf(fid,'\t\t},\n');
fprintf(fid,"\t\t'graph' : []\n");
fprintf(fid,"\t}\n");
fprintf(fid,"}\n");
fclose(fid);
