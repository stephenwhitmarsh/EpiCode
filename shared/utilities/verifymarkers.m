function wrong_markers = verifymarkers(cfg, MuseStruct)

% Use as :
%       wrong_markers = verifymarkers(cfg, MuseStruct)
% 
% verifymarkers check the marker list for the 2 kinds of Muse markers : 
% - for "one-off" markers, it verifies that there are no __START__ markers
%   associated 
% - for start/end markers, it verifies that there is the same number of
%   start and end markers 
% For the wrong markers detected, the timings are outputed by the fonction,
% and saved in a table in the directory : fullfile(cfg.datasavedir, 'wrong_markers')
% 
% Input : 
% - cfg.verifymarkers.one_off : list of "one-off" markers, ie {'marker1', 'marker2', 'marker3'}
% - cfg.verifymarkers.start   : list of "start" markers, ie {'marker1', 'marker2_start', 'marker3__START__'}
% - cfg.verifymarkers.end     : list of "end" markers, ie {'marker1', 'marker2_end', 'marker3__END__'}
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
%

cfg.verifymarkers         = ft_getopt(cfg, 'verifymarkers', []);
cfg.verifymarkers.one_off = ft_getopt(cfg.verifymarkers, 'one_off', []);
cfg.verifymarkers.start   = ft_getopt(cfg.verifymarkers, 'start', []);
cfg.verifymarkers.end     = ft_getopt(cfg.verifymarkers, 'end', []);

for ipart = 1:size(MuseStruct, 2)
    for idir = 1:length(MuseStruct{ipart})
        
        wrong_markers{ipart}{idir}           = table;
        wrong_markers{ipart}{idir}.prefix    = cfg.prefix(1:end-1);
        wrong_markers{ipart}{idir}.directory = cfg.directorylist{ipart}{idir};
        
        %% one off markers
        for markername = string(cfg.verifymarkers.one_off)
            markername_wrong = sprintf('%s__START__', markername);
            if isfield(MuseStruct{ipart}{idir}.markers,char(markername_wrong))
                wrong_markers{ipart}{idir}.(markername_wrong) = MuseStruct{ipart}{idir}.markers.(markername_wrong).synctime;
            end
        end
                
        %% start/end markers
        if length(cfg.verifymarkers.start) ~= length(cfg.verifymarkers.end)
            error('There should have the same number of markers in cfg.verifymarkers.start and in cfg.verifymarkers.end');
        end
        for imarker = 1:size(cfg.verifymarkers.start, 2)
            
            MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker})          = ft_getopt(MuseStruct{ipart}{idir}.markers, cfg.verifymarkers.start{imarker}, []);
            MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker})            = ft_getopt(MuseStruct{ipart}{idir}.markers, cfg.verifymarkers.end{imarker}, []);
            
            MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime = ft_getopt(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}), 'synctime', []);
            MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime   = ft_getopt(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}), 'synctime', []);
            
            if length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime) ~= length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime)
                wrong_markers{ipart}{idir}.(cfg.verifymarkers.start{imarker}) = sprintf('%d start markers', length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.start{imarker}).synctime));
                wrong_markers{ipart}{idir}.(cfg.verifymarkers.end{imarker})   = sprintf('%d end markers', length(MuseStruct{ipart}{idir}.markers.(cfg.verifymarkers.end{imarker}).synctime));
            end
        end
        
        %% save results
        
        if size(wrong_markers{ipart}{idir}, 2) > 2
            ft_warning('Wrong marker(s) found in %s, %s.\n', cfg.prefix(1:end-1), cfg.directorylist{ipart}{idir});
        else
            wrong_markers{ipart}{idir}.markers = 'all_markers_are_good';
        end
        
        fname = fullfile(cfg.datasavedir, 'wrong_markers', sprintf('%sp%d_dir%d_wrong_markers.xlsx', cfg.prefix, ipart, idir));
        isdir_or_mkdir(fileparts(fname));
        delete(fname);
        writetable(wrong_markers{ipart}{idir}, fname);
        
    end
end