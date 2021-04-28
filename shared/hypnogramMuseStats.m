function [marker, hypnogram, stats] = hypnogramMuseStats(cfg, MuseStruct, force)

% HYPNOGRAMSTATS creates statistics based on hypnogram and markers in MuseStruct
%
% use as
%   [MuseStruct, marker, hypnogram] = hypnogramStats(cfg, MuseStruct, force)

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

warning('off', 'all');

fname = fullfile(cfg.datasavedir, sprintf('%shypnogramStats.mat', cfg.prefix));

if nargin == 1 || force == false
    if exist(fname, 'file')
        fprintf('Reading %s\n', fname);
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname, 'stats', 'marker', 'hypnogram');
            catch ME
                err_count = err_count + 1;
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        return
    end
end
 

fprintf('Calculating stats\n');

% eventnr will increased over all marker events and all files
eventnr = 0;
marker  = table;

% include NO_SCORE
hyplabels = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"];

% excluse NO_SCORE
% hyplabels = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE"];

for markername = string(cfg.hyp.markers)

    % Go through different parts
    for ipart = 1 : size(cfg.directorylist, 2)

        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart}, 2)

            try
                StartRecord(ipart, idir) = MuseStruct{ipart}{idir}.markers.StartRecord.clock;
                StopRecord(ipart, idir)  = MuseStruct{ipart}{idir}.markers.StopRecord.clock;
            catch
            end

            if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(markername).synctime)
                continue
            end

            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(markername).synctime, 2)

                eventnr = eventnr + 1;
                marker.clock(eventnr)       = MuseStruct{ipart}{idir}.markers.(markername).clock(ievent);
                marker.name(eventnr)        = markername;
                marker.ipart(eventnr)       = ipart;
                marker.idir(eventnr)        = idir;

                % find overlap with hypnogram markers
                for hyplabel = hyplabels
                    if ~isfield(MuseStruct{ipart}{idir}.markers, [cell2mat(hyplabel), '__START__'])
                        continue
                    end
                    for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__START__']).synctime, 2)
                        x  = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent);
                        y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__START__']).synctime(i);
                        y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__END__']).synctime(i);
                        if (y1 < x) && (x < y2)
                            fprintf('Found "%s" in part %d, directory %d, overlapping with "%s" \n', markername, ipart, idir, cell2mat(hyplabel));
                            marker.hyplabel(eventnr) = hyplabel;
                        end
                    end
                end
            end
        end
    end
end

% try
%     marker.hyplabel(marker.hyplabel == "NO_SCORE")  = "AWAKE";
% catch
% end
% try
%     marker.hyplabel(ismissing(marker.hyplabel))     = "AWAKE";
% catch
% end

% find overlap of LFPs with hypnogram markers
hypnogram = table;
ihyp = 0;

% Go through different parts
for ipart = 1 : size(cfg.directorylist, 2)

    % Go through directory list
    for idir = 1 : size(cfg.directorylist{ipart}, 2)

        for hyplabel = hyplabels
            if isfield(MuseStruct{ipart}{idir}.markers, [cell2mat(hyplabel), '__START__'])
                for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__START__']).clock, 2)
                    ihyp = ihyp + 1;
                    hypnogram.starttime(ihyp)   = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__START__']).clock(i);
                    hypnogram.endtime(ihyp)     = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel), '__END__']).clock(i);
                    hypnogram.duration(ihyp)    = hypnogram.endtime(ihyp) - hypnogram.starttime(ihyp);
                    hypnogram.part(ihyp)        = ipart;
                    hypnogram.directory(ihyp)   = {MuseStruct{ipart}{idir}.directory};
                    hypnogram.hyplabel(ihyp)    = string(hyplabel);
                end
            end
        end
    end
end
%
% try
%     hypnogram.hyplabel(hypnogram.hyplabel == "NO_SCORE") = "AWAKE";
% catch
% end
hypnogram = sortrows(hypnogram);

% IEDs x Sleep stage per marker, separate for every night
for markername = string(cfg.hyp.markers)

    clear totaldur totalsum
    for ipart = unique(hypnogram.part)'

        s = 0;
        for hyplabel = hyplabels
            stats{ipart}.(markername).duration.(hyplabel)   = hours(sum(hypnogram.duration(hypnogram.hyplabel == hyplabel & hypnogram.part == ipart)));
            stats{ipart}.(markername).sum.(hyplabel)        =  sum(marker.hyplabel == hyplabel & strcmp(marker.name, markername) & marker.ipart == ipart);
%             stats{ipart}.(markername).avg.(hyplabel)        = mean(marker.hyplabel == hyplabel & strcmp(marker.name, markername) & marker.ipart == ipart);
            stats{ipart}.(markername).std.(hyplabel)        =  std(marker.hyplabel == hyplabel & strcmp(marker.name, markername) & marker.ipart == ipart);
            s = s + stats{ipart}.(markername).sum.(hyplabel);
        end
        for hyplabel = hyplabels
            if s == 0
                stats{ipart}.(markername).IEDrate.(hyplabel) = nan;
            else
                stats{ipart}.(markername).IEDrate.(hyplabel) = stats{ipart}.(markername).sum.(hyplabel) / stats{ipart}.(markername).duration.(hyplabel);
            end
        end
        for hyplabel = hyplabels
            stats{ipart}.(markername).IEDrateNorm.(hyplabel) = stats{ipart}.(markername).IEDrate.(hyplabel) / stats{ipart}.(markername).IEDrate.AWAKE;
        end
    end
end

% save results
save(fname, 'stats', 'marker', 'hypnogram');
