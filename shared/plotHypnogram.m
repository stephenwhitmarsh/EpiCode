function [t] = plotHypnogram(cfg, MuseStruct)

% PLOTHYPNOGRAM creates statistics based on hypnogram and markers in MuseStruct
%
% use as
%   [t] = plotHypnogram(cfg, MuseStruct)

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

cfg.visible = ft_getopt(cfg, 'visible', 'on');
hyplabels = {'NO_SCORE__START__', 'AWAKE__START__', 'PHASE_1__START__', 'PHASE_2__START__', 'PHASE_3__START__', 'REM__START__', 'NO_SCORE__END__', 'AWAKE__END__', 'PHASE_1__END__', 'PHASE_2__END__', 'PHASE_3__END__', 'REM__END__'};

% loop over different parts, i.e. lists of directories
for ipart = 1 : size(MuseStruct,2)

    % remove empty structs
    MuseStruct{ipart} = MuseStruct{ipart}(~cellfun('isempty', MuseStruct{ipart}));

    % will first append all markers over directories
    MuseStruct_append{ipart}            = [];
    MuseStruct_append{ipart}.markers    = [];

    % loop over directories
    for idir = 1 : size(MuseStruct{ipart},2)

        fprintf('Working on directory %d of %d, part %d\n', idir, size(MuseStruct{ipart},2), ipart);
        for markername = [string(cfg.hyp.markers), hyplabels]

            % if marker field doesn't exist yet, create it
            if ~isfield(MuseStruct_append{ipart}.markers, markername)
                MuseStruct_append{ipart}.markers.(markername) = [];
            end

            % if marker.clock field doesn't exist yet, create it
            if ~isfield(MuseStruct_append{ipart}.markers.(markername), 'clock')
                MuseStruct_append{ipart}.markers.(markername).clock = [];
            end

            % add directory
            if ~isfield(MuseStruct_append{ipart}.markers.(markername), 'directory')
                MuseStruct_append{ipart}.markers.(markername).directory = [];
            end

            if isfield(MuseStruct{ipart}{idir}.markers, markername)

                % append clock field
                if isfield(MuseStruct{ipart}{idir}.markers.(markername), 'clock')
                    MuseStruct_append{ipart}.markers.(markername).clock = ...
                        [MuseStruct_append{ipart}.markers.(markername).clock, ...
                        MuseStruct{ipart}{idir}.markers.(markername).clock];

                    % add directory
                    MuseStruct_append{ipart}.markers.(markername).directory = ...
                        [MuseStruct_append{ipart}.markers.(markername).directory; ...
                        repmat({MuseStruct{ipart}{idir}.directory}, size(MuseStruct{ipart}{idir}.markers.(markername).clock,2), 1)];

                    fprintf('Concatinated %d markers: %s\n', size(MuseStruct{ipart}{idir}.markers.(markername).clock, 2), markername);
                end
            end
        end
    end

    %     % remove empty markers
    %     fn = fieldnames(MuseStruct_append{ipart}.markers);
    %     for imarker = 1 : numel(fn)
    %         if ~isfield(MuseStruct_append{ipart}.markers.(fn{imarker}), 'clock')
    %             MuseStruct_append{ipart}.markers = rmfield(MuseStruct_append{ipart}.markers, fn{imarker});
    %         else
    %             if isempty(MuseStruct_append{ipart}.markers.(fn{imarker}).clock)
    %                 MuseStruct_append{ipart}.markers = rmfield(MuseStruct_append{ipart}.markers,fn{imarker});
    %             end
    %         end
    %     end

    % concatinate markers
    fn = fieldnames(MuseStruct_append{ipart}.markers);
    markerlabel     = [];
    starttime       = [];
    startsample     = [];
    directory       = [];
    for imarker = 1 : numel(fn)
        markerlabel     = [markerlabel;     repmat(convertCharsToStrings(fn{imarker}),numel(MuseStruct_append{ipart}.markers.(fn{imarker}).clock),1)];
        starttime       = [starttime;       MuseStruct_append{ipart}.markers.(fn{imarker}).clock'];
        directory       = [directory;       MuseStruct_append{ipart}.markers.(fn{imarker}).directory];
    end

    endtime = starttime; % same time if cant find end marker (below)

    t = table(markerlabel,starttime,endtime,directory);
    t = unique(t,'rows');
    t = sortrows(t,2);

    % find corresponding end-times of markers; code could be improved
    i = 1;
    while i < height(t)
        if contains(t.markerlabel(i),'__START__')
            indx                = find(contains(t.markerlabel,strcat(t.markerlabel{i}(1:end-9),'__END__')));
            endindx             = find(t.starttime(indx) > t.starttime(i),1,'first');
            if isempty(endindx)
                disp('sdf')
            end
            t.markerlabel(i)    = t.markerlabel{i}(1:end-9);
            t.endtime(i)        = t.starttime(indx(endindx));
            t.duration(i)       = t.endtime(i) - t.starttime(indx(endindx));
            t(indx(endindx),:)  = [];
        end
        i = i + 1;
    end

    % select hypnogram labels to plot
    hyp_tbl = t(contains(t.markerlabel, {'NO_SCORE', 'AWAKE', 'PHASE_1', 'PHASE_2', 'PHASE_3', 'REM'}),:);

    % select markers to plot
    mrk_tbl = t(any(t.markerlabel == string(cfg.hyp.markers), 2),:);

    if isempty(hyp_tbl)
        fprintf('No hypnogram found in part %d\n', ipart)
        continue
    end

    % alternative selection of start/end of hypnogram
    % startsindx = find(contains(t.markerlabel,'CriseStart'));
    % endsindx   = find(contains(t.markerlabel,'CriseEnd'));
    % startsindx = find(contains(t.markerlabel,cfg.pattern.startmarker));
    % endsindx   = find(contains(t.markerlabel,cfg.pattern.endmarker));
    % fprintf('\n%d Patterns found\n',numel(startsindx));

    % segment into patterns/seizures/hypnograms,
    % if there is more than 4 hours in between them
    hyp_startsindx  = [1; find(diff(hyp_tbl.endtime) > hours(4))+1];
    hyp_endsindx    = [find(diff(hyp_tbl.endtime) > hours(4)); height(hyp_tbl)];

    % select markers that occur within hypnogram
    hyp_starttime   = hyp_tbl.starttime(hyp_startsindx);
    hyp_endtime     = hyp_tbl.starttime(hyp_endsindx);
    mrk_night       = mrk_tbl(mrk_tbl.starttime >= hyp_starttime & mrk_tbl.endtime <= hyp_endtime,:);
    maxlength       = max(hyp_endtime - hyp_starttime);

    %% plotting
    fig = figure('visible', cfg.visible);
    subplot(numel(unique(mrk_tbl.markerlabel))+1,1,1); hold;

    fill([hyp_starttime, hyp_starttime + maxlength, hyp_starttime + maxlength, hyp_starttime], [0 0 1 1], [1 1 1], 'EdgeColor',[1 1 1]);
    X = [];
    Y = [];
    for im = 1 : height(hyp_tbl)
        if ~isempty(X)
            % if there's a gap, 'fill' with 0
            if hyp_tbl.starttime(im) ~= X(end)
                X = [X, X(end) hyp_tbl.starttime(im)];
                Y = [Y, 0,  0];
            end
        end
        X = [X, hyp_tbl.starttime(im), hyp_tbl.endtime(im)];

        % height in hypnogram is based on order of cfg.hyp.contains
        switch cell2mat(hyp_tbl.markerlabel(im))
            case 'NO_SCORE'
                y = 5;
            case 'AWAKE'
                y = 5;
            case 'REM'
                y = 4;
            case 'PHASE_1'
                y = 3;
            case 'PHASE_2'
                y = 2;
            case 'PHASE_3'
                y = 1;
        end
        Y = [Y, y, y];
    end

    for i = 1 : length(X)-1
        if Y(i) ~= 0 && Y(i+1) ~= 0
            if Y(i) == 4 && Y(i+1) == 4 % REM gets thicker line
                plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
            else
                plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
            end
        end
    end
    set(gca,'Layer','top');
    set(gca,'Ytick', 1 : 5, 'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'},'TickDir','out');
    axis tight;
    xl = xlim;

    % plot markers
    iplot = 2;
    for markername = unique(mrk_tbl.markerlabel)'
        subplot(numel(unique(mrk_tbl.markerlabel))+1, 1, iplot); hold;
        title(markername);
        sel = mrk_tbl.markerlabel == markername;
        histogram(mrk_tbl.starttime(sel), 'BinWidth', seconds(60), 'EdgeColor', 'none', 'FaceColor', [0, 0, 0])
        %         for im = sel'
        %             fill([mrk_tbl.starttime(im), mrk_tbl.endtime(im), mrk_tbl.endtime(im), mrk_tbl.starttime(im)], [0 0 1 1], [0 0 0]);
        %         end
        xlim(xl);
        iplot = iplot + 1;
        set(gca,'Xticklabels', []);
        ylabel('IEDs per minute');
    end

    % print to file
    set(fig,'PaperOrientation', 'landscape');
    set(fig,'PaperUnits', 'normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    set(fig,'Renderer', 'Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix, 'p', num2str(ipart), '_hypnogram.pdf']), '-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix, 'p', num2str(ipart), '_hypnogram.png']), '-r600');

    %     % print a 3 meter version for fun
    %     set(gcf,'PaperUnits','centimeters');
    %     set(gcf,'PaperSize',[300 10]);
    %     h.PaperUnits = 'centimeters';
    %     h.PaperPosition = [0 0 300 10];
    %     h.Units = 'centimeters';
    %     h.PaperSize=[300 10];
    %     h.Units = 'centimeters';
    %     print(h, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'part',num2str(ipart),'-hypnogram_3m.pdf']),'-r600');

    writetable(t,fullfile(cfg.datasavedir,[cfg.prefix, 'p', num2str(ipart), '_hypnogram.txt']));

end
disp('Done');
