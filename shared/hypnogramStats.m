function [marker, hypnogram] = hypnogramStats(cfg, MuseStruct, force)

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

warning('off','all');

cfg.visible = ft_getopt(cfg, 'visible', 'on');

fname = fullfile(cfg.datasavedir, sprintf('%shypnogramStats.mat',cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('**********************************\n');
    fprintf('** loading precomputed stats *****\n');
    fprintf('**********************************\n\n');
    load(fname, 'marker', 'hypnogram');
    return
end

fprintf('***********************\n');
fprintf('** creating stats *****\n');
fprintf('***********************\n\n');

% elect those markers to load
% markerlist = [];
% for i = 1 : size(cfg.name,2)
%     if ismember(cfg.name{i},cfg.hyp.markers)
%         markerlist = [markerlist, i];
%     end
% end

% eventnr will increased over all marker events and all files
eventnr = 0;
marker  = table;

for markername = string(cfg.hyp.markers)
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            try
                StartRecord(ipart,idir) = MuseStruct{ipart}{idir}.markers.StartRecord.clock;
                StopRecord(ipart,idir)  = MuseStruct{ipart}{idir}.markers.StopRecord.clock;
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
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(markername).synctime,2)
                
                eventnr = eventnr + 1;
                marker.stage(eventnr)  = -2;
                marker.clock(eventnr)  = MuseStruct{ipart}{idir}.markers.(markername).clock(ievent);
                marker.name{eventnr}   = markername;
                marker.ipart(eventnr)  = ipart;
                marker.idir(eventnr)   = idir;
                
                % find overlap with hypnogram markers
                for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                    if ~isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                        continue
                    end
                    for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                        x1 = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent);
                        x2 = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(ievent);
                        y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i);
                        y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i);
                        
                        if (y1 < x1) && (x1 < y2)
                            
                            fprintf('Found "%s" in part %d, directory %d, overlapping with "%s" \n', markername, ipart, idir, cell2mat(hyplabel));
                            switch cell2mat(hyplabel)
                                case 'PHASE_1'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 1;
                                    marker.stage(eventnr)  = 1;
                                case 'PHASE_2'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 2;
                                    marker.stage(eventnr)  = 2;
                                case 'PHASE_3'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 3;
                                    marker.stage(eventnr)  = 3;
                                case 'REM'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 4;
                                    marker.stage(eventnr)  = 4;
                                case 'AWAKE'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 0;
                                    marker.stage(eventnr)  = 0;
                                case 'NO_SCORE'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).stage = 0;
                                    marker.stage(eventnr)  = 0;
                                otherwise
                                    error('Unexpected label name in Hypnogram\n');
                            end
                        end
                    end
                end
            end
        end
    end
end

% find overlap of LFPs with hypnogram markers
hypnogram = table;
ihyp = 0;

% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)
    
    % Go through directory list
    for idir = 1 : size(cfg.directorylist{ipart},2)
        
        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
            if isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).clock,2)
                    ihyp = ihyp + 1;
                    hypnogram.starttime(ihyp)   = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).clock(i);
                    hypnogram.endtime(ihyp)     = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).clock(i);
                    hypnogram.duration(ihyp)    = hypnogram.endtime(ihyp) - hypnogram.starttime(ihyp);
                    hypnogram.part(ihyp)        = ipart;
                    hypnogram.directory(ihyp)   = {MuseStruct{ipart}{idir}.directory};
                    switch cell2mat(hyplabel)
                        case 'PHASE_1'
                            hypnogram.stage(ihyp)  = 1;
                            hypnogram.stagelabel{ihyp} = 'STAGE 1';
                        case 'PHASE_2'
                            hypnogram.stage(ihyp)  = 2;
                            hypnogram.stagelabel{ihyp} = 'STAGE 2';
                        case 'PHASE_3'
                            hypnogram.stage(ihyp)  = 3;
                            hypnogram.stagelabel{ihyp} = 'STAGE 3';
                        case 'REM'
                            hypnogram.stage(ihyp)  = 4;
                            hypnogram.stagelabel{ihyp} = 'REM';
                        case 'AWAKE'
                            hypnogram.stage(ihyp)  = 0;
                            hypnogram.stagelabel{ihyp} = 'AWAKE';
                        case 'NO_SCORE'
                            hypnogram.stage(ihyp)  = -1;
                            hypnogram.stagelabel{ihyp} = 'AWAKE';
                        otherwise
                            error('Unexpected label name in Hypnogram\n');
                    end
                end
            end
        end
        
    end
end

hypnogram = sortrows(hypnogram);

%% plotting IEDs x Sleep stage per template, SEPARATE for every night
for markername = string(cfg.hyp.markers)
    
    clear totaldur totalsum
    for ipart = unique(hypnogram.part)'
        for i = 0 : 4
            totaldur(ipart,i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i & hypnogram.part == ipart)));
        end
        for i = 0 : 4
            totalsum(ipart,i+1) = sum(marker.stage == i & strcmp(marker.name, markername) & marker.ipart == ipart);
        end
    end
    if sum(sum(totalsum)) == 0
        continue
    end
    IEDrate = totalsum ./ totaldur * 60;
    
    fig = figure('visible', cfg.visible);
    
    % Number of IEDs per night / stage
    subplot(4,1,1);
    hb = bar(0:4, totalsum,1);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,totalsum(ib,:)', num2str(totalsum(ib,:)'),'vert', 'bottom', 'horiz', 'center');
    end
    xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'});
    ylabel('Number of IEDs');
    title('Number of IEDs per sleep stage');
    box off
    legend({'Night 1', 'Night 2', 'Night 3'}, 'location', 'eastoutside');
    if max(max(totalsum)) > 0
        ylim([0, max(max(totalsum))*1.3]);
    end
    set(gca,'TickLength',[0 0]); % set(gca,'Yticklabels',[])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % hours spend in sleep stage
    subplot(4,1,2);
    hb = bar(0 : 4,totaldur/60/60,1);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,totaldur(ib,:)'/60/60, num2str(round(totaldur(ib,:)'/60/60,1)),'vert','bottom','horiz','center');
    end
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('Total duration per sleep stage (hrs)');
    ylabel('Hours');
    box off
    legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
    ylim([0, max(max(totaldur)/60/60)*1.3]);
    set(gca,'TickLength',[0 0])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % number of IEDs per minute
    subplot(4,1,3);
    hb = bar(0 : 4, IEDrate, 1);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,IEDrate(ib,:)',num2str(round(IEDrate(ib,:),1)'),'vert','bottom','horiz','center');
    end
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('IEDs per minute');
    ylabel('IEDs per minute');
    box off
    legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
    if max(max(IEDrate)) > 0
        ylim([0, max(max(IEDrate))*1.3]);
    end
    set(gca,'TickLength',[0 0])
    %         set(gca,'Yticklabels',[])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % IED rate normalized to wake
    IEDrateNorm = IEDrate ./ IEDrate(:,1);
    subplot(4,1,4);
    hb = bar(0 : 4, IEDrateNorm, 1);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,IEDrateNorm(ib,:)', num2str(round(IEDrateNorm(ib,:),1)'),'vert','bottom','horiz','center');
    end
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('IEDs per minute normalized to wake stage');
    ylabel('IED rate vs. wake');
    box off
    legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
    if max(max(IEDrateNorm)) > 0
        ylim([0, max(max(IEDrateNorm))*1.3]);
    end
    set(gca,'TickLength',[0 0])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % print to file
    set(fig, 'PaperOrientation','landscape');
    set(fig, 'PaperUnits','normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix,markername,'_x_sleep_stage.pdf')),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix,markername,'_x_sleep_stage.png')),'-r600');
    
    %% plotting IEDs x Sleep stage per template, COMBINED for every night
    clear totaldur totalsum
    for i = 0 : 4
        totaldur(i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i)));
    end
    for i = 0 : 4
        totalsum(i+1) = sum(marker.stage == i & strcmp(marker.name, markername));
    end
    IEDrate = totalsum ./ totaldur * 60;
    
    fig = figure('visible', cfg.visible);
    cm  = cool(5);
    x   = [2 3 4 5 1];
    
    % Number of IEDs per night / stage
    subplot(4,1,1); hold;
    for ib = 1:5
        hb = bar(x(ib), totalsum(ib), 1);
        xData = hb.XData+hb.XOffset;
        text(xData, totalsum(ib), num2str(totalsum(ib)),'vert', 'bottom', 'horiz', 'center');
        set(hb, 'FaceColor', cm(x(ib),:));
    end
    xlim([0,6]);
    xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
    ylabel('Number of IEDs');
    title('Number of IEDs per sleep stage');
    box off
    if max(totalsum) > 0
        ylim([0, max(totalsum)*1.3]);
    end
    set(gca,'TickLength',[0 0]); % set(gca,'Yticklabels',[])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % hours spend in sleep stage
    subplot(4,1,2); hold;
    for ib = 1:5
        hb = bar(x(ib), totaldur(ib) /60 /60, 1);
        xData = hb.XData + hb.XOffset;
        text(xData, totaldur(ib) /60 /60, num2str(round(totaldur(ib) /60 /60), 1),'vert', 'bottom', 'horiz', 'center');
        set(hb, 'FaceColor', cm(x(ib),:));
    end
    xlim([0,6]);
    xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
    title('Total duration per sleep stage (hrs)');
    ylabel('Hours');
    box off
    ylim([0, max(totaldur) /60 /60 * 1.3]);
    set(gca,'TickLength', [0 0])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % number of IEDs per minute
    subplot(4,1,3); hold;
    for ib = 1:5
        hb = bar(x(ib), IEDrate(ib), 1);
        xData = hb.XData + hb.XOffset;
        text(xData, IEDrate(ib), num2str(round(IEDrate(ib), 1)), 'vert', 'bottom', 'horiz', 'center');
        set(hb, 'FaceColor', cm(x(ib),:));
    end
    xlim([0,6]);
    xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
    title('IEDs per minute');
    ylabel('IEDs per minute');
    box off
    if max(IEDrate) > 0
        ylim([0, max(IEDrate)*1.3]);
    end
    set(gca,'TickLength',[0 0])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % IED rate normalized to wake
    IEDrateNorm = IEDrate ./ IEDrate(1);
    subplot(4,1,4); hold;
    for ib = 1:5
        hb = bar(x(ib), IEDrateNorm(ib), 1);
        xData = hb.XData + hb.XOffset;
        text(xData,IEDrateNorm(ib), num2str(round(IEDrateNorm(ib),1)),'vert','bottom','horiz','center');
        set(hb, 'FaceColor', cm(x(ib),:));
    end
    xlim([0,6]);
    xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
    title('IEDs per minute normalized to wake stage');
    ylabel('IED rate vs. wake');
    box off
    if max(IEDrateNorm) > 0
        ylim([0, max(IEDrateNorm) * 1.3]);
    end
    set(gca,'TickLength',[0 0])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % print to file
    set(fig, 'PaperOrientation', 'landscape');
    set(fig, 'PaperUnits', 'normalized');
    set(fig, 'PaperPosition', [0 0 1 1]);
    set(fig, 'Renderer','Painters');
    print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, markername,'_x_all_nights_x_sleep_stage.pdf')), '-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, markername,'_x_all_nights_x_sleep_stage.png')), '-r600');
    
end

%% plotting IEDs x Sleep stage for all templates) SEPARATE FOR EVERY NIGHT
clear totaldur totalsum
for ipart = unique(hypnogram.part)'
    for i = 0 : 4
        totaldur(ipart, i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i & hypnogram.part == ipart)));
    end
    for i = 0 : 4
        totalsum(ipart, i+1) = sum(marker.stage == i & marker.ipart == ipart);
    end
end
IEDrate = totalsum ./ totaldur * 60;

fig = figure('visible', cfg.visible);

% Number of IEDs per night / stage
subplot(4,1,1);
hb = bar(0:4, totalsum,1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,totalsum(ib,:)', num2str(totalsum(ib,:)'),'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'});
ylabel('Number of IEDs');
title('Number of IEDs per sleep stage');
box off
legend({'Night 1', 'Night 2', 'Night 3'}, 'location', 'eastoutside');
if max(max(totalsum)) > 0
    ylim([0, max(max(totalsum))*1.3]);
end
set(gca,'TickLength',[0 0]); % set(gca,'Yticklabels',[])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% hours spend in sleep stage
subplot(4,1,2);
hb = bar(0 : 4,totaldur/60/60,1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,totaldur(ib,:)' /60 /60, num2str(round(totaldur(ib,:)' /60 /60, 1)),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
title('Total duration per sleep stage (hrs)');
ylabel('Hours');
box off
legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
ylim([0, max(max(totaldur) /60 /60) * 1.3]);
set(gca,'TickLength', [0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% number of IEDs per minute
subplot(4,1,3);
hb = bar(0 : 4, IEDrate, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,IEDrate(ib,:)', num2str(round(IEDrate(ib,:),1)'), 'vert', 'bottom', 'horiz', 'center');
end
xticklabels({'Wake', 'Stage 1', 'Stage 2', 'Stage 3', 'REM'})
title('IEDs per minute');
ylabel('IEDs per minute');
box off
legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
if max(max(IEDrate)) > 0
    ylim([0, max(max(IEDrate))*1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% IED rate normalized to wake
IEDrateNorm = IEDrate ./ IEDrate(:,1);
subplot(4,1,4);
hb = bar(0 : 4, IEDrateNorm, 1);
for ib = 1:numel(hb)
    xData = hb(ib).XData+hb(ib).XOffset;
    text(xData,IEDrateNorm(ib,:)', num2str(round(IEDrateNorm(ib,:),1)'),'vert','bottom','horiz','center');
end
xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
title('IEDs per minute normalized to wake stage');
ylabel('IED rate vs. wake');
box off
legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
if max(max(IEDrateNorm)) > 0
    ylim([0, max(max(IEDrateNorm)) * 1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% print to file
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer','Painters');
print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'all_templates_x_sleep_stage.pdf')), '-r600');
print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'all_templates_x_sleep_stage.png')), '-r600');


%% plotting IEDs x Sleep stage for all templates) FOR ALL NIGHTS COMBINED
clear totaldur totalsum
for i = 0 : 4
    totaldur(i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i)));
end
for i = 0 : 4
    totalsum(i+1) = sum(marker.stage == i);
end
IEDrate = totalsum ./ totaldur * 60;

fig = figure('visible', cfg.visible);
cm  = cool(5);
x   = [2 3 4 5 1];

% Number of IEDs per night / stage
subplot(4,1,1); hold;
for ib = 1:5
    hb = bar(x(ib), totalsum(ib), 1);
    xData = hb.XData+hb.XOffset;
    text(xData, totalsum(ib), num2str(totalsum(ib)),'vert', 'bottom', 'horiz', 'center');
    set(hb, 'FaceColor', cm(x(ib),:));
end
xlim([0,6]);
xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
ylabel('Number of IEDs');
title('Number of IEDs per sleep stage');
box off
if max(totalsum) > 0
    ylim([0, max(totalsum)*1.3]);
end
set(gca,'TickLength',[0 0]); % set(gca,'Yticklabels',[])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% hours spend in sleep stage
subplot(4,1,2); hold;
for ib = 1:5
    hb = bar(x(ib), totaldur(ib) /60 /60, 1);
    xData = hb.XData + hb.XOffset;
    text(xData, totaldur(ib) /60 /60, num2str(round(totaldur(ib) /60 /60), 1),'vert', 'bottom', 'horiz', 'center');
    set(hb, 'FaceColor', cm(x(ib),:));
end
xlim([0,6]);
xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
title('Total duration per sleep stage (hrs)');
ylabel('Hours');
box off
ylim([0, max(totaldur) /60 /60 * 1.3]);
set(gca,'TickLength', [0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% number of IEDs per minute
subplot(4,1,3); hold;
for ib = 1:5
    hb = bar(x(ib), IEDrate(ib), 1);
    xData = hb.XData + hb.XOffset;
    text(xData, IEDrate(ib), num2str(round(IEDrate(ib), 1)), 'vert', 'bottom', 'horiz', 'center');
    set(hb, 'FaceColor', cm(x(ib),:));
end
xlim([0,6]);
xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
title('IEDs per minute');
ylabel('IEDs per minute');
box off
if max(IEDrate) > 0
    ylim([0, max(IEDrate)*1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% IED rate normalized to wake
IEDrateNorm = IEDrate ./ IEDrate(1);
subplot(4,1,4); hold;
for ib = 1:5
    hb = bar(x(ib), IEDrateNorm(ib), 1);
    xData = hb.XData + hb.XOffset;
    text(xData,IEDrateNorm(ib), num2str(round(IEDrateNorm(ib),1)),'vert','bottom','horiz','center');
    set(hb, 'FaceColor', cm(x(ib),:));
end
xlim([0,6]);
xticklabels({'', 'REM','Wake', 'Stage 1', 'Stage 2', 'Stage 3'});
title('IEDs per minute normalized to wake stage');
ylabel('IED rate vs. wake');
box off
if max(IEDrateNorm) > 0
    ylim([0, max(IEDrateNorm) * 1.3]);
end
set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% print to file
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);
set(fig, 'Renderer','Painters');
print(fig, '-dpdf', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'all_templates_x_all_nights_x_sleep_stage.pdf')), '-r600');
print(fig, '-dpng', fullfile(cfg.imagesavedir, strcat(cfg.prefix, 'all_templates_x_all_nights_x_sleep_stage.png')), '-r600');

% save results
save(fname, 'marker', 'hypnogram');
