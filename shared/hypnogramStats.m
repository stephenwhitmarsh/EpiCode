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

fname = fullfile(cfg.datasavedir,sprintf('%shypnogramStats.mat',cfg.prefix));

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
markerlist = [];
for i = 1 : size(cfg.name,2)
    if ismember(cfg.name{i},cfg.hyp.markers)
        markerlist = [markerlist, i];
    end
end

% eventnr will increased over all marker events and all files
eventnr = 0;

marker = table;

for imarker = markerlist
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            try
                StartRecord(ipart,idir) = MuseStruct{ipart}{idir}.markers.StartRecord.clock;
                StopRecord(ipart,idir)  = MuseStruct{ipart}{idir}.markers.StopRecord.clock;
            catch
            end
            
            if ~isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                continue
            end
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)
                
                eventnr = eventnr + 1;
                marker.stage(eventnr)  = -2;
                marker.clock(eventnr)  = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(ievent);
                marker.name{eventnr}   = cfg.muse.startend{imarker,1};
                marker.ipart(eventnr)  = ipart;
                marker.idir(eventnr)   = idir;
                
                % find overlap with hypnogram markers
                for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                    if ~isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                        continue
                    end
                    for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                        x1 = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent);
                        x2 = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(ievent);
                        y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i);
                        y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i);
                        
                        if (y1 < x1) && (x1 < y2)
                            
                            fprintf('Found "%s" in part %d, directory %d, overlapping with "%s" \n',cfg.name{imarker},ipart,idir,cell2mat(hyplabel));
                            switch cell2mat(hyplabel)
                                case 'PHASE_1'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 1;
                                    marker.stage(eventnr)  = 1;
                                case 'PHASE_2'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 2;
                                    marker.stage(eventnr)  = 2;
                                case 'PHASE_3'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 3;
                                    marker.stage(eventnr)  = 3;
                                case 'REM'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 4;
                                    marker.stage(eventnr)  = 4;
                                case 'AWAKE'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 0;
                                    marker.stage(eventnr)  = 0;
                                case 'NO_SCORE'
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 0;
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

hypnogram = table;
ihyp = 0;

% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)
    
    % Go through directory list
    for idir = 1 : size(cfg.directorylist{ipart},2)
        
        % find overlap with hypnogram markers
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
                            hypnogram.stage(ihyp)  = 0;
                            hypnogram.stagelabel{ihyp} = 'AWAKE';
                        otherwise
                            error('Unexpected label name in Hypnogram\n');
                    end
                end
            end
        end
        
    end
end

% plotting
hypnogram = sortrows(hypnogram);

for imarker = markerlist
    
    clear totaldur totalsum
    for ipart = unique(hypnogram.part)'
        for i = 0 : 4
            totaldur(ipart,i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i & hypnogram.part == ipart)));
        end
        for i = 0 : 4
            totalsum(ipart,i+1) = sum(marker.stage == i & strcmp(marker.name, cfg.name{imarker}) & marker.ipart == ipart);
        end  
    end
    if sum(sum(totalsum)) == 0
        continue
    end
    IEDrate = totalsum ./ totaldur * 60;
    
    h = figure;
    
    % Number of IEDs per night / stage
    subplot(4,1,1);
    hb = bar([0 : 4],totalsum,1);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,totalsum(ib,:)', num2str(totalsum(ib,:)'),'vert','bottom','horiz','center');
    end
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'});
    ylabel('Number of IEDs');
    title('Number of IEDs per sleep stage');
    box off
    legend({'Night 1', 'Night 2', 'Night 3'},'location','eastoutside');
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
        %XData property is the tick labels/group centers; XOffset is the offset
        %of each distinct group
        xData = hb(ib).XData+hb(ib).XOffset;
        text(xData,IEDrateNorm(ib,:)',num2str(round(IEDrateNorm(ib,:),1)'),'vert','bottom','horiz','center');
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
    %         set(gca,'Yticklabels',[])
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    
    % print to file
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    set(h,'Renderer','Painters');
    print(h, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{imarker},'_x_sleep_stage.pdf']),'-r600');
    print(h, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{imarker},'_x_sleep_stage.png']),'-r600');
    
end

% save results
save(fname, 'marker', 'hypnogram');
