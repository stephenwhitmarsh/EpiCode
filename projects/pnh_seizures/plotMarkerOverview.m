function plotMarkerOverview(MuseStruct)

%% concatinate marker directories

MuseStruct_append = rmfield(MuseStruct{1},'filenames');

% for first directory
mrknames = fieldnames(MuseStruct{1}.markers);
for imarker = 1 : numel(mrknames)
    MuseStruct_append.markers.(mrknames{imarker}).sample = [];
    if isfield(MuseStruct{1}.markers.(mrknames{imarker}).events,'sample')
        for ievent = 1 : length(MuseStruct{1}.markers.(mrknames{imarker}).events)
            MuseStruct_append.markers.(mrknames{imarker}).sample = ...
                [MuseStruct_append.markers.(mrknames{imarker}).sample, ...
                MuseStruct{1}.markers.(mrknames{imarker}).events(ievent).sample];
        end
    end
    
    MuseStruct_append.markers.(mrknames{imarker}).sec = [];
    if isfield(MuseStruct{1}.markers.(mrknames{imarker}).events,'time')
        for ievent = 1 : length(MuseStruct{1}.markers.(mrknames{imarker}).events)
            MuseStruct_append.markers.(mrknames{imarker}).sec = ...
                [MuseStruct_append.markers.(mrknames{imarker}).sec, ...
                MuseStruct{1}.markers.(mrknames{imarker}).events(ievent).time];
        end
    end
end

% adding next directories
for idir = 2 : size(MuseStruct,2)
    try % some might be empty
        mrknames = fieldnames(MuseStruct{idir}.markers);
        for imarker = 1 : numel(mrknames)
            if ~isfield(MuseStruct_append.markers,(mrknames{imarker}))
                MuseStruct_append.markers.(mrknames{imarker}) = MuseStruct{idir}.markers.(mrknames{imarker});
            end
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'clock')
                MuseStruct_append.markers.(mrknames{imarker}).clock = [];
            end
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}),'clock')
                MuseStruct_append.markers.(mrknames{imarker}).clock = ...
                    [MuseStruct_append.markers.(mrknames{imarker}).clock, ...
                    MuseStruct{idir}.markers.(mrknames{imarker}).clock];
            end
            
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'samples')
                MuseStruct_append.markers.(mrknames{imarker}).sample = [];
            end
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}).events,'sample')
                for ievent = 1 : length(MuseStruct{idir}.markers.(mrknames{imarker}).events)
                    MuseStruct_append.markers.(mrknames{imarker}).sample = ...
                        [MuseStruct_append.markers.(mrknames{imarker}).sample, ...
                        MuseStruct{idir}.markers.(mrknames{imarker}).events(ievent).sample];
                end
            end
            
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'sec')
                MuseStruct_append.markers.(mrknames{imarker}).sec = [];
            end
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}).events,'time')
                for ievent = 1 : length(MuseStruct{idir}.markers.(mrknames{imarker}).events)
                    MuseStruct_append.markers.(mrknames{imarker}).sec = ...
                        [MuseStruct_append.markers.(mrknames{imarker}).sec, ...
                        MuseStruct{idir}.markers.(mrknames{imarker}).events(ievent).time];
                end
            end
            
            
        end
    catch
    end
end

%% plotting

h = figure;

subplot(8,1,1); hold;
fill([MuseStruct_append.markers.StartRecord.clock(1),MuseStruct_append.markers.StopRecord.clock(end),MuseStruct_append.markers.StopRecord.clock(end),MuseStruct_append.markers.StartRecord.clock(1)],[0 0 1 1],[1 0 0],'EdgeColor','none','facealpha',1);
for i = 1 : length(MuseStruct_append.markers.StartRecord.clock)
    fill([MuseStruct_append.markers.StartRecord.clock(i),MuseStruct_append.markers.StopRecord.clock(i),MuseStruct_append.markers.StopRecord.clock(i),MuseStruct_append.markers.StartRecord.clock(i)],[0 0 1 1],[0 1 0],'EdgeColor',[0 0 1],'facealpha',1);
end
title('Data');
axis tight
axx = xlim;
xlim(axx)

subplot(8,1,2); hold;
for i = 1 : length(MuseStruct_append.markers.BAD__START__.clock)
    fill([MuseStruct_append.markers.BAD__START__.clock(i),MuseStruct_append.markers.BAD__END__.clock(i),MuseStruct_append.markers.BAD__END__.clock(i),MuseStruct_append.markers.BAD__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('Artifacts');
xlim(axx)

subplot(8,1,3); hold;
for i = 1 : length(MuseStruct_append.markers.RR__START__.clock)
    fill([MuseStruct_append.markers.RR__START__.clock(i),MuseStruct_append.markers.RR__END__.clock(i),MuseStruct_append.markers.RR__END__.clock(i),MuseStruct_append.markers.RR__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('Fast ripples');
xlim(axx)

subplot(8,1,4); hold;
for i = 1 : length(MuseStruct_append.markers.VF1__START__.clock)
    fill([MuseStruct_append.markers.VF1__START__.clock(i),MuseStruct_append.markers.VF1__END__.clock(i),MuseStruct_append.markers.VF1__END__.clock(i),MuseStruct_append.markers.VF1__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF1');
xlim(axx)

subplot(8,1,5); hold;
for i = 1 : length(MuseStruct_append.markers.VF3__START__.clock)
    fill([MuseStruct_append.markers.VF3__START__.clock(i),MuseStruct_append.markers.VF3__END__.clock(i),MuseStruct_append.markers.VF3__END__.clock(i),MuseStruct_append.markers.VF3__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF3');
xlim(axx)

subplot(8,1,6); hold;
for i = 1 : length(MuseStruct_append.markers.VF4__START__.clock)
    fill([MuseStruct_append.markers.VF4__START__.clock(i),MuseStruct_append.markers.VF4__END__.clock(i),MuseStruct_append.markers.VF4__END__.clock(i),MuseStruct_append.markers.VF4__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF4');
xlim(axx)

subplot(8,1,7); hold;
for i = 1 : length(MuseStruct_append.markers.PP__START__.clock)
    fill([MuseStruct_append.markers.PP__START__.clock(i),MuseStruct_append.markers.PP__END__.clock(i),MuseStruct_append.markers.PP__END__.clock(i),MuseStruct_append.markers.PP__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('PP');
xlim(axx)

subplot(8,1,8); hold;

[P.minutes, P.edges]    = discretize(MuseStruct_append.markers.P.clock,'minute');
[P.histcount,edges,bin] = histcounts(P.minutes,-0.5:1:P.minutes(end)+0.5);

bar(P.edges,P.histcount,1,'facecolor','k')
ylabel('Spikes per minute');
axis tight;
xlim(axx)

c = smooth(P.histcount);
plot(P.edges,c,'r','linewidth',1)

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'overview_visual.pdf','-r600');
