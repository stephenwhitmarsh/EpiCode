function plotPatterns(config,MuseStruct)

% remove empty structs
MuseStruct = MuseStruct(~cellfun('isempty',MuseStruct));
MuseStruct_append = []; %rmfield(MuseStruct{1},'filenames');
MuseStruct_append.markers = [];

% adding next directories
for idir = 1 : size(MuseStruct,2)
    fprintf('Working on directory %d of %d\n',idir,size(MuseStruct,2));
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
                MuseStruct_append.markers.(mrknames{imarker}).samples = [];
            end
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}).events,'sample')
                for ievent = 1 : length(MuseStruct{idir}.markers.(mrknames{imarker}).events)
                    MuseStruct_append.markers.(mrknames{imarker}).samples = ...
                        [MuseStruct_append.markers.(mrknames{imarker}).samples, ...
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

% remove empty markers
fn = fieldnames(MuseStruct_append.markers);
for imarker = 1 : numel(fn)
    if isempty(MuseStruct_append.markers.(fn{imarker}).samples)
        MuseStruct_append.markers = rmfield(MuseStruct_append.markers,fn{imarker});
    end
end

% concatinate markers
fn = fieldnames(MuseStruct_append.markers);
markerlabel = [];
starttime = [];
for imarker = 1 : numel(fn)
    markerlabel = [markerlabel; repmat(convertCharsToStrings(fn{imarker}),numel(MuseStruct_append.markers.(fn{imarker}).clock),1)];
    starttime = [starttime; MuseStruct_append.markers.(fn{imarker}).clock'];
end

endtime = starttime; % same time if cant find end marker
t = table(markerlabel,starttime,endtime);
t = unique(t,'rows');
t = sortrows(t,2);

i = 1;
while i < height(t)
    if contains(t.markerlabel(i),'__START__')     
        t.markerlabel(i) = t.markerlabel{i}(1:end-9);
        endsindx = find(contains(t.markerlabel,strcat(t.markerlabel(i),'__END__')));
        endsindx = endsindx(endsindx > i);
        endindx = endsindx(1);
        t.endtime(i) = t.starttime(endindx);
        t(endindx,:) = [];
    end
    i = i + 1;
end

% t(contains(t.markerlabel,config.hyp.notcontains),:) = [];
t = t(contains(t.markerlabel,config.hyp.contains),:);

% 
% startsindx = find(contains(t.markerlabel,'CriseStart'));
% endsindx   = find(contains(t.markerlabel,'CriseEnd'));
% 
% startsindx = find(contains(t.markerlabel,config.pattern.startmarker));
% endsindx   = find(contains(t.markerlabel,config.pattern.endmarker));
% fprintf('\n%d Patterns found\n',numel(startsindx));

startsindx  = [1; find(diff(t.endtime) > hours(2))+1];
endsindx    = [find(diff(t.endtime) > hours(2)); length(endsindx)];


ii = 1;
for i = 1 : numel(startsindx)
    if endsindx(i) - startsindx(i) > 1
        s{ii} = t(startsindx(i)+1:endsindx(i)-1,:);
        Cstart(ii) = t.starttime(startsindx(i));
        Cend(ii)    = t.starttime(endsindx(i));
        ii = ii + 1;
    end
end
t(contains(t.markerlabel,config.pattern.startmarker) | contains(t.markerlabel,config.pattern.endmarker),:) = [];
maxlength = max(Cend - Cstart);


%% plotting
colortable.label = unique(t.markerlabel);
colortable.color = linspecer(numel(unique(t.markerlabel)));

h = figure;
for is = 1 : size(s,2)
    subplot(size(s,2)+1,1,is); hold;
    fill([s{is}.starttime(1), Cstart(is) + maxlength, Cstart(is) + maxlength, s{is}.starttime(1)],[0 0 1 1],[1 1 1],'EdgeColor',[1 1 1]);
    X = [];
    Y = [];
    for im = 1 : height(s{is})
        c = colortable.color(strcmp(colortable.label,s{is}.markerlabel(im)),:);
        if ~isempty(X)
            if s{is}.starttime(im) ~= X(end)
                X = [X, X(end) s{is}.starttime(im)];
                Y = [Y, 0 0];
            end
        end
        X = [X, s{is}.starttime(im), s{is}.endtime(im)];
        y = find(contains(config.hyp.contains,t.markerlabel(im)));
        Y = [Y, y, y];
    end
    for i = 1 : length(X)-1
        if Y(i) ~= 0 && Y(i+1) ~= 0
            if strcmp(config.hyp.contains(Y(i)),'REM') && strcmp(config.hyp.contains(Y(i+1)),'REM')
                
                plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k','LineWidth',3);
            else
                plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k');
            end
        end
    end
    set(gca,'Layer','top');
    set(gca,'Ytick', 1 : length(config.hyp.contains),'Yticklabels',strrep(config.hyp.contains,'_',' '),'TickDir','out');
    axis tight;
end
% 
% subplot(size(s,2)+1,1,size(s,2)+1); hold;
% for marker = 1 : numel(colortable.label)
%     c = colortable.color(marker,:);
%     x1 = marker*2;
%     x2 = marker*2+1.5;
%     fill([x1,x2,x2,x1],[0 0 1 1],c,'EdgeColor',c,'facealpha',1);
%     text((x1+x2)/2,0.1,colortable.label(marker),'Rotation',90,'Fontsize',10);
%     set(gca,'Layer','top');
% end
% set(gca,'XTick',[]);
% axis tight

% print to file
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
set(h,'Renderer','Painters');
print(h, '-dpdf', fullfile(config.imagesavedir,[config.prefix,'overview_patterns.pdf']),'-r600');
% writetable(t,fullfile(config.imagesavedir,'seg_labels'));
disp('Done');
