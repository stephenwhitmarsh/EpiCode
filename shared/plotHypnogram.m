function plotHypnogram(config,MuseStruct)

% loop over different parts, i.e. lists of directories
for ipart = 1 : size(MuseStruct,2)
    
    % remove empty structs
    MuseStruct{ipart} = MuseStruct{ipart}(~cellfun('isempty',MuseStruct{ipart}));
    
    % will first append all markers over directories
    MuseStruct_append{ipart}            = [];
    MuseStruct_append{ipart}.markers    = [];
    
    % loop over directories
    for idir = 1 : size(MuseStruct{ipart},2)
        fprintf('Working on directory %d of %d\n',idir,size(MuseStruct{ipart},2));
        try % some might be empty
            mrknames = fieldnames(MuseStruct{ipart}{idir}.markers);
            for imarker = 1 : numel(mrknames)
                
                % if marker field doesn't exist yet, create it
                if ~isfield(MuseStruct_append{ipart}.markers,(mrknames{imarker}))
                    MuseStruct_append{ipart}.markers.(mrknames{imarker}) = MuseStruct{ipart}{idir}.markers.(mrknames{imarker});
                end
                
                % if marker.clock field doesn't exist yet, create it
                if ~isfield(MuseStruct_append{ipart}.markers.(mrknames{imarker}),'clock')
                    MuseStruct_append{ipart}.markers.(mrknames{imarker}).clock = [];
                end
                
                % append clock field
                if isfield(MuseStruct{ipart}{idir}.markers.(mrknames{imarker}),'clock')
                    MuseStruct_append{ipart}.markers.(mrknames{imarker}).clock = ...
                        [MuseStruct_append{ipart}.markers.(mrknames{imarker}).clock, ...
                        MuseStruct{ipart}{idir}.markers.(mrknames{imarker}).clock];
                end
                
            end
        catch
        end
    end
    
    % remove empty markers
    fn = fieldnames(MuseStruct_append{ipart}.markers);
    for imarker = 1 : numel(fn)
        if ~isfield(MuseStruct_append{ipart}.markers.(fn{imarker}),'clock')
            MuseStruct_append{ipart}.markers = rmfield(MuseStruct_append{ipart}.markers,fn{imarker});
        else
            if isempty(MuseStruct_append{ipart}.markers.(fn{imarker}).clock)
                MuseStruct_append{ipart}.markers = rmfield(MuseStruct_append{ipart}.markers,fn{imarker});
            end
        end
    end
    
    % concatinate markers
    fn = fieldnames(MuseStruct_append{ipart}.markers);
    markerlabel = [];
    starttime = [];
    startsample = [];
    for imarker = 1 : numel(fn)
        markerlabel     = [markerlabel;     repmat(convertCharsToStrings(fn{imarker}),numel(MuseStruct_append{ipart}.markers.(fn{imarker}).clock),1)];
        starttime       = [starttime;       MuseStruct_append{ipart}.markers.(fn{imarker}).clock'];
    end
    
    endtime = starttime; % same time if cant find end marker (below)
    
    t = table(markerlabel,starttime,endtime);
    t = unique(t,'rows');
    t = sortrows(t,2);
    
    % find corresponding end-times of markers; code could be improved
    i = 1;
    while i < height(t)
        if contains(t.markerlabel(i),'__START__')
            t.markerlabel(i)    = t.markerlabel{i}(1:end-9);
            endsindx            = find(contains(t.markerlabel,strcat(t.markerlabel{i},'__END__')));
            endsindx            = endsindx(endsindx > i);
            endindx             = endsindx(1);
            t.endtime(i)        = t.starttime(endindx);
            t(endindx,:)        = [];
        end
        i = i + 1;
    end
    
    % select hypnogram labels to plot
    hyp_tbl = t(contains(t.markerlabel,config.hyp.contains),:);
    mrk_tbl = t(contains(t.markerlabel,config.hyp.markers),:);
    
    % alternative selection of start/end of hypnogram
    % startsindx = find(contains(t.markerlabel,'CriseStart'));
    % endsindx   = find(contains(t.markerlabel,'CriseEnd'));
    % startsindx = find(contains(t.markerlabel,config.pattern.startmarker));
    % endsindx   = find(contains(t.markerlabel,config.pattern.endmarker));
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

    h = figure;
    subplot(numel(unique(mrk_tbl.markerlabel))+1,1,1); hold;
    
    fill([hyp_starttime, hyp_starttime + maxlength, hyp_starttime + maxlength,hyp_starttime],[0 0 1 1],[1 1 1],'EdgeColor',[1 1 1]);
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
        
        % height in hypnogram is based on order of config.hyp.contains
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
    set(gca,'Ytick', 1 : 5,'Yticklabels',{'STAGE 3','STAGE 2','STAGE 1','REM','WAKE'},'TickDir','out');
    axis tight;
    xl = xlim;
    
    % plot markers
    iplot = 2;
    for markername = unique(mrk_tbl.markerlabel)'
        
        subplot(numel(unique(mrk_tbl.markerlabel))+1,1,iplot); hold;
        
        title(markername);
        sel = find(mrk_tbl.markerlabel == markername);
        for im = sel'
            fill([mrk_tbl.starttime(im), mrk_tbl.endtime(im), mrk_tbl.endtime(im), mrk_tbl.starttime(im)],[0 0 1 1],[0 0 0]);
            
        end
        xlim(xl);
        iplot = iplot + 1;
        
    end
    
    % print to file
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    set(h,'Renderer','Painters');
    print(h, '-dpdf', fullfile(config.imagesavedir,[config.prefix,'part',num2str(ipart),'-hypnogram.pdf']),'-r600');
    
    %     % print a 3 meter version for fun
    %     set(gcf,'PaperUnits','centimeters');
    %     set(gcf,'PaperSize',[300 10]);
    %     h.PaperUnits = 'centimeters';
    %     h.PaperPosition = [0 0 300 10];
    %     h.Units = 'centimeters';
    %     h.PaperSize=[300 10];
    %     h.Units = 'centimeters';
    %     print(h, '-dpdf', fullfile(config.imagesavedir,[config.prefix,'part',num2str(ipart),'-hypnogram_3m.pdf']),'-r600');
    
    writetable(t,fullfile(config.datasavedir,[config.prefix,'part',num2str(ipart),'-hypnogram.txt']));
    
end
disp('Done');
