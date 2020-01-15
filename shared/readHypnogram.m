function [marker, hyp]  = readHypnogram_parts(cfg, MuseStruct_micro, MuseStruct_macro)


% elect those markers to load
markerlist = [];
for i = 1 : size(cfg.name,2)
    if ismember(cfg.name{i},cfg.name)
        markerlist = [markerlist, i];
    end
end

marker = table;
istage = 0;

for imarker = markerlist
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            try
                StartRecord(ipart,idir) = MuseStruct_macro{ipart}{idir}.markers.StartRecord.clock;
                StopRecord(ipart,idir)  = MuseStruct_macro{ipart}{idir}.markers.StopRecord.clock;
            catch
            end
            
            if isfield(MuseStruct_macro{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                if isfield(MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'offset')
                    for ievent = 1 : size(MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).offset,2)
                        
                        istage = istage + 1;
                        
                        marker.stage(istage)  = -2;
                        marker.clock(istage)  = MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(ievent);
                        marker.name{istage}   = cfg.muse.startend{imarker,1};
                        marker.ipart(istage)  = ipart;
                        marker.idir(istage)   = idir;
                        
                        % find overlap with hypnogram markers
                        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                            if isfield(MuseStruct_macro{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                                for i = 1 : size(MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).offset,2)
                                    x1 = MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).offset(ievent);
                                    x2 = MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).offset(ievent);
                                    y1 = MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).offset(i);
                                    y2 = MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).offset(i);
                                    
                                    %                                             if intersect(x1:x2,y1:y2)
                                    if (y1 < x1) && (x1 < y2)
                                        
                                        fprintf('Found "%s" overlapping with "%s" : adding to trialinfo: ',cfg.name{imarker},cell2mat(hyplabel));
                                        switch cell2mat(hyplabel)
                                            case 'PHASE_1'
                                                fprintf('%d\n',1);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 1;
                                                MuseStruct_micro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 1;
                                                marker.stage(istage)  = 1;
                                            case 'PHASE_2'
                                                fprintf('%d\n',2);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 2;
                                                MuseStruct_micro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 2;
                                                marker.stage(istage)  = 2;
                                            case 'PHASE_3'
                                                fprintf('%d\n',3);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 3;
                                                MuseStruct_micro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 3;
                                                marker.stage(istage)  = 3;
                                            case 'REM'
                                                fprintf('%d\n',4);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 4;
                                                MuseStruct_micro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 4;
                                                marker.stage(istage)  = 4;
                                            case 'AWAKE'
                                                fprintf('%d\n',0);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 0;
                                                MuseStruct_micro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = 0;
                                                marker.stage(istage)  = 0;
                                            case 'NO_SCORE'
                                                fprintf('%d\n',0);
                                                MuseStruct_macro{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).stage = -1;
                                                marker.stage(istage)  = -1;
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
        end
    end
end



hyp = table;
ihyp = 0;


% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)
    
    % Go through directory list
    for idir = 1 : size(cfg.directorylist{ipart},2)
        
        % find overlap with hypnogram markers
        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
            if isfield(MuseStruct_macro{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                for i = 1 : size(MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).offset,2)
                    ihyp = ihyp + 1;
                    hyp.starttime(ihyp)   = MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).clock(i);
                    hyp.endtime(ihyp)     = MuseStruct_macro{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).clock(i);
                    hyp.duration(ihyp)    = seconds(hyp.endtime(ihyp) - hyp.starttime(ihyp));
                    switch cell2mat(hyplabel)
                        case 'PHASE_1'
                            hyp.stage(ihyp)  = 1;
                        case 'PHASE_2'
                            hyp.stage(ihyp)  = 2;
                        case 'PHASE_3'
                            hyp.stage(ihyp)  = 3;
                        case 'REM'
                            hyp.stage(ihyp)  = 4;
                        case 'AWAKE'
                            hyp.stage(ihyp)  = 0;
                        case 'NO_SCORE'
                            hyp.stage(ihyp)  = -1;
                        otherwise
                            error('Unexpected label name in Hypnogram\n');
                    end
                end
            end
        end
        
    end
end

hyp.stage(hyp.stage == -1) = 0;

hyp = sortrows(hyp);

for imarker = markerlist
    
    for i = 0 : 4
        totaldur(i+1) = sum(hyp.duration(hyp.stage == i));
    end
    for i = 0 : 4
        totalsum(i+1) = sum(marker.stage == i & strcmp(marker.name, cfg.name{imarker}));
    end
    
    h= figure;
    subplot(4,1,1);
    bar([0 : 4],totalsum);
    text([0 : 4],totalsum,num2str(totalsum'),'vert','bottom','horiz','center');
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('Number of IEDs detected in each stage');
    box off
    ylim([0, max(totalsum)*1.3]);
    set(gca,'TickLength',[0 0])
    
    t = arrayfun(@(x) sprintf('%3.2f',x),totaldur/60/60,'uniformoutput', false);
    
    subplot(4,1,2);
    bar([0 : 4],totaldur/60/60);
    text([0 : 4],totaldur/60/60,t,'vert','bottom','horiz','center');
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('Total duration per stage (hrs)');
    box off
    ylim([0, max(totaldur/60/60)*1.3]);
    set(gca,'TickLength',[0 0])
    
    t = arrayfun(@(x) sprintf('%3.2f',x),totalsum./totaldur*60,'uniformoutput', false);
    
    subplot(4,1,3);
    bar([0 : 4],totalsum./totaldur*60);
    text([0 : 4],totalsum./totaldur*60,t,'vert','bottom','horiz','center');
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('IEDs per minute');
    box off
    ylim([0, max(totalsum./totaldur*60)*1.3]);
    set(gca,'TickLength',[0 0])
    
    y = totalsum./totaldur*60;
    y = y ./ y(1);
    t = arrayfun(@(x) sprintf('%3.2f',x),y,'uniformoutput', false);
    
    subplot(4,1,4);
    bar([0 : 4],y);
    text([0 : 4],y,t,'vert','bottom','horiz','center');
    xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
    title('IEDs per minute normalized to wake stage');
    box off
    ylim([0, max(y)*1.3]);
    set(gca,'TickLength',[0 0])
    
    % print to file
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    set(h,'Renderer','Painters');
    print(h, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{imarker},'_x_sleep_stage.pdf']),'-r600');
    % writetable(t,fullfile(config.imagesavedir,'seg_labels'));
    
    
    disp('Done');
end