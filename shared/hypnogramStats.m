function [MuseStruct, marker, hypnogram] = hypnogramStats(cfg, MuseStruct, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [MuseStruct, marker, hypnogram] = hypnogramStats(cfg, MuseStruct, force)
%
% Outputs and plots occurance of markers in different sleep stages
% 
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','all');

fname = fullfile(cfg.datasavedir,sprintf('%shypnogramStats.mat',cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('***************************************\n');
    fprintf('** loading precomputed MuseStruct *****\n');
    fprintf('***************************************\n\n');
    load(fname,'MuseStruct', 'marker', 'hypnogram');
else
    
    if force == true
        fprintf('**********************************************\n');
        fprintf('** forced redoing of MuseStruct creation *****\n');
        fprintf('**********************************************\n\n');
    else
        fprintf('****************************\n');
        fprintf('** creating MuseStruct *****\n');
        fprintf('****************************\n\n');
    end
    
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
                
                if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                    if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                        for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)
                            
                            eventnr = eventnr + 1;
                            
                            marker.stage(eventnr)  = -2;
                            marker.clock(eventnr)  = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(ievent);
                            marker.name{eventnr}   = cfg.muse.startend{imarker,1};
                            marker.ipart(eventnr)  = ipart;
                            marker.idir(eventnr)   = idir;
                            
                            % find overlap with hypnogram markers
                            for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                                if isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
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
    
    hypnogram = sortrows(hypnogram);
    
    for imarker = markerlist
        
        for i = 0 : 4
            totaldur(i+1) = seconds(sum(hypnogram.duration(hypnogram.stage == i)));
        end
        for i = 0 : 4
            totalsum(i+1) = sum(marker.stage == i & strcmp(marker.name, cfg.name{imarker}));
        end
        
        h = figure;
        subplot(4,1,1);
        bar([0 : 4],totalsum,'k');
        text([0 : 4],totalsum,num2str(totalsum'),'vert','bottom','horiz','center','FontSize',18);
        xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
        title('Number of IEDs detected in each stage');
        box off
        ylim([0, max(totalsum)*1.3]);
        set(gca,'TickLength',[0 0])
        set(gca,'Yticklabels',[])
        
        t = arrayfun(@(x) sprintf('%3.2f',x),totaldur/60/60,'uniformoutput', false);
        
        subplot(4,1,2);
        bar([0 : 4],totaldur/60/60,'k');
        text([0 : 4],totaldur/60/60,t,'vert','bottom','horiz','center','FontSize',18);
        xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
        title('Total duration per stage (hrs)');
        box off
        ylim([0, max(totaldur/60/60)*1.3]);
        set(gca,'TickLength',[0 0])
        set(gca,'Yticklabels',[])
        
        t = arrayfun(@(x) sprintf('%3.2f',x),totalsum./totaldur*60,'uniformoutput', false);
        
        subplot(4,1,3);
        bar([0 : 4],totalsum./totaldur*60,'k');
        text([0 : 4],totalsum./totaldur*60,t,'vert','bottom','horiz','center','FontSize',18);
        xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
        title('IEDs per minute');
        box off
        ylim([0, max(totalsum./totaldur*60)*1.3]);
        set(gca,'TickLength',[0 0])
        set(gca,'Yticklabels',[])
        
        y = totalsum./totaldur*60;
        y = y ./ y(1);
        t = arrayfun(@(x) sprintf('%3.2f',x),y,'uniformoutput', false);
        
        subplot(4,1,4);
        bar([0 : 4],y,'k');
        text([0 : 4],y,t,'vert','bottom','horiz','center','FontSize',18);
        xticklabels({'Wake','Stage 1','Stage 2','Stage 3','REM'})
        title('IEDs per minute normalized to wake stage');
        box off
        ylim([0, max(y)*1.3]);
        set(gca,'TickLength',[0 0])
        set(gca,'Yticklabels',[])
        
        % print to file
        set(h,'PaperOrientation','landscape');
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition', [0 0 1 1]);
        set(h,'Renderer','Painters');
        print(h, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,cfg.name{imarker},'_x_sleep_stage.pdf']),'-r600');
    end
    
    % save results
    save(fname,'MuseStruct', 'marker', 'hypnogram');
end
