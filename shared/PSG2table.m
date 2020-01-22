function [hypnogram] = PSG2table(cfg, MuseStruct, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [MuseStruct, marker, hypnogram] = hypnogramStats(cfg, MuseStruct, force)
%
% Outputs table of hypnogram
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all');

fname = fullfile(cfg.datasavedir,sprintf('%sPSG2table.mat',cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('**************************\n');
    fprintf('** loading Hypnogram *****\n');
    fprintf('**************************\n\n');
    load(fname,'hypnogram');
else
    
    if force == true
        fprintf('************************************\n');
        fprintf('** forced redoing of Hypnogram *****\n');
        fprintf('************************************\n\n');
    else
        fprintf('***************************\n');
        fprintf('** creating Hypnogram *****\n');
        fprintf('***************************\n\n');
    end
  
    % eventnr will increased over all marker events and all files
    eventnr = 0;
    
    hypnogram = table;
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            
            % Find hypnogram markers
            for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                
                if isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                    for i = 1 : length(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).clock)
                        eventnr = eventnr + 1;
                        hypnogram.night(eventnr)            = ipart;                        
                        hypnogram.stage(eventnr)            = hyplabel;
                        hypnogram.directory(eventnr)        = {MuseStruct{ipart}{idir}.directory};                        
                        hypnogram.startfilesec(eventnr)     = round(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i));
                        hypnogram.endfilesec(eventnr)       = round(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i));
                        hypnogram.duration(eventnr)         = hypnogram.endfilesec(eventnr) - hypnogram.startfilesec(eventnr);                             
                        hypnogram.startdatetime(eventnr)    = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).clock(i);
                        hypnogram.enddatetime(eventnr)      = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).clock(i);                  
                    end
                end
                
            end
        end
    end
    
    hypnogram = sortrows(hypnogram,'startdatetime');
    
    writetable(hypnogram,fullfile(cfg.datasavedir,[cfg.prefix,'hypnogram.txt']));
    writetable(hypnogram,fullfile(cfg.datasavedir,[cfg.prefix,'hypnogram.xls']));

    % save results
    save(fname,'hypnogram');
end
