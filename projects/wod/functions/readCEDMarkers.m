function [CEDStruct]  = readCEDMarkers(cfg, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readCEDMarkers
%
% Read all the event channels of CED data and return a structure similar to
% the one created by readMuseMarkers.m. So it can be used with other
% scripts of EpiCode ((c) Stephen Whitmarsh), or it can be written into 
% Muse marker file.
% 
% ## INPUT :
% cfg.prefix                    : name of the subject
% cfg.CEDrawdir                 : where are the Spike2 data
% cfg.directorylist{part}    : list of all data files, for each part
% force                         : if force == true, force reading again the
%                               events from the data file. If force == 
%                               false, the results are loaded from saved 
%                               file
% 
% ## OUTPUT
% CEDStruct                     : structure with all event timings from
%                               each event or marker channel of the Spike2 
%                               file.
%
% Need of CEDS64ML interface library (loaded with CEDS64LoadLib.m), and of 
% Spike2 software. Can only be ran on Windows (because of the library). To 
% use the CEDStruct in Linux, load it with the 'force = false' argument.
% 
% Paul Baudin
% paul.baudin@live.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = fullfile(cfg.datasavedir,sprintf('%sCEDStruct.mat',cfg.prefix));

if exist(fname,'file') && force == false
    fprintf('*************************************\n');
    fprintf('*** loading precomputed CEDStruct ***\n');
    fprintf('*************************************\n\n');
    load(fname,'CEDStruct');
else
    
    if exist(fname,'file') && force == true
        fprintf('********************************************\n');
        fprintf('*** forced redoing of CEDStruct creation ***\n');
        fprintf('********************************************\n\n');
    else
        fprintf('**************************\n');
        fprintf('*** creating CEDStruct ***\n');
        fprintf('**************************\n\n');
    end
   
    
    % Go through different parts
    fprintf('For %s\n', cfg.prefix(1:end-1));
    
    for ipart = 1 : size(cfg.directorylist,2)
        fprintf('Part %d\n', ipart);
        
        % Go through directory list
        for idir = 1 : size(cfg.directorylist{ipart},2)
            fprintf('Dir %d\n', idir);
            
            clear marks
            
            datapath = fullfile(cfg.CEDrawdir, cfg.directorylist{ipart}{idir});
            fprintf('Extracting marker timings from  %s \n',datapath);
             
            %Open file in Spike2.
            fid = CEDS64Open(datapath);
            if fid<0
                error('error while opening file %s. \nCheck that the path is correct, and that this file is not already opened in Spike2.',datapath);
            end
            
            %get channel nr
            channr = CEDS64MaxChan(fid);
            
            
            %Select all channels which are event or marker channels
            chanindex = [];
            for ichan = 1:channr
                [iType] = CEDS64ChanType(fid, ichan);
                if ismember(iType, [2, 3, 4, 5, 6, 7, 8])
                    chanindex = [chanindex, ichan];
                end
            end
            
            if isempty(chanindex)
                error('no event or marker channels found in %s', datapath);
            else
                fprintf('Found %d event or marker channels in %s\n', length(chanindex), datapath);
            end
            
            
            %Get channel names
            nametemp = []; 
            for imarker = 1:length(chanindex)
                [~, name{imarker}] = CEDS64ChanTitle(fid, chanindex(imarker));
                
                % cant make fieldnames with minusses
                if any(ismember('-',name{imarker}))
                    fprintf('Channel %s is renamed %s (cant make fieldnames with minusses)\n', name{imarker}, strrep(name{imarker},'-','_'));
                    name{imarker} = strrep(name{imarker},'-','_');
                end
                
                % cant make fieldnames with white spaces
                if any(ismember(' ',name{imarker}))
                    fprintf('Channel %s is renamed %s (cant make fieldnames with white spaces)\n', name{imarker}, strrep(name{imarker},' ','_'));
                    name{imarker} = strrep(name{imarker},' ','_');
                end
                
                %create name if is empty
                if isempty(name{imarker})
                    name{imarker} = sprintf('chan%d', chanindex(imarker));
                    fprintf('Channel %d is renamed %s (has no name)\n', chanindex(imarker), name{imarker});
                end
                
                %add 'x' to the field name if it begins with a number
                if ismember(name{imarker}(1), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
                    name{imarker} = insertBefore(name{imarker},name{imarker},'X');
                    fprintf('Channel %s is renamed %s (begins with a number)\n', name{imarker}(2:end), name{imarker});
                end
                
                %rename channel if it has the same name that a previous one
                if any(strcmp(name{imarker}, nametemp))
                    oldname = name{imarker};
                    name{imarker} = sprintf('%s_chan%d', name{imarker}, chanindex(imarker));
                    fprintf('Channel %d %s as the same name that a previous channel : renamed %s\n', chanindex(imarker), oldname, name{imarker});
                    if any(strcmp(oldname, name))
                        [~, idx] = ismember(oldname, name);
                        name{idx} = sprintf('%s_chan%d', name{idx}, chanindex(idx));
                        fprintf('First channel with the name %s is channel %d : renamed %s\n', oldname, chanindex(idx), name{idx});
                    end
                end
                
                %keep info for other channels
                nametemp{imarker} = name{imarker}; 
            end
            
            
            
            %read timings for each channel
            for imarker = 1:length(chanindex)
                [iread, events] = CEDS64ReadEvents(fid,chanindex(imarker),1000000,0);
                
                if iread > 1000000
                    error('Too many events in channel %d of file %s. \nIncrese the corresponding argument in CEDS64ReadEvents just above this line', ichan, datapath);
                end
                if iread < 0
                    error('Error while loading events of channel %d in file %s', chanindex(imarker), datapath);
                end
                
                if ~isempty(events)
                    marks{imarker} = (CEDS64TicksToSecs(fid,events))';
                else
                    marks{imarker} = [];
                end
               
                fprintf('Found %d occurences of %s\n', iread, name{imarker}); 
                
            end
            
            
            
            %recover "real time"
            [~, timeinteger]  = CEDS64TimeDate(fid);
            timems              = timeinteger(1);
            timesec             = timeinteger(2);
            timemin             = timeinteger(3);
            timehour            = timeinteger(4);
            timeday             = timeinteger(5);
            timemonth           = timeinteger(6);
            timeyear            = timeinteger(7);
            timestring          = sprintf('%d/%d/%d %d:%d:%d.%d',timeyear, timemonth, timeday, timehour, timemin, timesec, timems*10);
            starttime           = datetime(timestring, 'Format', 'yyyy/MM/dd HH:mm:ss.SSS');
            maxtime             = CEDS64TicksToSecs(fid, CEDS64MaxTime(fid));
            
            %Add CEDStruct infos. Extra fields are made to be consistant
            %with MuseStruct structure. So it could be used to write Muse
            %markers with writeMusemarkers.m
            CEDStruct{ipart}{idir}.starttime = starttime;
            CEDStruct{ipart}{idir}.directory = cfg.directorylist{ipart}{idir};
            CEDStruct{ipart}{idir}.endtime    = starttime + seconds(maxtime);

            for imarker = 1:length(chanindex)
                CEDStruct{ipart}{idir}.markers.(name{imarker}).events         = [];
                CEDStruct{ipart}{idir}.markers.(name{imarker}).comment        = sprintf('Spike2 chan nr %d', chanindex(imarker));
                CEDStruct{ipart}{idir}.markers.(name{imarker}).color          = '';%info not avalaible in CED librairy
                CEDStruct{ipart}{idir}.markers.(name{imarker}).editable       = 'Yes';
                CEDStruct{ipart}{idir}.markers.(name{imarker}).classid        = sprintf('+%d', imarker);
                CEDStruct{ipart}{idir}.markers.(name{imarker}).classgroupid   = '+3';
                for ievent = 1 : size(marks{imarker},2)
                    CEDStruct{ipart}{idir}.markers.(name{imarker}).trialnum                 = 1;
                    CEDStruct{ipart}{idir}.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent);
                    CEDStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent)            = seconds(marks{imarker}(ievent)) + CEDStruct{ipart}{idir}.starttime;
                end
            end
        end%idir
    end%ipart
    
    save(fname,'CEDStruct');
    CEDS64CloseAll();
end %exist && force

end