function [CEDStruct]  = readCEDmarkers(cfg, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readCEDevents
%
% Read all the event channels of CED data and return a structure similar to
% the one created by readMuseMarkers.m. So it can be used with other
% scripts of EpiCode ((c) Stephen Whitmarsh), or it can be written into 
% Muse marker file.
% 
% ## INPUT :
% cfg.prefix                    : name of the subject
% cfg.rawdir                 : where are the Spike2 data
% cfg.directorylist{part}    : list of all data files, for each part
% cfg.datasavedir            : where to save output file
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
            
            temp = dir(fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir}, '.smr*']));%.smr* because some data are .smr and other .smrx
            datapath = fullfile(cfg.rawdir, temp.name);
            fprintf('Extracting marker timings from  %s \n',datapath);
             
            %Open file in Spike2.
            fid = CEDS64Open(datapath);
            if fid<0
                error('error while opening file %s. \nCheck that the path is correct, and that this file is not already opened in Spike2.',datapath);
            end
            
            %get channel nr
            channr = CEDS64MaxChan(fid);
            
            
            %Select all channels which are event or marker channels
            chanindex       = [];
            textmarkindex   = [];
            for ichan = 1:channr
                [iType] = CEDS64ChanType(fid, ichan);
                if ismember(iType, [2, 3, 4, 5]) %FIXME : 6 and 7 : realmark, wavemark. See if need to implement for those chan types
                    chanindex = [chanindex, ichan];
                elseif iType == 8 %textmark
                    chanindex     = [chanindex, ichan];
                    textmarkindex = [textmarkindex, ichan];
                end
            end
            
            if isempty(chanindex)
                error('no event or marker channels found in %s', datapath);
            else
                fprintf('Found %d event or marker channels in %s\n', length(chanindex), datapath);
            end
            
            
            %Get channel names and rename it if needed
            nametemp = []; 
            for imarker = 1:length(chanindex)
                [~, name{imarker}] = CEDS64ChanTitle(fid, chanindex(imarker));
                name{imarker}      = renamechan_CED(name{imarker}, chanindex(imarker),name(1:imarker-1),true); %name{1:imarker-1} is [] if imarker = 1
            end

            
            %read timings for each channel
            for imarker = 1:length(chanindex)
                [iread, events] = CEDS64ReadEvents(fid,chanindex(imarker),1000000,0);
                
                if iread > 1000000
                    error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadEvents just above this line', ichan, datapath);
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
            CEDStruct{ipart}{idir}.endtime   = starttime + seconds(maxtime);

            for imarker = 1:length(chanindex)
                CEDStruct{ipart}{idir}.markers.(name{imarker}).comment        = sprintf('Spike2 chan nr %d', chanindex(imarker));
                if ismember(chanindex(imarker), textmarkindex)
                    [iread, marksinfos] = CEDS64ReadExtMarks(fid, chanindex(imarker),1000000,0);
                    if iread > 1000000
                        error('Too many events in channel %d of file %s. \nIncrease the corresponding argument in CEDS64ReadEvents just above this line', ichan, datapath);
                    end
                    %get code from marks info
                    for ievent = 1:size(marksinfos, 1)
%                         CEDStruct{ipart}{idir}.markers.(name{imarker}).textmark(ievent) = marksinfos(ievent).m_Data;
                        CEDStruct{ipart}{idir}.markers.(name{imarker}).code_1(ievent)   = marksinfos(ievent).m_Code1;
                        CEDStruct{ipart}{idir}.markers.(name{imarker}).code_2(ievent)   = marksinfos(ievent).m_Code2;
                        CEDStruct{ipart}{idir}.markers.(name{imarker}).code_3(ievent)   = marksinfos(ievent).m_Code3;
                        CEDStruct{ipart}{idir}.markers.(name{imarker}).code_4(ievent)   = marksinfos(ievent).m_Code4;
                    end
                end
                for ievent = 1 : size(marks{imarker},2)
                    CEDStruct{ipart}{idir}.markers.(name{imarker}).synctime(ievent)         = marks{imarker}(ievent);
                    CEDStruct{ipart}{idir}.markers.(name{imarker}).clock(ievent)            = seconds(marks{imarker}(ievent)) + CEDStruct{ipart}{idir}.starttime;
                end
            end
        end%idir
    end%ipart
    
    save(fname,'CEDStruct');
end %exist && force

end