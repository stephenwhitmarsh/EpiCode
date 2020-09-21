function [SpikeTrials] = readSpykingCircus_spikestage(cfg, SpikeRaw, MuseStruct, part_list, spikestage_list, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% readSpykingCircus_spikestage
%
% Define spike trials according to MuseStruct and cfg.spikestage
%
% Exemple of config needed :
%
% cfg.spikestage{1}.markerstart          = 'Crise_End';
% cfg.spikestage{1}.markerstartoffset    = 1; %seconds
% cfg.spikestage{1}.markerend            = 'SlowWave';
% cfg.spikestage{1}.indexEnd             = 1; %if 1 : take the i start
% marker and the i+1 end marker
% cfg.spikestage{1}.periodmin            = 30; %seconds
% cfg.spikestage{1}.triallength          = 5; %seconds
% cfg.spikestage{1}.subdivision          = 4; %minimum 2
% cfg.spikestage{1}.includeend           = true; %add the end of the period
% to the previous subdivisions
% cfg.spikestage{1}.markerendoffset      = -1; %seconds
% cfg.spikestage{1}.removeartefact       = 'period'; %period, trial, none
% cfg.spikestage{1}.stagenumber          = [1, 2, 3, 4, 5]; %array of
% integers. if single value : gives the same number to all the trials
%
% cfg.spikestage{2}.markerstart          = 'Baseline_Start';
% cfg.spikestage{2}.markerstartoffset    = 0; %seconds
% cfg.spikestage{2}.markerend            = 'Injection';
% cfg.spikestage{2}.indexEnd             = 0;
% cfg.spikestage{2}.periodmin            = 600; %seconds
% cfg.spikestage{2}.triallength          = 60; %seconds
% cfg.spikestage{2}.subdivision          = 'all'; %minimum 2
% cfg.spikestage{2}.includeend           = false;
% cfg.spikestage{2}.markerendoffset      = 0; %seconds
% cfg.spikestage{2}.removeartefact       = 'trial';
% cfg.spikestage{2}.stagenumber          = 0; %array of integers. to give
% the same number to all the trials, write one single value
%
% Paul Baudin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fname = fullfile(cfg.datasavedir,[cfg.prefix,'spiketrials_stages', cfg.circus.postfix,'.mat']);


if exist(fname,'file') && force == false
    load(fname,'SpikeTrials');
else
    
    for ipart = part_list
                
        temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
        hdr_fname   = fullfile(temp(1).folder,temp(1).name);
        hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        
        if strcmp(spikestage_list, 'all')
            spikestage_list = 1:length(cfg.spikestage);
        end
        
        itrial = 0;
        
        for ispikestage = spikestage_list
            
            clear markerstart_synctime markerend_synctime BAD_Start_synctime BAD_End_synctime
            %clear starttime startsample endtime endsample offset trialnr dirnr stage
            
            if ~(cfg.spikestage{ispikestage}.indexEnd == 1 || cfg.spikestage{ispikestage}.indexEnd == 0)
                error('indexEnd can only be 1 or 0.');
            end
            
            if cfg.spikestage{ispikestage}.subdivision < 2
                error('Need of minimum 2 subdivisions, set cfg.spikestage{ispikestage}.subdivision with an integer >= 2');
            end
            
            if strcmp(cfg.spikestage{ispikestage}.subdivision, 'all')
                stagenumber_forall = cfg.spikestage{ispikestage}.stagenumber(1);
            end
            
            %get all marker timings
            [~,markerstart_synctime]        = concatenateMuseMarkers(MuseStruct, ipart, cfg.spikestage{ispikestage}.markerstart);
            [~,markerend_synctime]          = concatenateMuseMarkers(MuseStruct, ipart, cfg.spikestage{ispikestage}.markerend);
            [~,BAD_Start_synctime]          = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__START__');
            [~,BAD_End_synctime]            = concatenateMuseMarkers(MuseStruct, ipart, 'BAD__END__');
            
            %safety check
            if ~(size(markerstart_synctime,2) == size(markerend_synctime,2))
                error('Not same amount of %s and %s. \nFor %s part %d',cfg.spikestage{ispikestage}.markerstart, cfg.spikestage{ispikestage}.markerend, cfg.prefix(1:end-1), ipart);
            end
            if size(BAD_Start_synctime,2) ~= size(BAD_End_synctime,2)
                error('Not same amount of BAD Start and End markers. \nFor %s part %d', cfg.prefix(1:end-1), ipart);
            end
            
            for iperiod = 1:size(markerstart_synctime,2) - cfg.spikestage{ispikestage}.indexEnd
                
                periodlength = markerend_synctime{1,iperiod + cfg.spikestage{ispikestage}.indexEnd } - markerstart_synctime{1,iperiod};
                
                %ignore too small periods
                if periodlength > cfg.spikestage{ispikestage}.periodmin
                    
                    %find if iperiod has artefact
                    hasartefact = false;
                    
                    if strcmp(cfg.spikestage{ispikestage}.removeartefact, 'period')
                        
                        periodinterval = round(markerstart_synctime{1,iperiod} * hdr.Fs : markerend_synctime{1,iperiod + cfg.spikestage{ispikestage}.indexEnd} * hdr.Fs);
                        
                        for iartefact = 1 : size(BAD_Start_synctime,2)
                            artefactinterval = round(BAD_Start_synctime{1, iartefact} * hdr.Fs : BAD_End_synctime{1, iartefact} *hdr.Fs);
                            if intersect(periodinterval,artefactinterval) %intervals are expressed in sample points
                                fprintf('Found artefact in part %d, spikestage %d, period %d : period not loaded (%s)\n',ipart, ispikestage ,iperiod, cfg.prefix(1:end-1));
                                hasartefact= true;
                                break;
                            end
                        end
                    end
                    
                    %ignore trial if has artefact
                    if ~hasartefact
                        
                        %find how many cuts to do
                        if strcmp(cfg.spikestage{ispikestage}.subdivision, 'all')
                            n_cuts = fix(periodlength/cfg.spikestage{ispikestage}.triallength); %inferior integer
                        else
                            n_cuts = cfg.spikestage{ispikestage}.subdivision;
                        end
                        
                        %check if not too many trials compared to period
                        total_cut_length =  n_cuts * cfg.spikestage{ispikestage}.triallength + cfg.spikestage{ispikestage}.markerstartoffset + cfg.spikestage{ispikestage}.markerendoffset;
                        if total_cut_length > periodlength
                            error('Sum of all cutted trials is bigger than the period, part %d, spikestage %d', ipart, ispikestage);
                        end
                        
                        %if asked in cfg : number of the stage is the same for all the cuts.
                        if exist('stagenumber_forall')
                            cfg.spikestage{ispikestage}.stagenumber = [];
                            cfg.spikestage{ispikestage}.stagenumber = ones(1,n_cuts) * stagenumber_forall;
                        end
                        
                        %error if not the same amout of stage number and of cuts :
                        if length(cfg.spikestage{ispikestage}.stagenumber) ~= n_cuts + cfg.spikestage{ispikestage}.includeend
                            error('Not the same amount of cuts and of stage numbers, part %d, spikestage %d, period %d', ipart, ispikestage, iperiod);
                        end
                        
                        
                        %cut the period in different trials
                        for i_cutperiod = 0 : n_cuts-1 %start at zero
                            
                            starttime_temp = markerstart_synctime{1,iperiod} + cfg.spikestage{ispikestage}.markerstartoffset + i_cutperiod*periodlength/n_cuts;
                            endtime_temp = markerstart_synctime{1,iperiod} + cfg.spikestage{ispikestage}.markerstartoffset + i_cutperiod*periodlength/n_cuts + cfg.spikestage{ispikestage}.triallength;
                            
                            %find if i_cutperiod has artefact
                            hasartefact = false;
                            
                            if strcmp(cfg.spikestage{ispikestage}.removeartefact, 'trial')
                                
                                trialinterval = round(starttime_temp * hdr.Fs : endtime_temp * hdr.Fs);
                                
                                for iartefact = 1 : size(BAD_Start_synctime,2)
                                    artefactinterval = round(BAD_Start_synctime{1, iartefact} * hdr.Fs : BAD_End_synctime{1, iartefact} *hdr.Fs);
                                    if intersect(trialinterval,artefactinterval) %intervals are expressed in sample points
                                        fprintf('Found artefact in part %d, spikestage %d, period %d, trial %d : trial not loaded (%s)\n',ipart, ispikestage ,iperiod, i_cutperiod, cfg.prefix(1:end-1));
                                        hasartefact= true;
                                        break;
                                    end
                                end
                            end
                            
                            %ignore trial if has artefact
                            if ~hasartefact
                                
                                itrial = itrial + 1;
                                
                                starttime(itrial)   = starttime_temp;
                                startsample(itrial) = starttime(itrial)* hdr.Fs;
                                endtime(itrial)     = endtime_temp;
                                endsample(itrial)   = endtime(itrial) * hdr.Fs;
                                offset(itrial)      = 0;
                                trialnr(itrial)     = itrial;
                                dirnr(itrial)       = markerstart_synctime{2,iperiod};
                                stage(itrial)       = cfg.spikestage{ispikestage}.stagenumber(i_cutperiod+1); %+1 because i_cutperiod starts at zero
                            
                            end
                            
                        end %i_cutperiod
                        
                        %last trial
                        if cfg.spikestage{ispikestage}.includeend
                            
                            itrial = itrial + 1;
                            
                            starttime(itrial)   = markerend_synctime{1,iperiod+cfg.spikestage{ispikestage}.indexEnd} + cfg.spikestage{ispikestage}.markerendoffset - cfg.spikestage{ispikestage}.triallength;
                            startsample(itrial) = starttime(itrial) * hdr.Fs;
                            endtime(itrial)     = markerend_synctime{1,iperiod+cfg.spikestage{ispikestage}.indexEnd} + cfg.spikestage{ispikestage}.markerendoffset;
                            endsample(itrial)   = endtime(itrial) * hdr.Fs;
                            offset(itrial)      = 0;
                            trialnr(itrial)     = itrial;
                            dirnr(itrial)       = markerstart_synctime{2,iperiod};
                            stage(itrial)       = cfg.spikestage{ispikestage}.stagenumber(end);
                            
                        end
                        
                        
                        
                    end
                end
            end %iperiod
            
        end %ispikestage
        
        % create Fieldtrip trl
        cfgtemp                             = [];
        cfgtemp.trl                         = round([startsample; endsample; offset]');
        cfgtemp.trl(:,4)                    = trialnr;
        cfgtemp.trl(:,6)                    = dirnr; %not column 5, to be consistent with Stephen's LFP trialinfo numbers
        cfgtemp.trl(:,7)                    = stage;
        cfgtemp.trl(:,8)                    = starttime;
        cfgtemp.trl(:,9)                    = endtime;
        cfgtemp.trl                         = cfgtemp.trl(startsample > 0 & endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs
        cfgtemp.trlunit                     = 'samples';
        cfgtemp.hdr                         = hdr;
        SpikeTrials{ipart}       = ft_spike_maketrials(cfgtemp,SpikeRaw{ipart});
        
        
        % to be consistent with Stephen's hypnogram trialinfo :
        SpikeTrials{ipart}.trialinfo             = table;
        SpikeTrials{ipart}.trialinfo.stage       = stage';
        SpikeTrials{ipart}.trialinfo.starttime   = starttime';
        SpikeTrials{ipart}.trialinfo.endtime     = endtime';
        SpikeTrials{ipart}.trialinfo.trialnr     = trialnr';
        
    end  %ipart
    
    save(fname,'SpikeTrials');
    
end %if exist fname && force == false



end




