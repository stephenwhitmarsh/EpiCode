function [data, MuseStruct] = removetrialLFP_marker(cfg, data, MuseStruct, partlist, markerlist, method, markerstart, markerend, indexstart, indexend, timefrombegin, timefromend, plotdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function removetrialLFP_marker
% 
% ### INPUT :
% cfg                : all the parameters specific for the subject, defined
%                      in a separated script.
% data               : LFP data obtained by the readLFP.m function. data is
%                      of the type data{ipart}{imarker}.
% MuseStruct         : structure obtained by the readMuseMarkers.m script.
%                      It contains all the time information of the marker 
%                      position in the data.
% partlist           : list of parts of data to process. Can be set to 
%                      'all'.
% markerlist         : list of data group of trials to process. ('imarker'  
%                      in data{ipart}{imarker})  Can be set to 'all'.
% method : 
% 'remove'           : remove all trials which intersect with the period
% 'keep'             : keep only trials which intersect with the period
% markerstart/end    : name of the marker which define the begin and the
%                      end of the period.
% indexstart/end     : if need to take the +1 marker, i.e. for a period
%                      defined between samemarker(i) and samemarker(i+1).
%                      Otherwise let it as 0.
% timefrombegin/end  : relative time from markerstart and markerend, in 
%                      seconds, to define the period.
% plotdata           : true or false. Either to plot the results or not.
% 
% ### OUTPUT :
% data               : returned without the selected trials. The field 
%                      modified are : 
%                      - data{ipart}{imarker}.time
%                      - data{ipart}{imarker}.trial
%                      - data{ipart}{imarker}.trialinfo
% MuseStruct         : no modification of the marker timings. Add a field 
%                      with the indexes of the removed trials in each
%                      MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker},1)
% 
% Paul Baudin, with the help of Stephen Whitmarsh
% paul.baudin@live.fr
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*****************************************************\n');
fprintf('Remove trials in %s :\n',cfg.prefix(1:end-1));
fprintf('- %s all trials which intersect the defined period.\n',method);
fprintf('- period begins with each marker %s(i+%d)%+gs\n', markerstart,indexstart, timefrombegin);
fprintf('- period ends with each marker %s(i+%d)%+gs\n', markerend,indexend, timefromend);
fprintf('*****************************************************\n\n');

initial_prefix = cfg.prefix; %for merge

if strcmp(partlist, 'all')
    partlist = 1:size(data, 2);
end

for ipart = partlist
    
    fprintf('For part %d\n', ipart);
    
    if strcmp(markerlist, 'all')
        markerlist = 1:size(data{ipart}, 2);
    end
    
    for imarker = markerlist
        
        fprintf('For LFP name %s\n', cfg.LFP.name{imarker});
        
        
        if ~isempty(data{ipart}{imarker})
            
            n_trials        = 0;
            trialremoved    = false(1,size(data{ipart}{imarker}.trial,2));
            
            % for ipart = 1:length(MuseStruct)
            for idir = 1:length(MuseStruct{ipart})
                nr_notloaded(idir) = 0;
                
                if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                    if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}), 'synctime')
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).trialremoved = false(1,length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime));
                        
                        % if trial begin before first sample or end after last
                        % sample : is not loaded with readLFP
                        
                        for itrialMuse = 1:length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                            %if trial of the dir (data.trialinfo(:,3)==idir) does note have
                            %the same index (data.trialinfo(:,1) as Muse LFP marker
                            %then annotate as trialremoved
                            if ~ismember(itrialMuse, data{ipart}{imarker}.trialinfo((data{ipart}{imarker}.trialinfo(:,3)==idir),1)) 
                                nr_notloaded(idir) = nr_notloaded(idir)+1;
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).trialremoved(itrialMuse) = true;
                                fprintf('%s %s part %d dir %d : trial %d not loaded\n',cfg.prefix(1:end-1),cfg.muse.startend{imarker,1}, ipart, idir, itrialMuse);
                            end
                        end
                        
                        % check if there is an equal amount of start and end markers
                        if isfield(MuseStruct{ipart}{idir}.markers,markerstart)
                            if isfield(MuseStruct{ipart}{idir}.markers.(markerstart), 'synctime')
                                if size(MuseStruct{ipart}{idir}.markers.(markerstart).synctime,2)-size(MuseStruct{ipart}{idir}.markers.(markerend).synctime,2) == 0
                                    
                                    
                                    %convert idir{itrial} to data{ipart}{imarker}{itrial} :
                                    for itrial = n_trials+1 : n_trials+sum(data{ipart}{imarker}.trialinfo(:,3)==idir)
                                        itrialMuse = data{ipart}{imarker}.trialinfo(itrial,1); %trial of loaded data
                                        
                                        %trial as defined in readLFP.m
                                        trialinterval_begin     = (MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrialMuse) + cfg.epoch.toi{imarker}(1) - cfg.epoch.pad{imarker})*data{ipart}{imarker}.fsample;
                                        trialinterval_end       = (MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(itrialMuse) + cfg.epoch.toi{imarker}(2) - cfg.epoch.pad{imarker})*data{ipart}{imarker}.fsample;
                                        trialinterval           = round(trialinterval_begin : trialinterval_end);
                                        
                                        
                                        
                                        for isearchtrial = 1 : size(MuseStruct{ipart}{idir}.markers.(markerstart).synctime,2)
                                            %check that we dont exceed the
                                            %index values
                                            if isearchtrial + indexstart <= length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                                                if isearchtrial + indexend <= length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime)
                                                    
                                                    searchtrialinterval_begin       = (MuseStruct{ipart}{idir}.markers.(markerstart).synctime(isearchtrial + indexstart) + timefrombegin) * data{ipart}{imarker}.fsample;
                                                    searchtrialinterval_end         = (MuseStruct{ipart}{idir}.markers.(markerend).synctime(isearchtrial + indexend) + timefromend) * data{ipart}{imarker}.fsample;
                                                    
                                                else %first marker exists but last marker exceeds filed size
                                                    %end is the end of the file.
                                                    searchtrialinterval_end = seconds(MuseStruct{ipart}{idir}.starttime - MuseStruct{ipart}{idir}.endtime) * data{ipart}{imarker}.fsample;
                                                end
                                                
                                                searchtrialinterval             = round(searchtrialinterval_begin : searchtrialinterval_end);
                                                
                                                if intersect(trialinterval,searchtrialinterval) %intervals are expressed in sample points
                                                    %fprintf('Found searched-trial in part %d trial %d (MuseStruct dir %d, marker %s n°%d)\n',ipart ,itrial,idir, cfg.muse.startend{imarker,1},itrialMuse);
                                                    trialremoved(itrial) = 1;
                                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).trialremoved(itrialMuse) = true;
                                                end
                                            end
                                        end
                                        
                                    end %itrial
                                else
                                    error('Not the same amout of markerstart (%s) and markerend (%s) found in %s \n',markerstart,markerend,MuseStruct{ipart}{idir}.directory);
                                    
                                end
                                
                            end
                        end
                        if strcmp(method,'keep')
                            MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).trialremoved = ~MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).trialremoved;
                        end
                    end
                end
                
                n_trials = n_trials + sum(data{ipart}{imarker}.trialinfo(:,3)==idir); %nr of trials per dir
            end %idir
            %end
            
            if strcmp(method, 'keep')
                trialremoved = ~trialremoved;
            end
            
            if sum(trialremoved) > 0
                fprintf('Remove %d searched-trials\n', sum(trialremoved));
            else
                fprintf('No searched-trials found : all trials are kept\n');
            end
            
                
                
            if plotdata
                %h automatic setting :
                for itrial = 1 : length(data{ipart}{imarker}.trial)
                    t_0 = -(cfg.epoch.toi{imarker}(1)-cfg.epoch.pad{imarker}(1))*data{ipart}{imarker}.fsample; % offset for which t = 0;
                    h_temp(itrial) = max(data{ipart}{imarker}.trial{itrial}(1,round(-0.5*data{ipart}{imarker}.fsample)+t_0: round(0.5*data{ipart}{imarker}.fsample)+t_0)); %max between -0.5s and 0.5s. Avoid noise. Available for EEG and EMG.
                end
                
                h = abs(mean(h_temp)*2);
                %         if isMicromed
                %             h = mean(h_temp)/2;
                %         elseif isBrainvision
                %             h = mean(h_temp)*2;
                %         end
                
                cfgtemp = [];
                cfgtemp.channel = cfg.LFP.electrodetoplot{imarker};
                data_plot = ft_selectdata(cfgtemp, data{ipart}{imarker});
                
                fig = figure;
                hold;
                for itrial = 1 : size(data{ipart}{imarker}.trial,2)
                    if ~trialremoved(itrial)
                        color = 'k';
                    else
                        color = 'r';
                    end
                    plot(data_plot.time{itrial}, data_plot.trial{itrial}+(n_trials+1)*h-itrial*h, 'color', color);
                end
                
                xlabel(sprintf('Time from %s (s)', cfg.muse.startend{imarker,1}),'Interpreter','none', 'Fontsize',15);
                ylabel('Number of trials', 'Fontsize',15);
                if sum(nr_notloaded) == 0
                    title(sprintf('%s chan %s \nAll trials are loaded \n%d searched-trials removed from %d trials \nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{imarker}, data_plot.label{1}, sum(trialremoved), length(data{ipart}{imarker}.trial), markerstart,indexstart, timefrombegin, markerend,indexend, timefromend),'Fontsize',20,'Interpreter','none');
                else
                    title(sprintf('%s chan %s \n%d trials are not loaded (too close from begin or end of data) \n%d searched-trials trials removed from %d trials\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{imarker}, data_plot.label{1}, sum(nr_notloaded), sum(trialremoved), length(data{ipart}{imarker}.trial), markerstart,indexstart, timefrombegin, markerend,indexend, timefromend),'Fontsize',20,'Interpreter','none');
                end
                set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
                tick_interval = round(n_trials/10);
                if tick_interval == 0
                    tick_interval = 1;
                end
                yticks(h : h*tick_interval : n_trials*h);
                yticklabels(n_trials : -tick_interval : 0);
                set(gca,'TickDir','out');
                axis tight
                xlim(cfg.epoch.toi{imarker});
                
                %save plot
                if ~(exist (cfg.imagesavedir)==7)
                    mkdir(cfg.imagesavedir);
                    fprintf('Create folder %s\n',cfg.imagesavedir);
                end
                
                if ~(exist (fullfile(cfg.imagesavedir,'remove_trials'))==7)
                    mkdir(fullfile(cfg.imagesavedir,'remove_trials'));
                    fprintf('Create folder %s\n',fullfile(cfg.imagesavedir,'remove_trials'));
                end
                
                %rename prefix in case of "merge" data{ipart}{imarker}
                if isfield(cfg, 'merge')
                    if cfg.merge == true
                        if ipart > 1 && ipart == length(cfg.directorylist) %last part = merge (except if only one part, nothing to merge)
                            cfg.prefix = [initial_prefix, 'MERGED-'];
                        else
                            cfg.prefix = [initial_prefix, cfg.directorylist{ipart}{:}, '-'];
                        end
                    end
                end
                
                fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix, cfg.LFP.name{imarker}, '_remove_trials_',method,'_',markerstart,'_',markerend]),'-r600');
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix, cfg.LFP.name{imarker}, '_remove_trials_',method,'_',markerstart,'_',markerend]),'-r600');
                close all
                
                
            end
            
            
            %% Remove trials : 
            data{ipart}{imarker}.time       = data{ipart}{imarker}.time(~trialremoved);
            data{ipart}{imarker}.trial      = data{ipart}{imarker}.trial(~trialremoved);
            data{ipart}{imarker}.trialinfo  = data{ipart}{imarker}.trialinfo(~trialremoved,:);
        end
    end %imarker
end %ipart

end

