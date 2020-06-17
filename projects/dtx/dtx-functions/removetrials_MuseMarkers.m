function [data, MuseStruct] = removetrials_MuseMarkers(cfg, data, MuseStruct, partlist, labellist)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function removetrials_MuseMarkers
%
% ### INPUT :
%
% # cfg fields to define trials (cf. readLFP.m or readSpykingCircus.m for
% more details) : cfg.muse.startend, cfg.epoch.toi
%
% # cfg fields to plot data if required : cfg.LFP.electrodetoplot,
% cfg.imagesavedir, cfg.prefix, cfg.name.
%
% # Specific cfg fields :
% cfg.method :
% 'remove'           : remove all trials which intersect with the period
% 'keep'             : keep only trials which intersect with the period
% cfg.markerstart/end: name of the marker which define the begin and the
%                      end of the period.
% cfg.indexstart/end : if need to take the +1 marker, i.e. for a period
%                      defined between samemarker(i) and samemarker(i+1).
%                      Otherwise let it as 0.
% cfg.timefrombegin/end  : relative time from cfg.markerstart and cfg.markerend, 
%                          in seconds, to define the period.
% cfg.plotdata       : 'yes' or 'no'. Either to plot the results or not.
%
% # Data structures and infos :
% data               : LFP data or spike data, respectively obtained from
%                      readLFP.m and readSpykingCircus.m
% MuseStruct         : structure obtained by readMuseMarkers.m. It contains 
%                      all the time information of the marker position in 
%                      the data.
% partlist           : list of parts of data to process. Can be 'all'.
% labellist          : list of group of trials to process ('ilabel' in
%                      data{ipart}{ilabel}).  Can be 'all'.
%
% ### OUTPUT :
% data               : returned without the selected trials.
% MuseStruct         : no modification of the marker timings. Add a field
%                      with the indexes of the removed trials in each
%                      MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel},1)
%
% Paul Baudin
% paul.baudin@live.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*****************************************************\n');
fprintf('Remove trials :\n');
fprintf('- %s all trials which intersect the defined period.\n' , cfg.method);
fprintf('- period begins with each marker %s(i+%d)%+gs\n'       , cfg.markerstart, cfg.indexstart, cfg.timefrombegin);
fprintf('- period ends with each marker %s(i+%d)%+gs\n'         , cfg.markerend  , cfg.indexend  , cfg.timefromend);
fprintf('*****************************************************\n');

if strcmp(partlist, 'all')
    partlist = 1:size(data, 2);
end

for ipart = partlist
    
    fprintf('For part %d\n', ipart);
    
    if strcmp(labellist, 'all')
        labellist = 1:size(data{ipart}, 2);
    end
    
    for ilabel = labellist
        
        fprintf('*** %s ***\n', cfg.name{ilabel});
        
        if ~isempty(data{ipart}{ilabel})
            
            %Trialinfos are different between LFP and spike data :
            [type, ~] = ft_datatype(data{ipart}{ilabel});
            if strcmp(type, 'raw')
                dircolumn       = 3;
                trialnrcolumn   = 1;
                Fs              = data{ipart}{ilabel}.fsample;
            elseif strcmp(type, 'spike')
                dircolumn       = 8;
                trialnrcolumn   = 2;
                Fs              = data{ipart}{ilabel}.hdr.Fs;
            else
                error('Data is unsupported format. Must be spike data or raw data\n');
            end
            
            if ~isempty(data{ipart}{ilabel})
                
                n_trials        = 0;
                trialremoved    = false(1,size(data{ipart}{ilabel}.trialinfo,1));
                
                for idir = 1:length(MuseStruct{ipart})
                    nr_notloaded(idir) = 0;
                    
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{ilabel,1})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}), 'synctime')
                            
                            %keep index of removed trials
                            MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).trialremoved = false(1,length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime));
                            
                            % if trial begin before first sample or end after last
                            % sample : is not loaded with readLFP or readSpykingCircus :
                            for itrialMuse = 1:length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime)
                                if ~ismember(itrialMuse, data{ipart}{ilabel}.trialinfo((data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir),trialnrcolumn))
                                    nr_notloaded(idir) = nr_notloaded(idir)+1;
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).trialremoved(itrialMuse) = true;
                                    fprintf('%s %s part %d dir %d : trial %d not loaded\n',cfg.prefix(1:end-1),cfg.muse.startend{ilabel,1}, ipart, idir, itrialMuse);
                                end
                            end
                            
                            % check if there is an equal amount of start and end markers
                            if isfield(MuseStruct{ipart}{idir}.markers,cfg.markerstart)
                                if isfield(MuseStruct{ipart}{idir}.markers.(cfg.markerstart), 'synctime')
                                    if size(MuseStruct{ipart}{idir}.markers.(cfg.markerstart).synctime,2)-size(MuseStruct{ipart}{idir}.markers.(cfg.markerend).synctime,2) == 0
                                        
                                        trials_idx = find(data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir)';
                                        
                                        for itrial = trials_idx
                                            itrialMuse = data{ipart}{ilabel}.trialinfo(itrial,trialnrcolumn); %trial of loaded data
                                            
                                            %Redefine trials as defined by readLFP and readSpykingCircus
                                            ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime(itrialMuse) * Fs);
                                            idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime * Fs) >= ss,1,'first');
                                            es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime(idx) * Fs);
                                            if ~isempty(es)
                                                trialinterval_begin     = ss + cfg.epoch.toi{ilabel}(1)*Fs;
                                                trialinterval_end       = es + cfg.epoch.toi{ilabel}(2)*Fs;
                                                trialinterval           = round(trialinterval_begin : trialinterval_end);
                                                
                                                %Define markerstart and markerend periods
                                                for isearchedinterval = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.markerstart).synctime,2)
                                                    %check that we don't exceed the index values
                                                    if isearchedinterval + cfg.indexstart <= length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime)
                                                        if isearchedinterval + cfg.indexend <= length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime)
                                                            
                                                            searchedinterval_begin       = (MuseStruct{ipart}{idir}.markers.(cfg.markerstart).synctime(isearchedinterval + cfg.indexstart) + cfg.timefrombegin) * Fs;
                                                            searchedinterval_end         = (MuseStruct{ipart}{idir}.markers.(cfg.markerend).synctime(isearchedinterval + cfg.indexend) + cfg.timefromend) * Fs;
                                                            
                                                        else %first marker exists but last marker exceeds file size
                                                            searchedinterval_end = seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime) * Fs;
                                                        end
                                                        
                                                        searchedinterval             = round(searchedinterval_begin : searchedinterval_end);
                                                        
                                                        if intersect(trialinterval,searchedinterval) %intervals are expressed in sample points
                                                            %fprintf('Found searched-trial in part %d trial %d (MuseStruct dir %d, marker %s n°%d)\n',ipart ,itrial,idir, cfg.muse.startend{ilabel,1},itrialMuse);
                                                            trialremoved(itrial) = 1;
                                                            MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).trialremoved(itrialMuse) = true;
                                                        end
                                                    end
                                                end %iserchedinterval
                                                
                                            end %~isempty(es)
                                        end %itrial
                                    else
                                        error('Not the same amout of cfg.markerstart (%s) and cfg.markerend (%s) found in %s \n',cfg.markerstart,cfg.markerend,MuseStruct{ipart}{idir}.directory);
                                        
                                    end
                                end
                            end
                            if strcmp(cfg.method,'keep')
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).trialremoved = ~MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).trialremoved;
                            end
                        end
                    end
                end %idir
                
                if strcmp(cfg.method, 'keep')
                    trialremoved = ~trialremoved;
                end
                
                if sum(trialremoved) > 0
                    fprintf('Remove %d searched-trials\n', sum(trialremoved));
                else
                    fprintf('No searched-trials found : all trials are kept\n');
                end
                
                
                %% plot data only if is LFP               
                if strcmp(cfg.plotdata, 'yes') && strcmp(type, 'raw')
                    
                    %h automatic setting :
                    for itrial = 1 : length(data{ipart}{ilabel}.trial)
                        t_h_1 = find(data{ipart}{ilabel}.time{itrial} >= -0.5, 1, 'first');
                        t_h_2 = find(data{ipart}{ilabel}.time{itrial} >= 0.5, 1, 'first');
                        h_temp(itrial) = max(data{ipart}{ilabel}.trial{itrial}(1,t_h_1:t_h_2)); %max between -0.5s and 0.5s.
                    end
                    
                    h = abs(mean(h_temp)) * 2;
                    
                    cfgtemp = [];
                    cfgtemp.channel = cfg.LFP.electrodetoplot{ilabel};
                    data_plot = ft_selectdata(cfgtemp, data{ipart}{ilabel});
                    
                    for idir = unique(data_plot.trialinfo(:,dircolumn))'
                        fig = figure;
                        hold;
                        
                        trials_idx = find(data_plot.trialinfo(:,dircolumn)==idir)';
                        for itrial = trials_idx
                            if ~trialremoved(itrial)
                                color = 'k';
                            else
                                color = 'r';
                            end
                            plot(data_plot.time{itrial}, data_plot.trial{itrial}+(n_trials+1)*h-itrial*h, 'color', color);
                        end
                        
                        xlabel(sprintf('Time from %s (s)', cfg.muse.startend{ilabel,1}),'Interpreter','none', 'Fontsize',15);
                        ylabel('Number of trials', 'Fontsize',15);
                        if nr_notloaded(idir) == 0
                            title(sprintf('%s chan %s \nAll trials are loaded \n%d searched-trials removed from %d trials \nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{ilabel}, data_plot.label{1}, sum(trialremoved), length(data{ipart}{ilabel}.trial), cfg.markerstart,cfg.indexstart, cfg.timefrombegin, cfg.markerend,cfg.indexend, cfg.timefromend),'Fontsize',20,'Interpreter','none');
                        else
                            title(sprintf('%s chan %s \n%d trials are not loaded (too close from begin or end of data) \n%d searched-trials trials removed from %d trials\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{ilabel}, data_plot.label{1}, nr_notloaded(idir), sum(trialremoved), length(data{ipart}{ilabel}.trial), cfg.markerstart,cfg.indexstart, cfg.timefrombegin, cfg.markerend,cfg.indexend, cfg.timefromend),'Fontsize',20,'Interpreter','none');
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
                        xlim(cfg.epoch.toi{ilabel});
                        
                        %save plot
                        if ~(exist(cfg.imagesavedir)==7)
                            mkdir(cfg.imagesavedir);
                            fprintf('Create folder %s\n',cfg.imagesavedir);
                        end
                        
                        if ~(exist(fullfile(cfg.imagesavedir,'remove_trials'))==7)
                            mkdir(fullfile(cfg.imagesavedir,'remove_trials'));
                            fprintf('Create folder %s\n',fullfile(cfg.imagesavedir,'remove_trials'));
                        end
                        
                        fig.Renderer = 'Painters'; % Else pdf is saved to bitmap
                        set(fig,'PaperOrientation','landscape');
                        set(fig,'PaperUnits','normalized');
                        set(fig,'PaperPosition', [0 0 1 1]);
                        print(fig, '-dpdf', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', cfg.name{ilabel}, '_remove_trials_',cfg.method,'_',cfg.markerstart,'_',cfg.markerend,'_',num2str(idir)]),'-r600');
                        print(fig, '-dpng', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', cfg.name{ilabel}, '_remove_trials_',cfg.method,'_',cfg.markerstart,'_',cfg.markerend,'_',num2str(idir)]),'-r600');
                        close all
                    end
                    
                end
                
                %% Remove trials :
                
                cfgtemp = [];
                cfgtemp.trials = ~trialremoved;
                
                if strcmp(type, 'raw')
                    data{ipart}{ilabel} = ft_selectdata(cfgtemp, data{ipart}{ilabel});
                elseif strcmp(type, 'spike')
                    data{ipart}{ilabel} = ft_spike_select(cfgtemp, data{ipart}{ilabel});
                end
                
                data{ipart}{ilabel}.trialremoved = trialremoved;
                
            end
        end %~isempty data{ipart}{ilabel}
    end %ilabel
end %ipart

end
