function [data, MuseStruct] = dtx_removeartefactsLFP(data, MuseStruct, cfg, plotdata)
%Remove artefacted trials in data{ipart}{imarker} according to "BAD" Muse markers
%MuseStruct is returned with a new 'hasartefact' field, but the markers are
%not removed.
%Add field : data.nr_removedtrials

% /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
%OLD : do not work anymore. Is based on marker name indicated in
%cfg.LFP.name, it was not a good strategy. See the script removetrialLFP_marker
% /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
initial_prefix = cfg.prefix; %for merge


for ipart = 1:size(data, 2)
    for imarker = 1:size(data{ipart}, 2)
        if ~isempty(data{ipart}{imarker})
            
            n_trials        = 0;
            hasartefact     = false(1,size(data{ipart}{imarker}.trial,2));
            
            % for ipart = 1:length(MuseStruct)
            for idir = 1:length(MuseStruct{ipart})
                nr_notloaded(idir) = 0;
                
                if isfield(MuseStruct{ipart}{idir}.markers,cfg.LFP.name{imarker})
                    if isfield(MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}), 'synctime')
                        MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}).hasartefact = false(1,length(MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}).synctime));

                        % if trial begin before first sample or end after last
                        % sample : is not loaded with readLFP
                        
                        for itrialMuse = 1:length(MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}).synctime)
                            %if trial of the dir (data.trialinfo(:,3)==idir) does note have
                            %the same index (data.trialinfo(:,1) as Muse LFP marker
                            %then annotate as artefact
                            if ~ismember(itrialMuse, data{ipart}{imarker}.trialinfo((data{ipart}{imarker}.trialinfo(:,3)==idir),1))
                                nr_notloaded(idir) = nr_notloaded(idir)+1;
                                MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}).hasartefact(itrialMuse) = true;
                                fprintf('%s %s part %d dir %d : trial %d not loaded\n',cfg.prefix(1:end-1),cfg.LFP.name{imarker}, ipart, idir, itrialMuse);
                            end
                        end
                        
                        % check if there is an equal amount of BAD start and end markers
                        if isfield(MuseStruct{ipart}{idir}.markers,'BAD__START__')
                            if size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) == 0
                                
%                                 if isNeuralynx
%                                     error('I have to set the header reading code, not done yet');%hdr = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['.ncs']));
%                                 elseif isMicromed
%                                     hdr = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir}, '.TRC']));
%                                 elseif isBrainvision
%                                     hdr = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir},'.eeg']));
%                                 end
                                
                                %convert idir{itrial} to data{ipart}{imarker}{itrial} :
                                for itrial = n_trials+1 : n_trials+sum(data{ipart}{imarker}.trialinfo(:,3)==idir)
                                    itrialMuse = data{ipart}{imarker}.trialinfo(itrial,1); %trial of loaded data
                                    
                                    trialinterval = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrialMuse)*data{ipart}{imarker}.fsample : MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(itrialMuse)*data{ipart}{imarker}.fsample);
                                    
                                    for iartefact = 1 : size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)
                                        artefactinterval = round(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(iartefact)*data{ipart}{imarker}.fsample : MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(iartefact)*data{ipart}{imarker}.fsample);
                                        if intersect(trialinterval,artefactinterval) %intervals are expressed in sample points
                                            fprintf('Found artefact in part %d trial %d for marker %s \n',ipart ,itrial, cfg.LFP.name{imarker});
                                            hasartefact(itrial) = 1;
                                            MuseStruct{ipart}{idir}.markers.(cfg.LFP.name{imarker}).hasartefact(itrialMuse) = true;
                                        end
                                    end
                                    
                                    
                                end
                                
                                
                            else
                                error('Not the same amout of BAD starts than ends markers found in %s \n',MuseStruct{ipart}{idir}.directory);
                            end
                        end
                    end
                end
                n_trials = n_trials + sum(data{ipart}{imarker}.trialinfo(:,3)==idir); %nr of trials per dir
            end %idir
            %end
            
            if sum(hasartefact) > 0
                fprintf('===> Found %d artefacts in part %d for marker %s\n', sum(hasartefact), ipart, cfg.LFP.name{imarker});
            else
                fprintf('No artefacts found in part %d for marker %s\n',ipart,cfg.LFP.name{imarker});
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
                    if ~hasartefact(itrial)
                        color = 'k';
                    else
                        color = 'r';
                    end
                    plot(data_plot.time{itrial}, data_plot.trial{itrial}+(n_trials+1)*h-itrial*h, 'color', color);
                end
                
                xlabel(sprintf('Time from %s (s)', cfg.LFP.name{imarker}),'Interpreter','none', 'Fontsize',15);
                ylabel('Number of trials', 'Fontsize',15);
                if sum(nr_notloaded) == 0
                    title(sprintf('%s chan %s \nAll trials are loaded \n%d artefacted trials removed from %d trials', cfg.LFP.name{imarker}, data_plot.label{1}, sum(hasartefact), length(data{ipart}{imarker}.trial)),'Fontsize',20,'Interpreter','none');
                else
                    title(sprintf('%s chan %s \n%d trials are not loaded (too close from begin or end of data) \n%d artefacted trials removed from %d trials', cfg.LFP.name{imarker}, data_plot.label{1}, sum(nr_notloaded), sum(hasartefact), length(data{ipart}{imarker}.trial)),'Fontsize',20,'Interpreter','none');
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
                    fprintf('Create folder %s',cfg.imagesavedir);
                end
                
                if ~(exist (fullfile(cfg.imagesavedir,'remove_artefacts'))==7)
                    mkdir(fullfile(cfg.imagesavedir,'remove_artefacts'));
                    fprintf('Create folder %s',fullfile(cfg.imagesavedir,'remove_artefacts'));
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
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'remove_artefacts',[cfg.prefix, cfg.LFP.name{imarker}, '_remove_artefacts']),'-r600');
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'remove_artefacts',[cfg.prefix, cfg.LFP.name{imarker}, '_remove_artefacts']),'-r600');
                close all
                
                
            end
            
            
            %% Remove artefacts : 
            data{ipart}{imarker}.time       = data{ipart}{imarker}.time(~hasartefact);
            data{ipart}{imarker}.trial      = data{ipart}{imarker}.trial(~hasartefact);
            data{ipart}{imarker}.trialinfo  = data{ipart}{imarker}.trialinfo(~hasartefact,:);
        end
    end %imarker
end %ipart

end

