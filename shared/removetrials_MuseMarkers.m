function [data, MuseStruct] = removetrials_MuseMarkers(cfg, data, MuseStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [data, MuseStruct] = removetrials_MuseMarkers(cfg, data, MuseStruct)
% 
% Search for time periods in data, defined by Muse markers. Then, search if
% those periods overlap trials of data. Remove or keep only the trials
% which intersect.
% Data can be Fieldtrip spike data or Fieldtrip raw data (lfp).
% Default options are set to remove trials which intersect BAD__START__ and
% BAD__END markers.
%
% ### INPUT :
%
% # Data structures :
% data               : LFP data or spike data, respectively obtained from
%                      readLFP.m and readSpikeTrials_MuseMarkers.m. Can be
%                      obtained by another function if the trialinfo is
%                      consistent.
% MuseStruct         : structure obtained by readMuseMarkers.m. It contains
%                      all the time information of the marker position in
%                      the data.
%
% # cfg fields (all with default values):
% cfg.remove.method :
%   'remove'         : remove all trials which intersect with the period
%   'keep'           : keep only trials which intersect with the period
%   Default = 'remove'.
% cfg.remove.markerstart    : name of the marker which define the beginning of the
%                             period. Default = 'BAD__START__'.
% cfg.remove.markerend      : name of the marker which define the the end of the
%                             period. Default = 'BAD__END__'.
% cfg.remove.indexstart     : if need to take the +1 marker, i.e. for a period
%                             defined between samemarker(i) and samemarker(i+1).
%                             Otherwise let it as 0. Default = 0.
% cfg.remove.indexend       : if need to take the +1 marker, i.e. for a period
%                             defined between samemarker(i) and samemarker(i+1).
%                             Otherwise let it as 0. Default = 0.
% cfg.remove.timefrombegin  : relative time from cfg.remove.markerstart in seconds, to
%                             define the period to search. Default = 0.
% cfg.remove.timefromend    : relative time from cfg.remove.markerend in seconds, to
%                             define the period to search. Default = 0.
% cfg.remove.searchdirection : 'before', 'after', 'all'. Default = 'all'.
% Where the begin of the searched trial has to be compared to the begin of
% the data trial
% cfg.remove.keepindexes    : 'yes' or 'no', whether to keep indexes of removed
%                             trials in MuseStruct. This option use cfg.remove.muse.startend,
%                             so the trials have to be defined by readLFP or
%                             readSpikeTrials_MuseMarkers. If not, set this as
%                             'no'. Default = 'no'.
% cfg.remove.plotdata       : 'yes' or 'no'. Either to plot the results or not.
%                             Default = 'no'.
% cfg.remove.part_list      : list of parts to analyse. Can be an array of
%                             integers, or 'all'. Default = 'all'.
% cfg.remove.label_list     : list of groups of trials to analyse. Can be an array
%                             of integers, or 'all'. Default = 'all'.
%
% # cfg fields to plot data if required (necessary if cfg.remove.plotdata = 'yes') :
% cfg.remove.electrodetoplot : name of the electrode used for the plot.
%                              Cell array with one cell per data label
% cfg.imagesavedir        : where to save the images
% cfg.prefix              : prefix appended to the nama of output images
% cfg.name                : used in the title of the image, and in the file
%                           name
% cfg.epoch.toi           : x limits for the plot
%
% ### OUTPUT :
% data               : returned without the selected trials.
% MuseStruct         : no modification of the marker timings. Add a field
%                      with the indexes of the removed trials in each
%                      MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel},1),
%                      if cfg.remove.keepindexes = 'yes'.
%
% Note : some trials are not loaded in the data structures, ie there is a
% marker in MuseStruct, but no associated trial in data. The reason for
% that could be :
% - the trial starts or ends outside the file
% - the alignment failed for this marker, so it is removed before loading
%   and cutting data into trials.
%
% Paul Baudin
% paul.baudin@live.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the default cfg options
cfg.remove                  = ft_getopt(cfg,'remove',[]);
cfg.remove.method           = ft_getopt(cfg.remove, 'method'        , 'remove');
cfg.remove.markerstart      = ft_getopt(cfg.remove, 'markerstart'   , 'BAD__START__');
cfg.remove.markerend        = ft_getopt(cfg.remove, 'markerend'     , 'BAD__END__');
cfg.remove.indexstart       = ft_getopt(cfg.remove, 'indexstart'    , 0);
cfg.remove.indexend         = ft_getopt(cfg.remove, 'indexend'      , 0);
cfg.remove.timefrombegin    = ft_getopt(cfg.remove, 'timefrombegin' , 0);
cfg.remove.timefromend      = ft_getopt(cfg.remove, 'timefromend'   , 0);
cfg.remove.searchdirection  = ft_getopt(cfg.remove, 'searchdirection', 'all');
cfg.remove.keepindexes      = ft_getopt(cfg.remove, 'keepindexes'   , 'no');
cfg.remove.plotdata         = ft_getopt(cfg.remove, 'plotdata'      , 'no');
cfg.remove.part_list        = ft_getopt(cfg.remove, 'part_list'     , 'all');
cfg.remove.label_list       = ft_getopt(cfg.remove, 'label_list'    , 'all');


if isempty(data), fprintf('removetrials_MuseMarkers : Data is empty, nothing is done\n'); return, end

fprintf('\n*****************************************************\n');
fprintf('Remove trials :\n');
fprintf('- %s all trials which intersect the defined period.\n' , cfg.remove.method);
fprintf('- period begins with each marker %s(i+%d)%+gs\n'       , cfg.remove.markerstart, cfg.remove.indexstart, cfg.remove.timefrombegin);
fprintf('- period ends with each marker %s(i+%d)%+gs\n'         , cfg.remove.markerend  , cfg.remove.indexend  , cfg.remove.timefromend);
fprintf('*****************************************************\n');

if strcmp(cfg.remove.part_list,'all')
    cfg.remove.part_list = 1:size(MuseStruct,2);
end

for ipart = cfg.remove.part_list
    
    fprintf('For part %d\n', ipart);
    
    if strcmp(cfg.remove.label_list, 'all')
        cfg.remove.label_list = 1:size(data{ipart}, 2);
    end
    
    for ilabel = cfg.remove.label_list
        
        fprintf('*** For ilabel = %d ***\n', ilabel);
        
        if isempty(data{ipart}{ilabel})
            continue
        end
        
        %Trialinfos are different between LFP and spike data :
        [type, ~] = ft_datatype(data{ipart}{ilabel});
        if strcmp(type, 'raw')
            dircolumn       = 3;
            trialnrcolumn   = 1;
            startcolumn     = 5;
            endcolumn       = 6;
            Fs              = data{ipart}{ilabel}.fsample;
        elseif strcmp(type, 'spike')
            dircolumn       = 8;
            trialnrcolumn   = 2;
            startcolumn     = 3;
            endcolumn       = 4;
            Fs              = data{ipart}{ilabel}.hdr.Fs;
        else
            error('Data is unsupported format. Must be spike data or raw data\n');
        end
        
        if ~isempty(data{ipart}{ilabel})
                       
            for idir = 1:length(MuseStruct{ipart})
                
                trialremoved{idir}    = false(1,sum(data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir));
                nr_notloaded(idir) = 0;
                
                %keep info on removed trials, if based on muse markers
                if strcmp(cfg.remove.keepindexes, 'yes')
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.remove.muse.startend{ilabel,1})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}), 'synctime')
                            
                            %keep index of removed trials
                            MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).trialremoved = false(1,length(MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).synctime));
                            
                            % if trial begin before first sample or end after last
                            % sample : is not loaded with readLFP or readSpykingCircus :
                            for itrialMuse = 1:length(MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).synctime)
                                if ~ismember(itrialMuse, data{ipart}{ilabel}.trialinfo((data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir),trialnrcolumn))
                                    nr_notloaded(idir) = nr_notloaded(idir)+1;
                                    MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).trialremoved(itrialMuse) = true;
                                    fprintf('%s %s part %d dir %d : trial %d not loaded\n',cfg.prefix(1:end-1),cfg.remove.muse.startend{ilabel,1}, ipart, idir, itrialMuse);
                                end
                            end
                        end
                    end
                end
                
                % check if there is an equal amount of start and end markers
                if ~isfield(MuseStruct{ipart}{idir}.markers,cfg.remove.markerstart)
                    continue
                end
                if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerstart), 'synctime')
                    continue
                end
                if size(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerstart).synctime,2)-size(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerend).synctime,2) ~= 0
                    error('Not the same amout of cfg.remove.markerstart (%s) and cfg.remove.markerend (%s) found in %s \n',cfg.remove.markerstart,cfg.remove.markerend,MuseStruct{ipart}{idir}.directory);
                end
                trials_idx = find(data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir)';
                
                %dir_onset for spike data, because spike data are based on concatenated files, whereas lfp data are not
                if strcmp(type, 'spike')
                    dir_onset = data{ipart}{ilabel}.trialinfo(find(data{ipart}{ilabel}.trialinfo(:,dircolumn)==idir,1,'first'),9);
                else
                    dir_onset = 0;
                end
                i=0;
                
                for itrial = trials_idx
                    i = i+1;
                    
                    %waitbar
                    try eval(['fprintf(repmat(''\b'', 1, length(waitbar)));']); end %error the first time. \b to remove previous text
                    waitbar = sprintf('dir %d/%d : trial %d/%d',idir, length(MuseStruct{ipart}),i, length(trials_idx)); 
                    fprintf('%s', waitbar);
                    
                    
                    trialinterval_begin     = data{ipart}{ilabel}.trialinfo(itrial, startcolumn) - dir_onset;
                    trialinterval_end       = data{ipart}{ilabel}.trialinfo(itrial, endcolumn) - dir_onset;
                    trialinterval           = round(trialinterval_begin : trialinterval_end);
                    
                    %Define markerstart and markerend periods
                    for isearchedinterval = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerstart).synctime,2)
                        %check that we don't exceed the index values
                        if isearchedinterval + cfg.remove.indexstart <= length(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerstart).synctime)
                            if isearchedinterval + cfg.remove.indexend <= length(MuseStruct{ipart}{idir}.markers.(cfg.remove.markerend).synctime)
                                
                                searchedinterval_begin       = (MuseStruct{ipart}{idir}.markers.(cfg.remove.markerstart).synctime(isearchedinterval + cfg.remove.indexstart) + cfg.remove.timefrombegin) * Fs;
                                searchedinterval_end         = (MuseStruct{ipart}{idir}.markers.(cfg.remove.markerend).synctime(isearchedinterval + cfg.remove.indexend) + cfg.remove.timefromend) * Fs;
                                
                            else %first marker exists but last marker exceeds file size
                                searchedinterval_end = seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime) * Fs;
                            end
                            
                            if strcmp(cfg.remove.searchdirection, 'before') && searchedinterval_begin > trialinterval_begin
                                continue
                            elseif strcmp(cfg.remove.searchdirection, 'after') && searchedinterval_begin < trialinterval_begin
                                continue
                            end
                                
                            
                            searchedinterval             = round(searchedinterval_begin : searchedinterval_end);
                            
                            if intersect(trialinterval,searchedinterval) %intervals are expressed in sample points
                                %                                     fprintf('Found searched-trial in part %d trial %d (MuseStruct dir %d, marker %s n°%d)\n',ipart ,itrial,idir, cfg.remove.muse.startend{ilabel,1},itrialMuse);
                                trialremoved{idir}(i) = true;
                                if strcmp(cfg.remove.keepindexes, 'yes') %keep indexes of removed trials
                                    itrialMuse = data{ipart}{ilabel}.trialinfo(itrial,trialnrcolumn);%trial of loaded data
                                    MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).trialremoved(itrialMuse) = true;
                                end
                            end
                        end
                    end %iserchedinterval
                    
                end %itrial
                
                if strcmp(cfg.remove.method,'keep')&&strcmp(cfg.remove.keepindexes, 'yes')
                    MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).trialremoved = ~MuseStruct{ipart}{idir}.markers.(cfg.remove.muse.startend{ilabel,1}).trialremoved;
                end
            end %idir
            
            
            if strcmp(cfg.remove.method, 'keep')
                for idir = 1:size(trialremoved,2)
                    trialremoved{idir} = ~trialremoved{idir};
                end
            end
            
            if sum([trialremoved{:}]) > 0
                fprintf('\nRemove %d trials\n', sum([trialremoved{:}]));
            else
                fprintf('\nNo searched-trials found : all trials are kept\n');
            end
            
            
            %% plot data only if is LFP
            if strcmp(cfg.remove.plotdata, 'yes') && strcmp(type, 'raw')
                
                %h automatic setting :
                for itrial = 1 : length(data{ipart}{ilabel}.trial)
                    t_h_1 = find(data{ipart}{ilabel}.time{itrial} >= -0.5, 1, 'first');
                    t_h_2 = find(data{ipart}{ilabel}.time{itrial} >= 0.5, 1, 'first');
                    h_temp(itrial) = max(data{ipart}{ilabel}.trial{itrial}(1,t_h_1:t_h_2)); %max between -0.5s and 0.5s.
                end
                
                h = abs(mean(h_temp)) * 2;
                
                cfgtemp = [];
                cfgtemp.channel = cfg.remove.electrodetoplot{ilabel};
                data_plot = ft_selectdata(cfgtemp, data{ipart}{ilabel});
                
                n_trials_total = 0;
                
                for idir = unique(data_plot.trialinfo(:,dircolumn))'
                    fig = figure;
                    hold;
                    n_trials_dir = size(trialremoved{idir},2);
                    
                    for itrial = 1:n_trials_dir
                        if ~trialremoved{idir}(itrial)
                            color = 'k';
                        else
                            color = 'r';
                        end
                        plot(data_plot.time{itrial+n_trials_total}, data_plot.trial{itrial+n_trials_total}+(n_trials_dir+1)*h-itrial*h, 'color', color);
                    end
                    n_trials_total = n_trials_total+size(trialremoved{idir},2);
                    
                    xlabel('Time (s)', 'Fontsize',15);
                    ylabel('Number of trials', 'Fontsize',15);
                    if nr_notloaded(idir) == 0
                        title(sprintf('%s chan %s \nAll trials are loaded \n%d trials removed from %d\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{ilabel}, data_plot.label{1}, sum(trialremoved{idir}), size(trialremoved{idir},2), cfg.remove.markerstart,cfg.remove.indexstart, cfg.remove.timefrombegin, cfg.remove.markerend,cfg.remove.indexend, cfg.remove.timefromend),'Fontsize',20,'Interpreter','none');
                    else
                        title(sprintf('%s chan %s \n%d trials are not loaded \n%d trials removed from %d\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', cfg.LFP.name{ilabel}, data_plot.label{1}, nr_notloaded(idir), sum(trialremoved{idir}), size(trialremoved{idir},2), cfg.remove.markerstart,cfg.remove.indexstart, cfg.remove.timefrombegin, cfg.remove.markerend,cfg.remove.indexend, cfg.remove.timefromend),'Fontsize',20,'Interpreter','none');
                    end
                    set(gca, 'YTickLabel', '','FontWeight','bold', 'Fontsize',15);
                    tick_interval = round(n_trials_dir/10);
                    if tick_interval == 0
                        tick_interval = 1;
                    end
                    yticks(h : h*tick_interval : n_trials_dir*h);
                    yticklabels(n_trials_dir : -tick_interval : 0);
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
                    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', cfg.name{ilabel}, '_remove_trials_',cfg.remove.method,'_',cfg.remove.markerstart,'_',cfg.remove.markerend,'_',num2str(idir)]),'-r600');
                    print(fig, '-dpng', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', cfg.name{ilabel}, '_remove_trials_',cfg.remove.method,'_',cfg.remove.markerstart,'_',cfg.remove.markerend,'_',num2str(idir)]),'-r600');
                    close all
                end
                
            end
            
            %% Remove trials :
            
            cfgtemp = [];
            cfgtemp.trials = [];
            for idir = 1:size(trialremoved,2)
                cfgtemp.trials = [cfgtemp.trials, ~trialremoved{idir}];
            end
            cfgtemp.trials = find(cfgtemp.trials);
            
            if strcmp(type, 'raw')
                data{ipart}{ilabel} = ft_selectdata(cfgtemp, data{ipart}{ilabel});
            elseif strcmp(type, 'spike')
                data{ipart}{ilabel} = ft_spike_select(cfgtemp, data{ipart}{ilabel});
            end
            
            data{ipart}{ilabel}.trialremoved = trialremoved;
            
        end
    end %ilabel
end %ipart

end

