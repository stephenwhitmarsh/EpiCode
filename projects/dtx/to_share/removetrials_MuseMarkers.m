function [data, MuseStruct] = removetrials_MuseMarkers(cfg, data, MuseStruct, force)

% function [data, MuseStruct] = rmtrialstrials_MuseMarkers(cfg, data, MuseStruct)
%
% Search for time periods in data, defined by Muse markers. Then, search if
% those periods overlap trials of data. rmtrials or keep only the trials
% which intersect.
% Data can be Fieldtrip spike data or Fieldtrip raw data (lfp).
% Default options are set to rmtrials trials which intersect BAD__START__ and
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
% cfg.rmtrials.method :
%   'rmtrials'         : remove all trials which intersect with the period
%   'keep'           : keep only trials which intersect with the period
%   Default = 'rmtrials'.
% cfg.rmtrials.markerstart    : name of the marker which define the beginning of the
%                             period. Default = 'BAD__START__'.
% cfg.rmtrials.markerend      : name of the marker which define the the end of the
%                             period. Default = 'BAD__END__'.
% cfg.rmtrials.indexstart     : if need to take the +1 marker, i.e. for a period
%                             defined between samemarker(i) and samemarker(i+1).
%                             Otherwise let it as 0. Default = 0.
% cfg.rmtrials.indexend       : if need to take the +1 marker, i.e. for a period
%                             defined between samemarker(i) and samemarker(i+1).
%                             Otherwise let it as 0. Default = 0.
% cfg.rmtrials.timefrombegin  : relative time from cfg.rmtrials.markerstart in seconds, to
%                             define the period to search. Default = 0.
% cfg.rmtrials.timefromend    : relative time from cfg.rmtrials.markerend in seconds, to
%                             define the period to search. Default = 0.
% cfg.rmtrials.searchdirection : 'before', 'after', 'all'. Default = 'all'.
% Where the begin of the searched trial has to be compared to the begin of
% the data trial
% cfg.rmtrials.keepindexes    : 'yes' or 'no', whether to keep indexes of rmtrialsd
%                             trials in MuseStruct. This option use cfg.rmtrials.muse.startend,
%                             so the trials have to be defined by readLFP or
%                             readSpikeTrials_MuseMarkers. If not, set this as
%                             'no'. Default = 'no'.
% cfg.rmtrials.plotdata       : 'yes' or 'no'. Either to plot the results or not.
%                             Default = 'no'.
% cfg.rmtrials.part_list      : list of parts to analyse. Can be an array of
%                             integers, or 'all'. Default = 'all'.
% cfg.rmtrials.label_list     : list of groups of trials to analyse. Can be an array
%                             of integers, or 'all'. Default = 'all'.
% cfg.rmtrials.write          : 'yes' (default) or 'no', wheter to write
%                             output data on disk or not.
%
% # cfg fields to plot data if required (necessary if cfg.rmtrials.plotdata = 'yes') :
% cfg.rmtrials.electrodetoplot : name of the electrode used for the plot.
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
%                      with the indexes of the rmtrialsd trials in each
%                      MuseStruct{ipart}{idir}.markers.(startmarker,1),
%                      if cfg.rmtrials.keepindexes = 'yes'.
%
% Note : some trials are not loaded in the data structures, ie there is a
% marker in MuseStruct, but no associated trial in data. The reason for
% that could be :
% - the trial starts or ends outside the file
% - the alignment failed for this marker, so it is rmtrialsd before loading
%   and cutting data into trials.
%


fname = fullfile(cfg.datasavedir, [cfg.prefix, 'SpikeTrials_MuseMarkers_WithoutArtefacts.mat']);

if exist(fname) && force == false
    fprint('Loading precomputed removal of artefacts\n%s\n', fname);
    data = load(fname);
    return
end

% get the default cfg options
cfg.rmtrials                  = ft_getopt(cfg,'rmtrials',[]);
cfg.rmtrials.method           = ft_getopt(cfg.rmtrials, 'method'         , 'rmtrials');
cfg.rmtrials.markerstart      = ft_getopt(cfg.rmtrials, 'markerstart'    , 'BAD__START__');
cfg.rmtrials.markerend        = ft_getopt(cfg.rmtrials, 'markerend'      , 'BAD__END__');
cfg.rmtrials.indexstart       = ft_getopt(cfg.rmtrials, 'indexstart'     , 0);
cfg.rmtrials.indexend         = ft_getopt(cfg.rmtrials, 'indexend'       , 0);
cfg.rmtrials.timefrombegin    = ft_getopt(cfg.rmtrials, 'timefrombegin'  , 0);
cfg.rmtrials.timefromend      = ft_getopt(cfg.rmtrials, 'timefromend'    , 0);
cfg.rmtrials.searchdirection  = ft_getopt(cfg.rmtrials, 'searchdirection', 'all');
cfg.rmtrials.keepindexes      = ft_getopt(cfg.rmtrials, 'keepindexes'    , 'no');
cfg.rmtrials.plotdata         = ft_getopt(cfg.rmtrials, 'plotdata'       , 'no');
cfg.rmtrials.part_list        = ft_getopt(cfg.rmtrials, 'part_list'      , 'all');
cfg.rmtrials.label_list       = ft_getopt(cfg.rmtrials, 'label_list'     , 'all');
cfg.rmtrials.write            = ft_getopt(cfg.rmtrials, 'write'          , 'yes');


if isempty(data), fprintf('rmtrialstrials_MuseMarkers : Data is empty, nothing is done\n'); return, end

fprintf('\n*****************************************************\n');
fprintf('remove trials :\n');
fprintf('- %s all trials which intersect the defined period.\n' , cfg.rmtrials.method);
fprintf('- period begins with each marker %s(i+%d)%+gs\n'       , cfg.rmtrials.markerstart, cfg.rmtrials.indexstart, cfg.rmtrials.timefrombegin);
fprintf('- period ends with each marker %s(i+%d)%+gs\n'         , cfg.rmtrials.markerend  , cfg.rmtrials.indexend  , cfg.rmtrials.timefromend);
fprintf('*****************************************************\n');

if strcmp(cfg.rmtrials.part_list,'all')
    cfg.rmtrials.part_list = 1:size(MuseStruct,2);
end


for ipart = cfg.rmtrials.part_list
    
    fprintf('For part %d\n', ipart);
    
    if strcmp(cfg.rmtrials.label_list, 'all')
        cfg.rmtrials.label_list = 1:size(data{ipart}, 2);
    end
    
    for markername = string(fieldnames(data{ipart}))'
        
        fprintf('*** For %s ***\n', markername);
        
        if isempty(data{ipart}.(markername))
            continue
        end
                
        %check datatype, find Fs, and correct trialinfo (to ahvecompatibility between spike and raw fieldtrip structures)
        type = ft_datatype(data{ipart}.(markername));
        if strcmp(type, 'raw')
            Fs                                             = data{ipart}.(markername).fsample;
            data{ipart}.(markername).trialinfo.trialnr_dir = data{ipart}.(markername).trialinfo.trialnr;
        elseif strcmp(type, 'spike')
            Fs                                             = data{ipart}.(markername).hdr.Fs;
        else
            error('Data type ''%s'' is unsupported. Must be spike data or raw data\n', type);
        end
        
            ft_progress('init','text');
            
            for idir = 1:length(MuseStruct{ipart})
                
                trialdir              = data{ipart}.(markername).trialinfo.idir==idir;
                trialremoved{idir}    = false(1,sum(trialdir));
                nr_notloaded(idir)    = 0;
                
                %keep info on rmtrialsd trials, if based on muse markers
                if strcmp(cfg.rmtrials.keepindexes, 'yes')
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startmarker.(markername))
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)), 'synctime')
                            
                            %keep index of rmtrialsd trials
                            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialremoved = false(1,length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime));
                            
                            % if trial begin before first sample or end after last
                            % sample : is not loaded with readLFP or readSpykingCircus :
                            for itrialMuse = 1:length(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime)
                                if ~ismember(itrialMuse, data{ipart}.(markername).trialinfo.trialnr_dir(trialdir))
                                    nr_notloaded(idir) = nr_notloaded(idir)+1;
                                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialremoved(itrialMuse) = true;
                                    fprintf('%s %s part %d dir %d : trial %d not loaded\n',cfg.prefix(1:end-1),cfg.muse.startmarker.(markername), ipart, idir, itrialMuse);
                                end
                            end
                        end
                    end
                end
                
                % check if there is an equal amount of start and end markers
                if ~isfield(MuseStruct{ipart}{idir}.markers,cfg.rmtrials.markerstart)
                    continue
                end
                if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerstart), 'synctime')
                    continue
                end
                if size(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerstart).synctime,2)-size(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerend).synctime,2) ~= 0
                    error('Not the same amout of cfg.rmtrials.markerstart (%s) and cfg.rmtrials.markerend (%s) found in %s \n',cfg.rmtrials.markerstart,cfg.rmtrials.markerend,MuseStruct{ipart}{idir}.directory);
                end
                trials_idx = find(trialdir)';
                
                %dir_onset for spike data, because spike data are based on concatenated files, whereas lfp data are not
                if strcmp(type, 'spike')
                    dir_onset = data{ipart}.(markername).trialinfo.begsample(trials_idx(1));
                else
                    dir_onset = 0;
                end
                i=0;
                
                for itrial = trials_idx
                    i = i+1; %count the trials for waitbar
                    
                    %waitbar
                    ft_progress(idir/length(MuseStruct{ipart}), 'dir %d/%d : trial %d/%d', idir, length(MuseStruct{ipart}),i, length(trials_idx));
                    
                    trialinterval_begin     = data{ipart}.(markername).trialinfo.begsample(itrial) - dir_onset;
                    trialinterval_end       = data{ipart}.(markername).trialinfo.endsample(itrial) - dir_onset;
                    
                    %Define markerstart and markerend periods
                    for isearchedinterval = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerstart).synctime,2)
                        %check that we don't exceed the index values for begin
                        if isearchedinterval + cfg.rmtrials.indexstart > length(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerstart).synctime)
                            continue
                        end
                        searchedinterval_begin       = (MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerstart).synctime(isearchedinterval + cfg.rmtrials.indexstart) + cfg.rmtrials.timefrombegin) * Fs;
                        
                        %check that we don't exceed the index values for end
                        if isearchedinterval + cfg.rmtrials.indexend <= length(MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerend).synctime)
                            searchedinterval_end         = (MuseStruct{ipart}{idir}.markers.(cfg.rmtrials.markerend).synctime(isearchedinterval + cfg.rmtrials.indexend) + cfg.rmtrials.timefromend) * Fs;
                        else %first marker exists but last marker exceeds file size
                            searchedinterval_end = seconds(MuseStruct{ipart}{idir}.endtime - MuseStruct{ipart}{idir}.starttime) * Fs;
                        end
                        
                        if strcmp(cfg.rmtrials.searchdirection, 'before') & searchedinterval_begin > trialinterval_begin
                            continue
                        elseif strcmp(cfg.rmtrials.searchdirection, 'after') & searchedinterval_begin < trialinterval_begin
                            continue
                        end
                        
                        if searchedinterval_begin > trialinterval_begin & searchedinterval_begin > trialinterval_end
                            continue
                        elseif searchedinterval_end < trialinterval_begin & searchedinterval_end < trialinterval_end
                            continue
                        else
                            trialremoved{idir}(i) = true;
                            if strcmp(cfg.rmtrials.keepindexes, 'yes') %keep indexes of rmtrialsd trials
                                itrialMuse = data{ipart}.(markername).trialinfo.trialnr_dir(itrial);%trial of loaded data
                                MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialremoved(itrialMuse) = true;
                            end
                        end
                        
                    end %iserchedinterval
                    
                end %itrial
                
                
                if strcmp(cfg.rmtrials.method,'keep')&&strcmp(cfg.rmtrials.keepindexes, 'yes')
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialremoved = ~MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialremoved;
                end
            end %idir
            ft_progress('close');
            
            if strcmp(cfg.rmtrials.method, 'keep')
                for idir = 1:size(trialremoved,2)
                    trialremoved{idir} = ~trialremoved{idir};
                end
            end
            
            if sum([trialremoved{:}]) > 0
                fprintf('\nremove %d trials\n', sum([trialremoved{:}]));
            else
                fprintf('No searched-trials found : all trials are kept\n');
            end
            
            
            %% plot data only if is LFP
            if strcmp(cfg.rmtrials.plotdata, 'yes') && strcmp(type, 'raw')
                
                %h automatic setting :
                for itrial = 1 : length(data{ipart}.(markername).trial)
                    t_h_1 = find(data{ipart}.(markername).time{itrial} >= -0.5, 1, 'first');
                    t_h_2 = find(data{ipart}.(markername).time{itrial} >= 0.5, 1, 'first');
                    h_temp(itrial) = max(data{ipart}.(markername).trial{itrial}(1,t_h_1:t_h_2)); %max between -0.5s and 0.5s.
                end
                
                h = abs(mean(h_temp)) * 2;
                
                cfgtemp = [];
                cfgtemp.channel = cfg.rmtrials.electrodetoplot.(markername);
                data_plot = ft_selectdata(cfgtemp, data{ipart}.(markername));
                
                n_trials_total = 0;
                
                for idir = unique(data_plot.trialinfo.idir)'
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
                        title(sprintf('%s chan %s \nAll trials are loaded \n%d trials removed from %d\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', markername, data_plot.label{1}, sum(trialremoved{idir}), size(trialremoved{idir},2), cfg.rmtrials.markerstart,cfg.rmtrials.indexstart, cfg.rmtrials.timefrombegin, cfg.rmtrials.markerend,cfg.rmtrials.indexend, cfg.rmtrials.timefromend),'Fontsize',20,'Interpreter','none');
                    else
                        title(sprintf('%s chan %s \n%d trials are not loaded \n%d trials remove from %d\nSearch between %s(i+%d)%+gs and %s(i+%d)%+gs', markername, data_plot.label{1}, nr_notloaded(idir), sum(trialremoved{idir}), size(trialremoved{idir},2), cfg.rmtrials.markerstart,cfg.rmtrials.indexstart, cfg.rmtrials.timefrombegin, cfg.rmtrials.markerend,cfg.rmtrials.indexend, cfg.rmtrials.timefromend),'Fontsize',20,'Interpreter','none');
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
                    xlim(cfg.epoch.toi.(markername));
                    
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
                    print(fig, '-dpdf', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', convertStringsToChars(markername), '_remove_trials_',cfg.rmtrials.method,'_',cfg.rmtrials.markerstart,'_',cfg.rmtrials.markerend,'_',num2str(idir)]),'-r600');
                    print(fig, '-dpng', fullfile(cfg.imagesavedir,'remove_trials',[cfg.prefix,'part',num2str(ipart),'_', convertStringsToChars(markername), '_remove_trials_',cfg.rmtrials.method,'_',cfg.rmtrials.markerstart,'_',cfg.rmtrials.markerend,'_',num2str(idir)]),'-r600');
                    close all
                end
                
            end
            
            %% rmtrials trials :
            
            cfgtemp = [];
            cfgtemp.trials = [];
            for idir = 1:size(trialremoved,2)
                cfgtemp.trials = [cfgtemp.trials, ~trialremoved{idir}];
            end
            cfgtemp.trials = find(cfgtemp.trials);
            
            if strcmp(type, 'raw')
                data{ipart}.(markername) = ft_selectdata(cfgtemp, data{ipart}.(markername));
            elseif strcmp(type, 'spike')
                data{ipart}.(markername) = ft_spike_select(cfgtemp, data{ipart}.(markername));
            end
            
            data{ipart}.(markername).trialremoved = trialremoved;
            
    end %ilabel
end %ipart

%save data
if strcmp(cfg.rmtrials.write, 'yes')
    save(fname, 'data', '-v7.3');
end

end