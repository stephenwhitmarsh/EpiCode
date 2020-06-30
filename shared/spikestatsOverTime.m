function [stats, legend_out] = spikestatsOverTime(cfg, spikedata)

% [stats, legend] = spikestatsOverTime(cfg)
% Compute stats over time, trial by trial, for one unit. 
%
% ### Necessary input
% spikedata                  = spike data epoched in FieldTrip trial data structure
% cfg.statstime.spikechannel = label of the unit to analyse
% cfg.statstime.cutlength   = Time window for computing stats. If
%                       cfg.normtime = 'yes', must be between 0 and 1.
%                       Otherwise, it is in seconds. Ignored if cfg.plot =
%                       'trialavg' or 'scatter', because in those cases
%                       only one value is computed per trial.
% cfg.statstime.method      = can be 'isi', 'freq', 'cv2', 'cv', 'fanofactor', 
%                       'burstindex' or 'amplitude'
% cfg.statstime.timelock    = 'begin', 'end', 'no', whether to align data from the
%                       begin or the end of each trial. If 'no', no
%                       alignment is performed and real time is used in
%                       x axis.
% cfg.statstime.plot        = 'raw' (plot all data), 'raw+avg' (plot all data,
%                       and avg), 'avg' (plot only the avg of all the
%                       trials), 'trialavg' (plot avg of each trial),
%                       'scatter' (plot one point per trial), 'movmean',
%                       'movmean+avg'.
%
% ### Optional cfg fields
% cfg.statstime.trial_list = list of trials to analyse. Can be an array of
%                       integers with trials numbers, 'last', or 'all'. 
%                       Default = 'all'.
% cfg.statstime.removebursts= 'yes' or 'no', whether to detect bursts (several
%                       spikes with ISI <10ms), and to keep only the first 
%                       spike. Default = 'no'.
% cfg.statstime.removeempty= 'yes' or 'no', whether to ignore trials with no
%                       spike. If no, and trial has no spike, ISI is put at 
%                       0 and freq is put at Inf. For other cfg.methods, if 
%                       less than 2 spikes, the value is set to NaN 
%                       regardless to cfg.removeempty. Default = 'no'.
% cfg.statstime.removeoutlier= 'yes' or 'no', whether to remove outliers wit
%                       rmoutliers.m. Default = 'no'.
% cfg.statstime.normtime  = 'yes' or 'no', whether to normalize time or not.
%                       Default = 'no'.
% cfg.statstime.normvalues= Only two options for now : 'no' or 'begin' (first
%                       sample value is set to 1). Default = 'no'.
% cfg.statstime.color     = optional field to change the default colors of the
%                       plot. Can be 'defaults', or a color indications
%                       (letter or array of 3 numbers). Default = 'default'.
% cfg.statstime.saveplot  = 'yes' or 'no', whether respectively to save and 
%                       close the plot or to output it. Default = 'no'.
% 
% # Necessary cfg fields if cfg.saveplot = 'yes' :
% cfg.name            = if cfg.saveplot, name of the analysis (for title of 
%                       plot and of saved file)
% cfg.imagesavedir    = if cfg.saveplot, if cfg.saveplot = 'yes', where to  
%                       save the plot.
% cfg.prefix          = if cfg.saveplot = 'yes', prefix attached to the
%                       name of the saved-image
% 
% ### OUTPUT
% stats               = values, avg and std of computed method.
% legend              = target to add a legend afterwards
%
%

%get defaults cfg parameters
cfg.trial_list          = ft_getopt(cfg, 'trial_list'   , 'all');
cfg.removebursts        = ft_getopt(cfg, 'removebursts' , 'no');
cfg.removeempty         = ft_getopt(cfg, 'removeempty'  , 'no');
cfg.removeoutlier       = ft_getopt(cfg, 'removeoutlier', 'no');
cfg.normtime            = ft_getopt(cfg, 'normtime'     , 'no');
cfg.normvalues          = ft_getopt(cfg, 'normvalues'   , 'no');
cfg.color               = ft_getopt(cfg, 'color'        , 'default');
cfg.saveplot            = ft_getopt(cfg, 'saveplot'     , 'no');

if strcmp(cfg.trial_list, 'all')
    cfg.trial_list = 1:size(spikedata.trialinfo,1);
elseif strcmp(cfg.trial_list, 'last')
    cfg.trial_list = size(spikedata.trialinfo,1);
end

%Find unit index
unit_idx = find(strcmp(spikedata.label,cfg.spikechannel));

% Prepare figure
if strcmp(cfg.saveplot, 'yes')
    fig = figure;
end
hold on;

%remove bursts if required
if strcmp(cfg.removebursts, 'yes')
    % get timings and ISIs
    t               = spikedata.time{unit_idx};
    isi_all         = diff(t);
    % counting bursts as in Colder et al. 1996, & Staba et al. 2002
    timeburst       = 0.01; %time in seconds below which isi is considered as burst
    indx            = isi_all < timeburst;
    burstindx       = zeros(size(indx));
    toremove        = [];
    isi_intraburst    = [];
    
    for burstlength = 1 : 10 %nr of consecutive ISI<timeburst
        pattern     = [false, true(1,burstlength), false];
        bindx       = strfind(indx, pattern);
        
        if ~isempty(bindx)
            burstindx(bindx+1) = burstlength; % note +1 because pattern starts with zero
            % add to list to correct for bursts
            for ii = 1 : size(bindx,2)
                % remove all but first spike (at +1)
                toremove = [toremove, bindx(ii)+2:bindx(ii)+2+burstlength-1]; % burstlength = 1; 0 1 0 -> 0 1 x 0. %bindx begins by zero
                % add ISI within bursts
                isi_intraburst = [isi_intraburst, isi_all(bindx(ii)+1:bindx(ii)+burstlength)];
            end
        end
        stats.burstsum(burstlength) = sum(length(bindx));
    end
    % remove subsequenct APs after first AP of a burst
    spikedata.time{unit_idx}(toremove)      = [];
    spikedata.trial{unit_idx}(toremove)     = [];
    spikedata.timestamp{unit_idx}(toremove) = [];
    spikedata.samples{unit_idx}(toremove)   = [];
    spikedata.amplitude{unit_idx}(toremove) = [];
end

%compute isi
spikeisi = ft_spike_isi([], spikedata);

% Initialize values array
if strcmp(cfg.plot, 'trialavg')  || strcmp(cfg.plot, 'scatter')
    max_cuts = 1;
    cfg.cutlength = 0; %ignore cfg.cutlenght because no need to cut data
else
    if strcmp(cfg.normtime, 'yes')
        max_cuts = fix(1/cfg.cutlength)+1;
        if cfg.cutlength > 1, error('With cfg.normtime = ''yes'', cfg.cutlength must be <= 1.'); end
    else
        triallength_max = max(spikedata.trialtime(:,2) - spikedata.trialtime(:,1));
        max_cuts = fix(triallength_max/cfg.cutlength)+1; %+1 because last cut is an incmplete cut
        if cfg.cutlength > triallength_max, error('cfg.cutlength is longer than the max trial length'); end
    end
end
values = nan(size(spikedata.trialinfo,1),max_cuts);

%compute values
for itrial = cfg.trial_list
    clear x freq cv2 cv fanofactor burstindex isi_smooth amplitude
    spike_index = [];
    spike_index = spikedata.trial{unit_idx} == spikedata.trialinfo(itrial,7);
    
    if strcmp(cfg.removeempty, 'no') || sum(spike_index)>1 %if there is more than 1 spike in the trial (otherwise trial is ignored)
        
        %find spike time and isi of the trial
        spike_isi_trial               = spikeisi.isi{unit_idx}(spike_index);
        spike_time_trial              = spikedata.time{unit_idx}(spike_index); %time of the ISI in the trial
        if sum(spike_index) == 0 || sum(spike_index) == 1
            spike_isi_trial           = Inf;
        end
        
        %cut trials if required
        trialbegin              = spikedata.trialtime(itrial,1);
        trialend                = spikedata.trialtime(itrial,2);
        triallength             = trialend - trialbegin;
        if strcmp(cfg.plot, 'trialavg')  || strcmp(cfg.plot, 'scatter')
            ncuts = 0;
            lastcutlength = triallength;
        else
            ncuts                   = fix(triallength /cfg.cutlength);
            lastcutlength           = (triallength /cfg.cutlength - ncuts) * cfg.cutlength; %last uncomplete cut
        end
        
        %normalize time if required
        if strcmp(cfg.normtime,'yes')
            spike_time_trial        = (spike_time_trial - trialbegin) ./ triallength;
            triallength             = 1;
            ncuts                   = fix(1/cfg.cutlength);
            lastcutlength           = (1 /cfg.cutlength - ncuts) * cfg.cutlength; %last uncomplete cut
        end
                       
        %average values for each cut
        for icut=1:ncuts+1 %+1 : last uncomplete cut
            
            spike_cut_index                            = [];
            if icut == ncuts +1 %last cut
                spike_cut_index                            = (spike_time_trial > cfg.cutlength*(icut-1) & spike_time_trial < cfg.cutlength*(icut-1)+lastcutlength);
            else
                spike_cut_index                            = (spike_time_trial > cfg.cutlength*(icut-1) & spike_time_trial < cfg.cutlength*(icut));
            end
            
            if strcmp(cfg.removeempty, 'no') || sum(spike_cut_index) > 1
                
                isi_temp                = [];
                isi_temp                = spike_isi_trial(spike_cut_index);
                isi_smooth(icut)        = nanmean(isi_temp);
                
                if strcmp(cfg.removeempty, 'no') && ismember(sum(spike_cut_index), [0, 1])
                    isi_smooth(icut)    = Inf;
                end
                
                %freq
                freq(icut)              = 1 / isi_smooth(icut);
                
                %amplitude
                amplitude(icut)         = nanmean(spikedata.amplitude{unit_idx}(spike_cut_index));
                
                %CV2, CV, fanofactor
                if sum(spike_cut_index)>3
                    
                    cv(icut)            = nanstd(isi_temp) / nanmean(isi_temp);
                    fanofactor(icut)    = nanstd(isi_temp)^2 / nanmean(isi_temp);
                    burstindex(icut)    = sum(isi_temp<0.010)/sum(isi_temp<0.100);
                    
                    cv2_temp = [];
                    for i = 1:length(isi_temp)-1
                        cv2_temp(i) = 2*abs(isi_temp(i)-isi_temp(i+1))/(isi_temp(i)+isi_temp(i+1));
                    end
                    cv2(icut)           = nanmean(cv2_temp);
                    
                else
                    cv2(icut)               = NaN;
                    cv(icut)                = NaN;
                    fanofactor(icut)        = NaN;
                    burstindex(icut)        = NaN;
                end
                
            else
                isismooth(icut)     = NaN;
                freq(icut)          = NaN;
                cv2(icut)           = NaN;
                cv(icut)            = NaN;
                fanofactor(icut)    = NaN;
                burstindex(icut)    = NaN;
                amplitude(icut)     = NaN;
            end
        end %icut
        
        % select values according ton input method
        switch cfg.method
            case 'isi'
                values(itrial,1:ncuts+1) = isismooth; %+1 because of last uncomplete cut
            case 'freq'
                values(itrial,1:ncuts+1) = freq;
            case 'cv2'
                values(itrial,1:ncuts+1) = cv2;
            case 'cv'
                values(itrial,1:ncuts+1) = cv;
            case 'fanofactor'
                values(itrial,1:ncuts+1) = fanofactor;
            case 'burstindex'
                values(itrial,1:ncuts+1) = burstindex;
            case 'amplitude'
                values(itrial,1:ncuts+1) = amplitude;
            otherwise
                error('cfg.method ''%s'' is not supported by this function', cfg.method);
        end
    end %sum(spike_index)>1
    
end %itrial

% remove outliers if required
if strcmp(cfg.removeoutlier, 'yes')
    values(isoutlier(values)) = nan;
end

% normalize values if required. Other normalizations method may be usefull
% to add later
switch cfg.normvalues
    case 'no'
        %nothing is done
    case 'begin'
        for irow = 1:size(values,1)
            firstsample = find(~isnan(values(irow,:)),1,'first');
            if ~isempty(firstsample)
                values(irow,:) = values(irow,:) ./ values(irow, firstsample);
            end
        end
    otherwise
        error('cfg.normvalues %s is not supported by this function', cfg.normvalues);
end

% align values to the end if required
if strcmp(cfg.timelock, 'end')
    for itrial = 1:size(values,1)
        if sum(~isnan(values(itrial,:))) > 1
            last_idx(itrial)                                                        = find(~isnan(values(itrial,:)),1,'last');
            values_temp(itrial,1:size(values,2)-last_idx(itrial))                   = nan(1,size(values,2)-last_idx(itrial));
            values_temp(itrial,size(values,2)-last_idx(itrial)+1:size(values,2))    = values(itrial,1:last_idx(itrial));
        else
            last_idx(itrial)      = nan;
            values_temp(itrial,:) = nan(1,size(values,2));
        end
    end
    values = values_temp;
    clear values_temp
end

%PLOT


switch cfg.timelock
    
    case 'no'
        %plot non-aligned values, abscisse is real time
        for itrial = 1:size(values,1)
            starttime = spikedata.trialinfo(itrial, 3) / spikedata.hdr.Fs;
            x = (0 : max_cuts-1) * cfg.cutlength + starttime; %-1 because starts at zero
            if strcmp(cfg.color, 'default')
                color = 'k';
            else
                color = cfg.color;
            end
            
            switch cfg.plot
                case 'raw'
                    legend_out = plot(x, values(itrial,:), 'Color',color, 'LineWidth',2);
                case 'movmean'
                    legend_out = plot(x, movmean(values(itrial,:),4), 'Color',color, 'LineWidth',2);
                case 'scatter'
                    legend_out = scatter(starttime, nanmean(values(itrial,:)), 7,'o','filled', 'MarkerEdgeColor',color, 'MarkerFaceColor', color);
                case 'trialavg'
                    endtime = spikedata.trialinfo(itrial, 4) / spikedata.hdr.Fs;
                    legend_out = plot([starttime endtime], [nanmean(values(itrial,:)) nanmean(values(itrial,:))], 'Color',color, 'LineWidth',2);
                case {'raw+avg', 'movmean+avg'}
                    error('with cfg.timelock = ''no'', the option cfg.plot = ''...+avg'' is not availabale');
                case 'avg'
                    error('with cfg.timelock = ''no'', the option cfg.plot = ''avg'' is not availabale');
                otherwise
                    error('cfg.plot = %s is not supported', cfg.plot);
            end
        end
        
    case {'begin', 'end'}
        %plot aligned values, with darker grey for the last
        x = (0 : max_cuts-1) * cfg.cutlength; %-1 because starts at zero
        if strcmp(cfg.timelock, 'end'), x = x-x(end); end %end abscisse at zero if required
        
        c_greys = 0.9 : -0.9/size(values,1) : 0; %Color darker for the last trials
        
        for itrial = 1:size(values,1)
            if strcmp(cfg.color, 'default')
                color = [c_greys(itrial), c_greys(itrial), c_greys(itrial)];
            else
                color = cfg.color;
            end
            
            switch cfg.plot
                case {'raw', 'raw+avg'}
                    legend_out = plot(x, values(itrial,:), 'Color', color);
                case {'movmean', 'movmean+avg'}
                    legend_out = plot(x, movmean(values(itrial,:),4), 'Color', color);
                case 'avg'
                    %do nothing here, avg is plotted later
                case 'trialavg'
                    no_nan_idx = find(~isnan(values(itrial,:)));
                    if ~isempty(no_nan_idx)
                        legend_out = plot([x(no_nan_idx(1)) x(no_nan_idx(end))], [nanmean(values(itrial,:)) nanmean(values(itrial,:))], 'Color', color);
                    end
                case 'scatter'
                    legend_out = scatter(0, nanmean(values(itrial,:)), 7,'o','filled', 'MarkerEdgeColor',color, 'MarkerFaceColor', color);
                otherwise
                    error('cfg.plot = %s is not supported', cfg.plot);
            end
        end
    otherwise
        error('cfg.timelock = %s is not supported', cfg.timelock);
end

%plot avg, if required
values_avg = nanmean(values,1);
values_std = nanstd(values,0,1);%not plot, compute for output
if strcmp(cfg.plot, 'raw+avg') || strcmp(cfg.plot, 'avg') || strcmp(cfg.plot, 'movmean+avg')
    avg = plot(x,values_avg,'g','LineWidth',2);
    if strcmp(cfg.plot, 'avg'), legend_out = avg; end
end

%plot non-timelocked extremity if necessary
if ~strcmp(cfg.normtime, 'yes') && ~strcmp(cfg.timelock, 'no') && (strcmp(cfg.plot, 'raw') || strcmp(cfg.plot, 'raw+avg'))
    if strcmp(cfg.timelock, 'end')
        for itrial = 1:size(values,1)
            if last_idx(itrial) > 0
                first_idx = find(~isnan(values(itrial,:)),3);
                if length(first_idx) > 1
                    %correct bug in plot : there is no point if it is a single number followed by a nan.
                    try
                        if first_idx(2) == first_idx(1) +1,     first_idx = first_idx(1);
                        elseif first_idx(3) == first_idx(2) +1, first_idx = first_idx(2);
                        else,                                   first_idx = first_idx(3);
                        end
                        
                        scatter(x(first_idx),values(itrial, first_idx), '.r');
                    end
                end
            end
        end
    elseif strcmp(cfg.timelock, 'begin')
        for itrial = 1:size(values,1)
            clear idx
            idx = find(~isnan(values(itrial,:)));
            if ~isempty(idx)
                scatter(x(idx(end)),values(itrial, idx(end)), '.r');
            end
        end
    end
end

axis tight

% if strcmp(cfg.normtime, 'yes'), norm = ' (normalized)'; else,norm = []; end
% xlabel(sprintf('Time from begining of trial%s (s)',norm));
ylabel(cfg.method);
% title(sprintf('%s, %s : %d trials', cfg.name, spikedata.label{unit_idx},size(spikedata.trialinfo,1)), 'Fontsize',18, 'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');

%% Output stats
stats.cfg       = cfg;
stats.values    = values;
stats.time      = x;
stats.avg       = values_avg;
stats.std       = values_std;


%% save plot
if strcmp(cfg.saveplot, 'yes')
    set(gca,'Fontsize',15);
    
    if ~(exist (cfg.imagesavedir)==7)
        mkdir(cfg.imagesavedir);
        fprintf('Create forlder %s',cfg.imagesavedir);
    end
    
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    
    if strcmp(cfg.normtime, 'yes'), norm = '_timenorm'; else,norm = []; end
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'_',cfg.name,'_',spikedata.label{unit_idx},'_spikestatstrial_',cfg.method,norm,'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'_',cfg.name,'_',spikedata.label{unit_idx},'_spikestatstrial_',cfg.method,norm,'.png']),'-r600');
    close all
    
end

end

