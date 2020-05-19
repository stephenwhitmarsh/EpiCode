function stats = spikestatsOverTime(cfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats = spikestatsOverTime(cfg)
% Compute stats over time, trial by trial, for one unit. Time can be
% normalized or not. The computed stat is selected in cfg.
% 
% ### INPUT
% cfg.spikedata       = spike data epoched in FieldTrip trial data structure
% cfg.spikeisi        = output of ft_spike_isi applied to spikedata
% cfg.spikechannel    = label of the unit to analyse
% cfg.normtime        = 'yes' or 'no', whether to normalize time or not
% cfg.cutlength       = Time window for computing stats on ISI. If 
%                       normtime = 'yes', must be between 0 and 1. 
%                       Otherwise, it is in seconds.
% cfg.method          = can be 'isi', 'freq', 'cv2', 'cv', 'fanofactor', or
%                       'burstindex'
% cfg.removeempty     = 'yes' or 'no', whether to ignore trials with no
%                       spike (for freq or ISI, otherwise need at least 2
%                       spikes)
% cfg.removeoutlier   = 'yes' or 'no', whether to remove outliers wit
%                       rmoutliers.m
% cfg.saveplot        = 'yes' or 'no', whether to save and close the
%                       plot or to output it (respectively)
% cfg.name            = name of the analysis (for title of plot and of
%                       saved file)
% cfg.imagesavedir    = if saveplot = 'yes', where to save the plot.
% cfg.prefix          = if saveplot = 'yes', prefix attached to the name of
%                       the saved-image
% 
% ### OUTPUT
% stats               = values, avg and std of computed method.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize values array
if strcmp(cfg.normtime, 'yes')
    max_cuts = fix(1/cfg.cutlength);
    values = nan(size(cfg.spikedata.trialinfo,1),max_cuts);
    if cfg.cutlength > 1, error('With cfg.normtime = ''yes'', cfg.cutlength must be <= 1.'); end
else
    triallength_max = max(cfg.spikedata.trialtime(:,2) - cfg.spikedata.trialtime(:,1));
    max_cuts = fix(triallength_max/cfg.cutlength);
    values = nan(size(cfg.spikedata.trialinfo,1),max_cuts);
    if cfg.cutlength > triallength_max, error('cfg.cutlength is longer than the max trial length'); end
end

%Find unit index
unit_idx = find(strcmp(cfg.spikedata.label,cfg.spikechannel));

% Prepare figure
if strcmp(cfg.saveplot, 'yes')
    fig = figure;
end
hold on;


for itrial = 1:size(cfg.spikedata.trialinfo,1)
    clear x freq cv2 cv fanofactor burstindex isi_smooth 
    spike_index = [];
    spike_index = cfg.spikedata.trial{unit_idx} == cfg.spikedata.trialinfo(itrial,7);
    
    if strcmp(cfg.removeempty, 'no') || sum(spike_index)>1 %if there is more than 1 spike in the trial (otherwise trial is ignored)
        
        %find spike and isi of the trial, trial length, and nr of cuts to do
        spike_isi               = cfg.spikeisi.isi{unit_idx}(spike_index);
        spike_time              = cfg.spikedata.time{unit_idx}(spike_index); %time of the ISI in the trial
        if sum(spike_index) == 0 || sum(spike_index) == 1
            spike_isi           = Inf;
        else
            spike_isi               = spike_isi(2:end); %remove first NaN
            spike_time              = spike_time(2:end); %remove first NaN
        end
        trialbegin              = cfg.spikedata.trialtime(itrial,1);
        trialend                = cfg.spikedata.trialtime(itrial,2);
        triallength             = trialend - trialbegin;
        ncuts                   = fix(triallength /cfg.cutlength);
        
        %normalize time if required
        if strcmp(cfg.normtime,'yes')
            spike_time              = (spike_time - trialbegin) ./ triallength;
            triallength             = 1;
            ncuts                   = fix(1/cfg.cutlength);
        end
        
        %average values for each cut
        for icut=1:ncuts 
            
            isi_smooth_index                            = [];
            isi_smooth_index                            = (spike_time > cfg.cutlength*(icut-1) & spike_time < cfg.cutlength*(icut));
            
            if strcmp(cfg.removeempty, 'no') || sum(isi_smooth_index) >0  
                isi_temp                = [];
                isi_temp                = spike_isi(isi_smooth_index);
                isi_smooth(icut)        = nanmean(isi_temp);
                
                if strcmp(cfg.removeempty, 'no') && sum(isi_smooth_index) ==0  
                    isi_smooth(icut)    = Inf;
                end
                
                %abscisse value for plot
                x(icut)                 = cfg.cutlength*icut;
                
                %freq
                freq(icut)              = 1 / isi_smooth(icut);
                
                %CV2, CV, fanofactor
                if sum(isi_smooth_index)>2
                    
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
            end
        end
        
        % select values according ton input method
        switch cfg.method
            case 'isi'
                values(itrial,1:ncuts) = isismooth;
            case 'freq'
                values(itrial,1:ncuts) = freq;
            case 'cv2'
                values(itrial,1:ncuts) = cv2;
            case 'cv'
                values(itrial,1:ncuts) = cv;
            case 'fanofactor'
                values(itrial,1:ncuts) = fanofactor;
            case 'burstindex'
                values(itrial,1:ncuts) = burstindex;
            otherwise
                error('cfg.method %s is not supported by this function', cfg.method);
        end
        
    end %sum(spike_index)>1

end %itrial

if strcmp(cfg.removeoutlier, 'yes')
    values = rmoutliers(values);
end

% plot values
x = cfg.cutlength : cfg.cutlength : max_cuts * cfg.cutlength;
plot(x, values);

%plot avg +/- SD
values_avg = nanmean(values,1);
values_std = nanstd(values,0,1);

x = cfg.cutlength : cfg.cutlength : cfg.cutlength*size(values,2);

avg = plot(x,values_avg,'k','LineWidth',2);

% y = [values_avg - values_std; values_std; values_std]';
% filled_SD = area(x,y);%,'FaceAlpha',[0, 0.5, 0.5]);
% filled_SD(1).FaceAlpha = 0; filled_SD(2).FaceAlpha = 0.3; filled_SD(3).FaceAlpha = 0.3;
% filled_SD(1).EdgeColor = 'none'; filled_SD(2).EdgeColor = 'none'; filled_SD(3).EdgeColor = 'none';
% filled_SD(2).FaceColor = 'g'; filled_SD(3).FaceColor = 'g';
% legend([avg, filled_SD(2)],'Average','SD');

axis tight
ax = axis;
ylim([0, ax(4)]);

if strcmp(cfg.normtime, 'yes'), norm = ' (normalized)'; else,norm = []; end
% xlabel(sprintf('Time from begining of trial%s (s)',norm));
ylabel(cfg.method);
% title(sprintf('%s, %s : %d trials', cfg.name, cfg.spikedata.label{unit_idx},size(cfg.spikedata.trialinfo,1)), 'Fontsize',18, 'Interpreter','none');
set(gca,'FontWeight','bold' );
set(gca,'TickDir','out');

%% Output stats
stats.method    = cfg.method;
stats.value     = values;
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
    print(fig, '-dpdf', fullfile(cfg.imagesavedir,[cfg.prefix,'_',cfg.name,'_',cfg.spikedata.label{unit_idx},'_spikestatstrial_',cfg.method,norm,'.pdf']),'-r600');
    print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'_',cfg.name,'_',cfg.spikedata.label{unit_idx},'_spikestatstrial_',cfg.method,norm,'.png']),'-r600');
    close all

end

end

