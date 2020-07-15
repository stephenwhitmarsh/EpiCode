function [stats] = spikestatsOverTime(cfg, spikedata,force)

% spikedata = fieldtrip spike data, epoched into trials (can be only 1
% trial)
% outpout : stats{ipart}{ilabel}.(param){i_unit}{i_trial} : one value per window
% time = time of the end of the window

% note : if time window is too big compared to trial length, return average
% value for the whole trial 

%get defaults cfg parameters
cfg.statstime                     = ft_getopt(cfg, 'statstime', []);
cfg.statstime.timewin             = ft_getopt(cfg.statstime, 'timewin'      , 10);
cfg.statstime.slidestep           = ft_getopt(cfg.statstime, 'slidestep'    , 1);
cfg.statstime.removebursts        = ft_getopt(cfg.statstime, 'removebursts' , 'no');
cfg.statstime.removeempty         = ft_getopt(cfg.statstime, 'removeempty'  , 'no');
cfg.statstime.suffix              = ft_getopt(cfg.statstime, 'suffix'       , []);
cfg.statstime.label_list          = ft_getopt(cfg.statstime, 'label_list'   , 'all');
cfg.statstime.write               = ft_getopt(cfg.statstime, 'write'        , false);

% load precomputed stats if required
fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikestatsOverTime', cfg.statstime.suffix, '.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed spike stats over time\n');
    load(fname,'stats');
    return
end

for ipart = size(spikedata, 2)
    
    if strcmp(cfg.statstime.label_list, 'all')
        cfg.statstime.label_list = 1:size(spikedata{ipart},2);
    end
    
    for ilabel = cfg.statstime.label_list
        
        stats{ipart}{ilabel}.cfg = cfg;
        stats{ipart}{ilabel}.label = spikedata{ipart}{ilabel}.label;
        stats{ipart}{ilabel}.analysis_name = cfg.name{ilabel};
        
        %remove bursts if required
        if strcmp(cfg.statstime.removebursts, 'yes')
            for i_unit = 1:size(spikedata{ipart}{ilabel}.label, 2)
                % get timings and ISIs
                t               = spikedata{ipart}{ilabel}.time{i_unit};
                isi_all         = diff(t);            % counting bursts as in Colder et al. 1996, & Staba et al. 2002
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
                    stats{ipart}{ilabel}.burstsum{i_unit}(burstlength) = sum(length(bindx));
                end
                % remove subsequenct APs after first AP of a burst
                spikedata{ipart}{ilabel}.time{i_unit}(toremove)      = [];
                spikedata{ipart}{ilabel}.trial{i_unit}(toremove)     = [];
                spikedata{ipart}{ilabel}.timestamp{i_unit}(toremove) = [];
                spikedata{ipart}{ilabel}.samples{i_unit}(toremove)   = [];
                spikedata{ipart}{ilabel}.amplitude{i_unit}(toremove) = [];
            end
        end
        
        %compute isi
        spikeisi = ft_spike_isi([], spikedata{ipart}{ilabel});
        
        % Compute stats on sliding time window along all data.
        for i_unit = 1:size(spikedata{ipart}{ilabel}.label, 2)
            ft_progress('init', 'text',     sprintf('Part %d/%d, Label %d/%d, Unit %d/%d (%d trials):', ipart,size(spikedata, 2),ilabel, size(spikedata{ipart}, 2), i_unit, size(spikedata{ipart}{ilabel}.label, 2), size(spikedata{ipart}{ilabel}.trialinfo, 1)));
            trials_win_count = 0;
            
            for itrial = 1:size(spikedata{ipart}{ilabel}.trialinfo, 1)
                %initialize values for each unit and each trial
                %Startsample and EndSample in trialinfo : columns 3 and 4.
                starttrial      = spikedata{ipart}{ilabel}.trialinfo(itrial, 3) / spikedata{ipart}{ilabel}.hdr.Fs;
                endtrial        = spikedata{ipart}{ilabel}.trialinfo(itrial, 4) / spikedata{ipart}{ilabel}.hdr.Fs;
                trial_length    = endtrial - starttrial;
                i_window        = 0;
                nr_windows      = length(cfg.statstime.timewin : cfg.statstime.slidestep : trial_length);
                nr_win_all_trials = trials_win_count + nr_windows * (size(spikedata{ipart}{ilabel}.trialinfo, 1) - itrial);
                    
                %if time window is too big compared to trial length, return
                %average value for the whole trial
                if nr_windows == 0
                    trial_length = cfg.statstime.timewin;
                end
                
                %go trhough each window
                for t_end = cfg.statstime.timewin : cfg.statstime.slidestep : trial_length
                    
                    % compute window time and display wait bar
                    t_start = t_end - cfg.statstime.timewin;
                    ts = t_start + spikedata{ipart}{ilabel}.trialtime(itrial, 1);
                    te = t_end + spikedata{ipart}{ilabel}.trialtime(itrial, 1);
                    i_window = i_window+1;
                    trials_win_count = trials_win_count + 1;
                    ft_progress(trials_win_count/(nr_win_all_trials+i_window));
                    
                    %find spikes in this window
                    spike_index = (spikedata{ipart}{ilabel}.trial{i_unit}==itrial) & (spikedata{ipart}{ilabel}.time{i_unit} > ts) & (spikedata{ipart}{ilabel}.time{i_unit} <=  te);
                    
                    %store time for later analysis
                    stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,1) = t_start+starttrial;
                    stats{ipart}{ilabel}.window_times{i_unit}{itrial}(i_window,2) = t_end+starttrial;
                    stats{ipart}{ilabel}.time{i_unit}{itrial}(i_window) = t_end + starttrial; %recover original time
                    
                    % if no or one spike in this window, set all to NaN (except freq and ISI if removempty == 'no')
                    if strcmp(cfg.statstime.removeempty, 'yes') && sum(spike_index)<=1
                        stats{ipart}{ilabel}.isismooth{i_unit}{itrial}(i_window)     = NaN;
                        stats{ipart}{ilabel}.freq{i_unit}{itrial}(i_window)          = NaN;
                        stats{ipart}{ilabel}.cv2{i_unit}{itrial}(i_window)           = NaN;
                        stats{ipart}{ilabel}.cv{i_unit}{itrial}(i_window)            = NaN;
                        stats{ipart}{ilabel}.fanofactor{i_unit}{itrial}(i_window)    = NaN;
                        stats{ipart}{ilabel}.burstindex{i_unit}{itrial}(i_window)    = NaN;
                        stats{ipart}{ilabel}.amplitude{i_unit}{itrial}(i_window)     = NaN;
                        continue
                    end
                    
                    %find spike time and isi of the window
                    spike_isi_win               = spikeisi.isi{i_unit}(spike_index);
                    if sum(spike_index) <= 1
                        spike_isi_win           = Inf;
                    end
                    
                    %compute average value for each window
                    stats{ipart}{ilabel}.isi_smooth{i_unit}{itrial}(i_window)        = nanmean(spike_isi_win);
                    stats{ipart}{ilabel}.freq{i_unit}{itrial}(i_window)              = 1 / nanmean(spike_isi_win);
                    stats{ipart}{ilabel}.amplitude{i_unit}{itrial}(i_window)         = nanmean(spikedata{ipart}{ilabel}.amplitude{i_unit}(spike_index));
                    
                    %CV2, CV, fanofactor, burstindex : do not compute value if less than 3 spikes
                    if sum(spike_index)<=3
                        stats{ipart}{ilabel}.cv2{i_unit}{itrial}(i_window)           = NaN;
                        stats{ipart}{ilabel}.cv{i_unit}{itrial}(i_window)            = NaN;
                        stats{ipart}{ilabel}.fanofactor{i_unit}{itrial}(i_window)    = NaN;
                        stats{ipart}{ilabel}.burstindex{i_unit}{itrial}(i_window)    = NaN;
                        continue
                    end
                    
                    stats{ipart}{ilabel}.cv{i_unit}{itrial}(i_window)            = nanstd(spike_isi_win) / nanmean(spike_isi_win);
                    stats{ipart}{ilabel}.fanofactor{i_unit}{itrial}(i_window)    = nanstd(spike_isi_win)^2 / nanmean(spike_isi_win);
                    stats{ipart}{ilabel}.burstindex{i_unit}{itrial}(i_window)    = sum(spike_isi_win<0.010)/sum(spike_isi_win<0.100);
                    
                    cv2_temp = [];
                    for i = 1:length(spike_isi_win)-1
                        cv2_temp(i) = 2*abs(spike_isi_win(i)-spike_isi_win(i+1))/(spike_isi_win(i)+spike_isi_win(i+1));
                    end
                    stats{ipart}{ilabel}.cv2{i_unit}{itrial}(i_window)          = nanmean(cv2_temp);
                    
                end %t_end
            end %itrial
            ft_progress('close');
        end %i_unit
    end %ilabel
end %ipart

save(fname, 'stats', '-v7.3');

end

