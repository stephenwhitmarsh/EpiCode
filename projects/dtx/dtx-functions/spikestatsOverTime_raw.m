function [stats] = spikestatsOverTime_raw(cfg, spikedata, hdr)

% spikedata = raw spike data, not epoched into trials
% stats : one field per unit
% stats{ipart}.time : data relative to the raw spikedata (not to the latency)

% après : remove artefacts, plot pour chaque unit

%get defaults cfg parameters
cfg.statstime                     = ft_getopt(cfg, 'statstime', []);
cfg.statstime.latency             = ft_getopt(cfg.statstime, 'latency'      , 'all');
cfg.statstime.timewin             = ft_getopt(cfg.statstime, 'timewin'      , 10);
cfg.statstime.slidestep           = ft_getopt(cfg.statstime, 'slidestep'    , 1);
cfg.statstime.removebursts        = ft_getopt(cfg.statstime, 'removebursts' , 'no');
cfg.statstime.removeempty         = ft_getopt(cfg.statstime, 'removeempty'  , 'no');

%Find latency
if strcmp(cfg.statstime.latency, 'all')
    cfg.statstime.latency = [0-hdr.nSamplesPre hdr.nSamples] ./ hdr.Fs;
end


for ipart = size(spikedata, 2)
    
    stats{ipart}.cfg = cfg;
    stats{ipart}.hdr = hdr;
    stats{ipart}.label = spikedata{ipart}.label;
   
    %convert raw spikedata{ipart} into 1-trial spikedata{ipart}.
    cfgtemp                     = [];
    cfgtemp.trl                 = [cfg.statstime.latency, 0, cfg.statstime.latency] .* hdr.Fs;
    cfgtemp.trlunit             = 'samples';
    cfgtemp.hdr                 = hdr;
    cfgtemp.timestampspersecond = hdr.TimeStampPerSample * hdr.Fs;
    spikedata{ipart}            = ft_spike_maketrials(cfgtemp, spikedata{ipart});
    
    %remove bursts if required
    if strcmp(cfg.statstime.removebursts, 'yes')
        for i_unit = 1:size(spikedata{ipart}.label, 2)
            % get timings and ISIs
            t               = spikedata{ipart}.time{i_unit};
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
                stats{ipart}.burstsum{i_unit}(burstlength) = sum(length(bindx));
            end
            % remove subsequenct APs after first AP of a burst
            spikedata{ipart}.time{i_unit}(toremove)      = [];
            spikedata{ipart}.trial{i_unit}(toremove)     = [];
            spikedata{ipart}.timestamp{i_unit}(toremove) = [];
            spikedata{ipart}.samples{i_unit}(toremove)   = [];
            spikedata{ipart}.amplitude{i_unit}(toremove) = [];
        end
    end
    
    %compute isi
    spikeisi = ft_spike_isi([], spikedata{ipart});
    
    % Compute stats on sliding time window along all data.
    % the x values are relative to the input spikedata{ipart}, not to the latency, so
    % we keep track of where is each measurement on the original data
    for i_unit = 1:size(spikedata{ipart}.label, 2)
        
        ft_progress('init', 'text',     sprintf('Unit %d/%d', i_unit, size(spikedata{ipart}.label, 2)));
        
        %initialize values for each unit
        t_start     = 0;
        t_end       = cfg.statstime.timewin;
        i_window    = 0;
        data_length = spikedata{ipart}.trialtime(2) - spikedata{ipart}.trialtime(1);
        nr_windows  = length(cfg.statstime.timewin : cfg.statstime.slidestep : data_length);
        
        %go trhough each window
        for t_end = cfg.statstime.timewin : cfg.statstime.slidestep : data_length
            
            % compute t_start time and display wait bar
            t_start = t_end - cfg.statstime.timewin;
            i_window = i_window+1;
            ft_progress(i_window/nr_windows, 'Computing stats on window %d from %d', i_window, nr_windows);
            
            %find spikes in this window
            spike_index = [];
            spike_index = spikedata{ipart}.time{i_unit} > t_start & spikedata{ipart}.time{i_unit} <=  t_end;
            
            %store time for later analysis
            stats{ipart}.window_times{i_unit}(i_window,1) = t_start;
            stats{ipart}.window_times{i_unit}(i_window,2) = t_end;
            stats{ipart}.time{i_unit}(i_window) = t_end + spikedata{ipart}.trialtime(1); %recover original time

            
            % if no or one spike in this window, set all to NaN (except freq and ISI if removempty == 'no'
            if strcmp(cfg.statstime.removeempty, 'yes') && sum(spike_index)<=1
                stats{ipart}.isismooth{i_unit}(i_window)     = NaN;
                stats{ipart}.freq{i_unit}(i_window)          = Nan;
                stats{ipart}.cv2{i_unit}(i_window)           = NaN;
                stats{ipart}.cv{i_unit}(i_window)            = NaN;
                stats{ipart}.fanofactor{i_unit}(i_window)    = NaN;
                stats{ipart}.burstindex{i_unit}(i_window)    = NaN;
                stats{ipart}.amplitude{i_unit}(i_window)     = NaN;
                continue
            end
            
            %find spike time and isi of the window
            spike_isi_win               = spikeisi.isi{i_unit}(spike_index);
            spike_time_win              = spikedata{ipart}.time{i_unit}(spike_index); %time of the ISI in the trial
            if sum(spike_index) <= 1
                spike_isi_win           = Inf;
            end
            
            %compute average value for each window
            stats{ipart}.isi_smooth{i_unit}(i_window)        = nanmean(spike_isi_win);
            stats{ipart}.freq{i_unit}(i_window)              = 1 / nanmean(spike_isi_win);
            stats{ipart}.amplitude{i_unit}(i_window)         = nanmean(spikedata{ipart}.amplitude{i_unit}(spike_index));
            
            %CV2, CV, fanofactor, burstindex : do not compute value if less than 3 spikes
            if sum(spike_index)<=3
                stats{ipart}.cv2{i_unit}(i_window)           = NaN;
                stats{ipart}.cv{i_unit}(i_window)            = NaN;
                stats{ipart}.fanofactor{i_unit}(i_window)    = NaN;
                stats{ipart}.burstindex{i_unit}(i_window)    = NaN;
                continue
            end
            
            stats{ipart}.cv{i_unit}(i_window)            = nanstd(spike_isi_win) / nanmean(spike_isi_win);
            stats{ipart}.fanofactor{i_unit}(i_window)    = nanstd(spike_isi_win)^2 / nanmean(spike_isi_win);
            stats{ipart}.burstindex{i_unit}(i_window)    = sum(spike_isi_win<0.010)/sum(spike_isi_win<0.100);
            
            cv2_temp = [];
            for i = 1:length(spike_isi_win)-1
                cv2_temp(i) = 2*abs(spike_isi_win(i)-spike_isi_win(i+1))/(spike_isi_win(i)+spike_isi_win(i+1));
            end
            stats{ipart}.cv2{i_unit}(i_window)          = nanmean(cv2_temp);
            
        end %t_end
        ft_progress('close');
    end %i_unit
    
end %ipart

end

