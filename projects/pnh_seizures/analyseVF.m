function analyseVF(MuseStruct)

addpath /Users/stephen.whitmarsh/WhitmarshEpilepsy/
addpath /Users/stephen.whitmarsh/fieldtrip/
ft_defaults
 
%% load or create MuseStruct

if isempty(MuseStruct)
    if exist('MuseStruct.mat','file')
        load('MuseStruct.mat','MuseStruct');
    else
        patient_directory       = '/Volumes/epimicro/Donnees-analyses/Stephen/pat_02230_0674/eeg';
        directory_searchstring  = '02230_2015-02-*';
        data_searchstring       = '*m1*.ncs';
        MuseStruct              = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);
        save('MuseStruct','MuseStruct');
    end
end

%% parameters

startend                = { 'VF1__START__','VF1__END__'; ...
    'RR__START__','RR__END__'; 'P','P'; 'PP__START__','PP__END__'};             % Muse Marker names to analyze
label                   = { 'VF1','VF2','P','PP'};                              % Own names 
prestim                 = [1,1,1,2.0];                                          % list of onset timing with respect to start-marker (s)
poststim                = [0,1,1,4.0];                                          % list of offset timing with respect to end-marker (s)
prepeak                 = [0.5,0.8,0.2,1];                                      % list of onset timing to align to (first) peak (s)
postpeak                = [0.8,0.8,0.8,1.5];                                    % list of offset timing to align to (first) peak (s)
baseline                = {'no','no','no',[-1,-0.5]};                           % whether to baseline correct after peak detection, if so, the baseline period (s)
alignmode               = {'first','first','max','first'};                      % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
avgmode                 = {'first','first','all','first'};                      % whether to average only the first peak, or all peaks {'first' or 'all}
filenrs                 = [1,4,6,7,8];                                          % do not bother with artefacted channels or reference channel (nr. 2 is reference in pat_02230)
resamplefs              = 640;                                                  % resample data to speed up and reduce memory (Hz)
stepsize                = [0.2, 0.05, 0.01, 0.01];                              % for histogram of duration of events in overview (s)
lpfreq                  = [4,4,40,40];                                          % lowpass filter freq to smooth peak detection (Hz)
peak_search_toi_ac      = [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];             % active period in which to search for peaks 
peak_search_toi_bl      = [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];             % baseline period in which to search for peaks
datasavedir             = '/Users/stephen.whitmarsh/WhitmarshEpilepsy/data';    % where to write data
imagesavedir            = '/Users/stephen.whitmarsh/WhitmarshEpilepsy/images';  % where to print images

%% Make overview of data as segmented by markers aligned to (first) peak

for imarker = 1 : size(startend,1)
    
    fname = fullfile(datasavedir,['dat_align_',label{imarker},'.mat']);
    
    if exist(fname,'file')
        
        fprintf('Found data in %s, loading... \n',fname);
        load(fname,'dat_align','dat_filt_chanavg_align','locs_ac_sel_trl','locs_ac_sel_avg','t1_ac_indx');
        
    else
        
        fprintf('Aligned data not found - making');
        hasmarker = zeros(length(MuseStruct),1);
        
        for idir = 1:length(MuseStruct)
            
            hdr = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(1).name));
            
            if isfield(MuseStruct{idir},'markers')
                if isfield(MuseStruct{idir}.markers,(startend{imarker,1}))
                    if ~isempty(MuseStruct{idir}.markers.(startend{imarker,1}).events)
                        
                        % create Fieldtrip trl
                        Startsample     = [];
                        Endsample       = [];
                        for ievent = 1 : size(MuseStruct{idir}.markers.(startend{imarker,1}).events,2)
                            Startsample             = [Startsample; MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).sample - prestim(imarker)*hdr.Fs];
                            Endsample               = [Endsample;   MuseStruct{idir}.markers.(startend{imarker,2}).events(ievent).sample + poststim(imarker)*hdr.Fs];
                            Offset                  = -ones(size(Endsample)) * prestim(imarker) * hdr.Fs;
                        end
                        
                        % load trials for each channel
                        for ifile = filenrs
                            cfg = [];
                            cfg.dataset             = fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(ifile).name);
                            cfg.trl                 = [Startsample, Endsample, Offset];
                            cfg.trl                 = cfg.trl(Startsample > 0 & Endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs
                            
                            filedat{ifile}          = ft_preprocessing(cfg);
                            filedat{ifile}.label{1} = filedat{ifile}.label{1}(end-6:end); % truncate label
                        end
                        
                        % concatinate channels
                        dirdat{idir}                = ft_appenddata([],filedat{filenrs});
                        clear filedat
                        
                        % downsample data and baseline correction
                        cfg                     = [];
                        cfg.resamplefs          = resamplefs;
                        if strcmp(baseline(imarker),'no')
                            cfg.demean          = 'no';
                        else
                            cfg.demean          = 'yes';
                            cfg.baselinewindow  = baseline{imarker};
                        end
                        dirdat{idir}                = ft_resampledata(cfg,dirdat{idir} );
                        
                        % flag for averaging
                        hasmarker(idir) = 1;
                    end
                end
            end
        end
        
        % concatinate data over trials
        dat = ft_appenddata([],dirdat{find(hasmarker)});
        clear dirdat
        
        % filtering data for peak detection
        dat_filt = dat;
        for itrial = 1 : size(dat.trial,2)
            for ichan = 1 : size(dat.label,1)
                dat_filt.trial{itrial}(ichan,:) = bandpassFilter(dat.trial{itrial}(ichan,:),resamplefs,1.0,lpfreq(imarker));
            end
        end
        
        % peak threshold detection based on baseline period
        cfg = [];
        cfg.avgoverchan     = 'yes';
        dat_filt_chanavg    = ft_selectdata(cfg,dat_filt); clear dat_filt
        max_peaks_ac_list   = [];
        
        for itrial = 1 : size(dat.trial,2)
            t1_bl_indx(itrial)                 = find(dat_filt_chanavg.time{itrial} > peak_search_toi_bl(imarker,1),1,'first');
            t2_bl_indx(itrial)                 = find(dat_filt_chanavg.time{itrial} < peak_search_toi_bl(imarker,2),1,'last');
            t1_ac_indx(itrial)                 = find(dat_filt_chanavg.time{itrial} > peak_search_toi_ac(imarker,1),1,'first');
            t2_ac_indx(itrial)                 = find(dat_filt_chanavg.time{itrial} < (dat_filt_chanavg.time{itrial}(end)-poststim(imarker)+peak_search_toi_ac(imarker,2)),1,'last'); % with respect to end of trial
            
            [peaks_bl{itrial},~]               = findpeaks(dat_filt_chanavg.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
            max_peaks_bl(itrial)               = mean(peaks_bl{itrial});
            [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_chanavg.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));
            max_peaks_ac(itrial)               = max(peaks_ac{itrial});
            max_peaks_ac_list                  = [max_peaks_ac_list, peaks_ac{itrial}];
            
        end
        
        for itrial = 1 : size(dat.trial,2)
            
            % threshold based on median of max peaks in all baseline periods
            peaks_sel_avg               = peaks_ac{itrial} > mean(max_peaks_bl) * 1;
            peaks_ac_sel_avg{itrial}    = peaks_ac{itrial}(peaks_sel_avg);
            locs_ac_sel_avg{itrial}     = locs_ac{itrial}(peaks_sel_avg);
            
            % threshold based on median of max peaks in trail-by-trial baseline
            peaks_sel_trl               = peaks_ac{itrial} > max(peaks_bl{itrial}) * 1;
            peaks_ac_sel_trl{itrial}    = peaks_ac{itrial}(peaks_sel_trl);
            locs_ac_sel_trl{itrial}     = locs_ac{itrial}(peaks_sel_trl);
            
        end
        
        % find trials that have no peaks detected (within 500 ms)
        haspeak = ones(1,size(dat.trial,2));
        for itrial = 1 : size(dat.trial,2)
            if isempty(locs_ac_sel_avg{itrial})
                haspeak(itrial) = 0;
                fprintf('Could not find peak in trial number %d. Will remove trial...\n',itrial);
            end
        end
        
        % align trials to peak
        dat_align               = dat;
        dat_filt_chanavg_align  = dat_filt_chanavg;
        for itrial = 1 : size(dat.trial,2)
            if ~isempty(locs_ac_sel_avg{itrial}) && haspeak(itrial) == 1
                
                % find relevant peak according to method
                if strcmp(alignmode(imarker),'nearest')
                    [~, ip] = min(abs(dat_filt_chanavg_align.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1))); % peak-time closes to zero
                elseif strcmp(alignmode(imarker),'max')
                    [~, ip] = max(dat_filt_chanavg_align.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1)); % max peak
                elseif strcmp(alignmode(imarker),'first')
                    [~, ip] = find(dat_filt_chanavg_align.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1)>0,1,'first'); % first peak after 0
                end
                
                % align time axis to relevant peak
                dat_align.time{itrial} = dat_align.time{itrial} - dat_align.time{itrial}(locs_ac_sel_avg{itrial}(ip)+t1_ac_indx(itrial)-1);
                dat_filt_chanavg_align.time{itrial} = dat_filt_chanavg_align.time{itrial} - dat_filt_chanavg_align.time{itrial}(locs_ac_sel_avg{itrial}(ip)+t1_ac_indx(itrial)-1);
            end
        end
             
        % select only trials that had a peak detected
        cfg = [];
        cfg.trials = find(haspeak);
        dat_align = ft_selectdata(cfg,dat_align);
        
        % save to disk
        save(fname,'dat_align','dat_filt_chanavg_align','locs_ac_sel_trl','locs_ac_sel_avg','t1_ac_indx');
        clear dat dat_filt_chanavg_align  
    end
     
    % map all trials onto time-independant grid (average over channels)
    mintime                 = 0;
    maxtime                 = 0;
    maxsamples              = 0;
    for itrial = 1 : size(dat_align.trial,2)
        if dat_align.time{itrial}(1) < mintime
            mintime = dat_align.time{itrial}(1);
        end
        if dat_align.time{itrial}(end) > maxtime
            maxtime = dat_align.time{itrial}(end);
        end
        if length(dat_align.time{itrial}) > maxsamples
            maxsamples = length(dat_align.time{itrial});
        end   
    end
    
    time      = mintime : 1/resamplefs : maxtime;
    trialgrid = zeros(size(dat_align.trial,2),length(time));
    
    i = 1;
    for itrial = 1 : size(dat_align.trial,2)
        istart = find(time > dat_align.time{itrial}(1),1,'first');
        l = length(dat_align.trial{itrial});
        n = mean(dat_align.trial{itrial},1);
        n = n - mean(n);
        n = n ./ max(abs(n));
        trialgrid(i,istart:istart+l-1) = n;
        i = i + 1;
    end
    
    % frequency analysis of whole trials peaks
    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.foilim      = [0.5,15];
    cfg.foilim      = [0,15];
    cfg.pad         = 'nextpow2';
    FFT             = ft_freqanalysis(cfg,dat_align);
    
    % time frequency analysis around peaks 
    cfg             = [];
    cfg.method      = 'mtmconvol';
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.foi         = 10:1:300;
    cfg.t_ftimwin   = ones(size(cfg.foi))*0.5;
    cfg.toi         = mintime:0.01:maxtime;
    TFR             = ft_freqanalysis(cfg,dat_align);
    
    % plot overview
    fig = figure;
    
    % plot peak detection results
    cfg                 = [];
    cfg.avgoverchan     = 'yes';
    dat_align_chanavg   = ft_selectdata(cfg,dat_align);
    
    subplot(4,1,1);
    
    hold;
    n = 10; % number of timecourses
    h = 150;
    for itrial = 1 : n
        plot(dat_filt_chanavg_align.time{itrial},dat_filt_chanavg_align.trial{itrial} + itrial*h,'linewidth',1,'color','r');
        plot(dat_filt_chanavg_align.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)),dat_filt_chanavg_align.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)) + itrial*h,'.r','markersize',h/7);
        plot(dat_align_chanavg.time{itrial},dat_align_chanavg.trial{itrial} + itrial*h,'color','k');
        if haspeak(itrial) == 0
            plot(dat_align_chanavg.time{itrial},dat_align_chanavg.trial{itrial} + itrial*h,'color','c');
        end
    end
    title('Peak detection examples');
    ylabel('Trials');
    xlabel('Time (s)');
    axis tight
    
    % plot grid of trial voltages
    subplot(5,1,2);
    image(trialgrid*255);
    colormap hot(255)
    
    mini = ceil(min(time));
    maxi = floor(max(time));
    tickindx = [];
    for i = mini : maxi
        [~, indx] = min(abs(time-i));
        tickindx = [tickindx, indx];
    end
    
    set(gca, 'XTick', tickindx);
    set(gca, 'XTickLabel', sprintf('%1.0f\n',time(tickindx)));
    title('Raw trial amplitudes over time');
    xlabel('Time (s)');
    ylabel('Trials');
    axis tight
    ax = axis;
    hcb = colorbar;
    title(hcb,'Microvolts')
    
    % plot TFR
    subplot(5,1,3);
    cfg = [];
    cfg.baseline = 'yes';
    cfg.baselinetype = 'relative';
    cfg.zlim = [0.9 1.1];
    ft_singleplotTFR(cfg,TFR);
    title('Time-frequency-representation');
    xlabel('Time (s)');
    ylabel('Freq');
    title(' ');
    hcb = colorbar;
    title(hcb,'Rel. change')
    
    % plot FFT
    subplot(5,1,4); hold;
    plot(FFT.freq,mean(FFT.powspctrm,1),'k');
    [ymax,imax] = max(mean(FFT.powspctrm,1));
    line([FFT.freq(imax),FFT.freq(imax)],[0,ymax],'color','r','linewidth',2);
    txt = ['\leftarrow ', num2str(FFT.freq(imax),3),'Hz'];
    text(FFT.freq(imax),ymax,txt,'color','r');
    title('Spectral analysis');
    xlabel('Hz');
    ylabel('Power');
    
    % plot histogram of durations
    subplot(5,1,5); hold;
    triallength = zeros(size(dat_align.trial,2),1);
    for itrial = 1 : size(dat_align.trial,2)
        triallength(itrial) = dat_align.time{itrial}(end) - poststim(imarker);
    end
    
    [VF1.histcount,VF1.edges,VF1.bin] = histcounts(triallength,0:stepsize(imarker):max(triallength)+1);
    bar(VF1.edges(2:end)-stepsize(imarker)/2,VF1.histcount,1,'facecolor','k')
    title('Distribution of durations');
    ylabel('Nr. of observations');
    xlabel('Duration (s)');
    [ymax,imax] = max(VF1.histcount);
    line([VF1.edges(imax)+stepsize(imarker)/2,VF1.edges(imax)+stepsize(imarker)/2],[0,ymax],'color','r','linewidth',2);
    txt = ['\leftarrow ', num2str(VF1.edges(imax)+stepsize(imarker)/2) ,'s'];
    text(VF1.edges(imax)+stepsize(imarker),ymax,txt,'color','r')
    
    [yleft,ileft] = find(VF1.histcount>0,1,'first');
    txt = ['\downarrow ', num2str(VF1.edges(ileft)+stepsize(imarker)/2) ,'s'];
    text(VF1.edges(ileft)+stepsize(imarker),yleft+20,txt,'color','r')
    
    [yright,iright] = find(VF1.histcount>0,1,'last');
    txt = ['\downarrow ', num2str(VF1.edges(iright)+stepsize(imarker)/2) ,'s'];
    text(VF1.edges(iright)+stepsize(imarker)/2,yright+20,txt,'color','r')
    
    % print to file
    set(fig,'PaperOrientation','landscape');
    set(fig,'PaperUnits','normalized');
    set(fig,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(imagesavedir,['overview_subplots_',label{imarker},'.pdf']),'-r300');
    clear trialgrid dat_chanavg  temp
    
    %% Make overview of average of peaks
    
    fname = fullfile(datasavedir,['trialdat_append_',label{imarker},'.mat']);
    if exist(fname,'file')       
        fprintf('Found data in %s, loading...',fname);
        load(fname,'trialdat_append');     
    else

        i = 1;
        for itrial = 1 : size(dat_align.trial,2)
            fprintf('Found %d peaks in trialnr. %d\n',size(locs_ac_sel_avg{itrial},2),itrial);
            
            if size(locs_ac_sel_avg{itrial},2) > 0
                if strcmp(avgmode(imarker),'all')
                    for ipeak = 1 : size(locs_ac_sel_avg{itrial},2)
                        cfg = [];
                        cfg.trials          = itrial;
                        center              = find(dat_align.time{itrial} > 0,1,'first');
                        cfg.begsample       = center - resamplefs*prepeak(imarker);
                        cfg.endsample       = center + resamplefs*postpeak(imarker);
                        
                        % don't read before or after data
                        if cfg.begsample > 0 && cfg.endsample <= size(dat_align.trial{itrial},2)
                            trialdat{i}         = ft_redefinetrial(cfg,dat_align);
                            i = i + 1;
                        end
                    end
                end
                if strcmp(avgmode(imarker),'first')
                    cfg = [];
                    cfg.trials          = itrial;
                    center              = find(dat_align.time{itrial} > 0,1,'first');
                    cfg.begsample       = center - resamplefs*prepeak(imarker);
                    cfg.endsample       = center + resamplefs*postpeak(imarker);
                    
                    % don't read before or after data
                    if cfg.begsample > 0 && cfg.endsample <= size(dat_align.trial{itrial},2)
                        trialdat{i}         = ft_redefinetrial(cfg,dat_align);
                        i = i + 1;
                    end
                end
            else
                fprintf('Did not find a peak in trialnr. %d\n',itrial);
            end
        end
        
        trialdat_append = ft_appenddata([],trialdat{:});
        clear trialdat
        trialdat_append = rmfield(trialdat_append,'cfg');
        
        % save to disk
        save(fname,'trialdat_append','locs_*','t1*');
        
    end
    
    cfg                 = [];
    cfg.vartrllength    = 2;
    avg                 = ft_timelockanalysis(cfg,trialdat_append);
    avg                 = rmfield(avg,'cfg');
    
    % time frequency analysis around peaks
    cfg             = [];
    cfg.method      = 'mtmconvol';
    cfg.output      = 'pow';
    cfg.taper       = 'hanning';
    cfg.foi         = 10:1:300;
    cfg.t_ftimwin   = ones(size(cfg.foi))*0.5;
    cfg.toi         = -prepeak(imarker):0.01:postpeak(imarker);
    freq            = ft_freqanalysis(cfg,trialdat_append);
       
    % plot overview of peaks
    figure;
    for ichan = 1:size(avg.label,1) % do not bother with artefacted channels or reference channel; 1 : size(P.fdir{idir},1)
        
        subplot(3,5,ichan); hold;
        [~,~,W,P] = findpeaks(avg.avg(ichan,:),avg.time,'MinPeakProminence',4,'Annotate','extents');
        findpeaks(avg.avg(ichan,:),avg.time,'MinPeakProminence',4,'Annotate','extents');
        legend('off');
        
        [~, ci] = max(W);
        title(sprintf('W=%.0fms, P=%.2f',W(ci)*1000,P(ci)));
        axis([-prepeak(imarker),postpeak(imarker), -100 150]);
        
        subplot(3,5,5+ichan); hold;
        patch([avg.time, avg.time(end:-1:1)],[avg.avg(ichan,:) - sqrt(avg.var(ichan,:)), avg.avg(ichan,end:-1:1) + sqrt(avg.var(ichan,end:-1:1))],[0 0 0],'facealpha',0.2,'edgecolor','none');
        plot(avg.time,avg.avg(ichan,:),'linewidth',1,'color',[0 0 0]);
        
        title([avg.label{ichan} ' (n=' num2str(size(trialdat_append.trial,2)) ')']);
        xlabel('ms');
        axis tight
        axis([-prepeak(imarker),postpeak(imarker), -100 150]);
        ax = axis;
        line([0,0],[ax(3) ax(4)],'LineStyle',':','color','k');
        xlabel('Time (s)');
        ylabel('Amplitude (microV)');
        
        subplot(3,5,10+ichan);
        cfg                 = [];
        cfg.channel         = 1;
        cfg.baseline        = 'yes';
        cfg.baselinetype    = 'relative';
        cfg.zlim            = 'maxmin';
        cfg.xlim            = [-prepeak(imarker), postpeak(imarker)];
        cfg.colormap        = hot;
        cfg.title           = ' ';
        cfg.colorbar        = 'no';
        ft_singleplotTFR(cfg,freq);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        
        % add scaled LFP line
        hold;
        ax = axis;
        scaled = (avg.avg(ichan,:) + abs(min(avg.avg(ichan,:))));
        scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
        plot(avg.time,scaled,'linewidth',1,'color',[0.1 0.5 1]);
        axis([-prepeak(imarker), postpeak(imarker), ax(3), ax(4)]);
        
    end
    
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    print(fig, '-dpdf', fullfile(imagesavedir,['overview_average_',label{imarker},'.pdf']),'-r300');
    
    clear freq trialdat_append* avg FFT dat_align
end
