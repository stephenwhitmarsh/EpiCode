function analyseVF_macro(MuseStruct_micro,MuseStruct_macro)

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/analysis/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/
ft_defaults

%% load or create MuseStruct

if isempty(MuseStruct_micro)
    if exist('MuseStruct_micro.mat','file')
        load('MuseStruct_micro','MuseStruct_micro');
    else
        patient_directory       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/';
        directory_searchstring  = '02230_2015-02-*';
        data_searchstring       = '*m1pNs*.ncs';
        MuseStruct_micro        = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);
        save('MuseStruct_micro','MuseStruct_micro');
    end
end

if isempty(MuseStruct_macro)
    if exist('MuseStruct_macro.mat','file')
        load('MuseStruct_macro','MuseStruct_macro');
    else
        patient_directory       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/';
        directory_searchstring  = '02230_2015-02-*';
        data_searchstring       = '*_1pNs*.ncs';
        MuseStruct_macro        = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);
        save('MuseStruct_macro','MuseStruct_macro');
    end
end

%% parameters

startend                = { 'VF1__START__','VF1__END__'; ...
    'RR__START__','RR__END__'; 'P','P'; 'PP__START__','PP__END__'};             % Muse Marker names to analyze
label                   = { 'VF1','VF2','P','PP'};                              % Own names
prestim                 = [1,1,1,2.0];                                          % list of onset timing with respect to start-marker (s)
poststim                = [1,1,1,4.0];                                          % list of offset timing with respect to end-marker (s)
prepeak                 = [0.5,0.8,0.2,1];                                      % list of onset timing to align to (first) peak (s)
postpeak                = [0.8,0.8,0.8,1.5];                                    % list of offset timing to align to (first) peak (s)
baseline                = {'no','no','no',[-1,-0.5]};                           % whether to baseline correct after peak detection, if so, the baseline period (s)
alignmode               = {'first','first','max','first'};                      % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
avgmode                 = {'first','first','all','first'};                      % whether to average only the first peak, or all peaks {'first' or 'all}
micro_labels            = {'m1pNs_1','m1pNs_4','m1pNs_6','m1pNs_7','m1pNs_8'};                                      % do not bother with artefacted channels or reference channel (nr. 2+8 is reference in pat_02230)
macro_labels            = {'_1pNs_1','_1pNs_2','_1pNs_3','_1pNs_4','_1pNs_5','_1pNs_6','_1pNs_7','_1pNs_8'};
resamplefs              = 640;                                                  % resample data to speed up and reduce memory (Hz), and align micro and macro
stepsize                = [0.2, 0.05, 0.01, 0.01];                              % for histogram of duration of events in overview (s)
lpfreq                  = [4,4,40,40];                                          % lowpass filter freq to smooth peak detection (Hz)
peak_search_toi_ac      = [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];             % active period in which to search for peaks
peak_search_toi_bl      = [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];             % baseline period in which to search for peaks
datasavedir             = '/network/lustre/iss01/charpier/stephen.whitmarsh/analysis/data';    % where to write data
imagesavedir            = '/network/lustre/iss01/charpier/stephen.whitmarsh/analysis/images';  % where to print images

%% Make overview of data as segmented by markers aligned to (first) peak

for imarker = 1 : size(startend,1)
    
    fname = fullfile(datasavedir,['dat_align_',label{imarker},'.mat']);
    
    %     if exist(fname,'file')
    %
    %         fprintf('Found data in %s, loading...',fname);
    %         load(fname,'dat_align','dat_filt_chanavg_align','locs_ac_sel_trl','locs_ac_sel_avg','t1_ac_indx');
    %
    %     else
    
    fprintf('Aligned data not found - making');
    
    hasmarker = zeros(length(MuseStruct_micro),1);
    
    for idir = 1:length(MuseStruct_micro)
        if isfield(MuseStruct_micro{idir},'markers')
            if isfield(MuseStruct_micro{idir}.markers,(startend{imarker,1}))
                if ~isempty(MuseStruct_micro{idir}.markers.(startend{imarker,1}).events)
                    
                    % select MICRO files
                    micro_filenrs = [];
                    for ifile = 1 : size(MuseStruct_micro{idir}.filenames,1)
                        for ilabel = 1 : size(micro_labels,2)
                            if ~isempty(strfind(MuseStruct_micro{idir}.filenames(ifile).name,micro_labels{ilabel}))
                                micro_filenrs = [micro_filenrs, ifile];
                            end
                        end
                    end
                    
                    % select MACRO files
                    macro_filenrs = [];
                    for ifile = 1 : size(MuseStruct_macro{idir}.filenames,1)
                        for ilabel = 1 : size(macro_labels,2)
                            if ~isempty(strfind(MuseStruct_macro{idir}.filenames(ifile).name,macro_labels{ilabel}))
                                macro_filenrs = [macro_filenrs, ifile];
                            end
                        end
                    end
                    
                    % load trials for selected MICRO channels
                    hdr_micro = ft_read_header(fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames(ifile).name));
                    for ifile = micro_filenrs
                        
                        % create Fieldtrip trl
                        Startsample     = [];
                        Endsample       = [];
                        for ievent = 1 : size(MuseStruct_micro{idir}.markers.(startend{imarker,1}).events,2)
                            Startsample         = [Startsample; MuseStruct_micro{idir}.markers.(startend{imarker,1}).events(ievent).sample - prestim(imarker)*hdr_micro.Fs];
                            Endsample           = [Endsample;   MuseStruct_micro{idir}.markers.(startend{imarker,2}).events(ievent).sample + poststim(imarker)*hdr_micro.Fs];
                            Offset              = -ones(size(Endsample)) * prestim(imarker) * hdr_micro.Fs;
                        end
                        
                        cfg                             = [];
                        cfg.trl                         = [Startsample, Endsample, Offset];
                        cfg.trl                         = cfg.trl(Startsample > 0 & Endsample < hdr_micro.nSamples,:); % so not to read before BOF or after EOFs
                        cfg.dataset                     = fullfile(MuseStruct_micro{idir}.directory,MuseStruct_micro{idir}.filenames(ifile).name);
                        filedat_micro{ifile}            = ft_preprocessing(cfg);
                        filedat_micro{ifile}.label{1}   = filedat_micro{ifile}.label{1}(end-6:end); % truncate label
                    end
                    
                    % load trials for selected MACRO channels
                    hdr_macro = ft_read_header(fullfile(MuseStruct_macro{idir}.directory,MuseStruct_macro{idir}.filenames(ifile).name));
                    for ifile = macro_filenrs
                        
                        % create Fieldtrip trl
                        Startsample     = [];
                        Endsample       = [];
                        for ievent = 1 : size(MuseStruct_macro{idir}.markers.(startend{imarker,1}).events,2)
                            Startsample         = [Startsample; MuseStruct_macro{idir}.markers.(startend{imarker,1}).events(ievent).sample - prestim(imarker)*hdr_macro.Fs];
                            Endsample           = [Endsample;   MuseStruct_macro{idir}.markers.(startend{imarker,2}).events(ievent).sample + poststim(imarker)*hdr_macro.Fs];
                            Offset              = -ones(size(Endsample)) * prestim(imarker) * hdr_macro.Fs;
                        end
                        
                        cfg                             = [];
                        cfg.trl                         = [Startsample, Endsample, Offset];
                        cfg.trl                         = cfg.trl(Startsample > 0 & Endsample < hdr_macro.nSamples,:); % so not to read before BOF or after EOFs
                        cfg.dataset                     = fullfile(MuseStruct_macro{idir}.directory,MuseStruct_macro{idir}.filenames(ifile).name);
                        filedat_macro{ifile}            = ft_preprocessing(cfg);
                        filedat_macro{ifile}.label{1}   = filedat_macro{ifile}.label{1}(end-6:end); % truncate label
                    end
                    
                    % concatinate channels, separately for MICRO/MACRO
                    dirdat_micro{idir}                = ft_appenddata([],filedat_micro{micro_filenrs});
                    dirdat_macro{idir}                = ft_appenddata([],filedat_macro{macro_filenrs});
                    clear filedat
                    
                    % downsample data and baseline correction
                    cfg                         = [];
                    cfg.resamplefs              = resamplefs;
                    if strcmp(baseline(imarker),'no')
                        cfg.demean              = 'no';
                    else
                        cfg.demean              = 'yes';
                        cfg.baselinewindow      = baseline{imarker};
                    end
                    dirdat_micro{idir}          = ft_resampledata(cfg,dirdat_micro{idir} );
                    dirdat_macro{idir}          = ft_resampledata(cfg,dirdat_macro{idir} );
                    
                    % flag for averaging
                    hasmarker(idir) = 1;
                end
            end
        end
    end
    
    % concatinate data over trials
    dat_micro = ft_appenddata([],dirdat_micro{find(hasmarker)});
    dat_macro = ft_appenddata([],dirdat_macro{find(hasmarker)});
    clear dirdat*
    
    % for backup/debug
    save(fullfile(datasavedir,['dat_raw_',label{imarker},'.mat']),'dat_micro','dat_macro');
    
    % filtering data for peak detection
    dat_micro_filt = dat_micro;
    for itrial = 1 : size(dat_micro.trial,2)
        for ichan = 1 : size(dat_micro.label,1)
            dat_micro_filt.trial{itrial}(ichan,:) = bandpassFilter(dat_micro.trial{itrial}(ichan,:),resamplefs,1.0,lpfreq(imarker));
        end
    end
    
    % do peak detection on average over micro channels
    cfg = [];
    cfg.avgoverchan        = 'yes';
    dat_micro_filt_chanavg = ft_selectdata(cfg,dat_micro_filt); clear dat_filt
    
    % peak threshold detection based on baseline period
    for itrial = 1 : size(dat_micro.trial,2)
        t1_bl_indx(itrial)                 = find(dat_micro_filt_chanavg.time{itrial} > peak_search_toi_bl(imarker,1),1,'first');
        t2_bl_indx(itrial)                 = find(dat_micro_filt_chanavg.time{itrial} < peak_search_toi_bl(imarker,2),1,'last');
        t1_ac_indx(itrial)                 = find(dat_micro_filt_chanavg.time{itrial} > peak_search_toi_ac(imarker,1),1,'first');
        t2_ac_indx(itrial)                 = find(dat_micro_filt_chanavg.time{itrial} < (dat_micro_filt_chanavg.time{itrial}(end)-poststim(imarker)+peak_search_toi_ac(imarker,2)),1,'last'); % with respect to end of trial
        
        [peaks_bl{itrial},~]               = findpeaks(dat_micro_filt_chanavg.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
        if ~isempty(peaks_bl{itrial})
            max_peaks_bl(itrial)               = max(peaks_bl{itrial});
        else
            max_peaks_bl(itrial)               = nan;
        end
        
        [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_micro_filt_chanavg.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));
        if ~isempty(peaks_ac{itrial})
            max_peaks_ac(itrial)               = max(peaks_ac{itrial});
        else
            max_peaks_ac(itrial)     = nan;
        end
    end
    
    for itrial = 1 : size(dat_micro.trial,2)
        
        % threshold based on median of max peaks in all baseline periods
        peaks_sel_avg               = peaks_ac{itrial} > nanmean(max_peaks_bl) * 0.5;
        peaks_ac_sel_avg{itrial}    = peaks_ac{itrial}(peaks_sel_avg);
        locs_ac_sel_avg{itrial}     = locs_ac{itrial}(peaks_sel_avg);
        
        % threshold based on median of max peaks in trail-by-trial baseline
        peaks_sel_trl               = peaks_ac{itrial} > max(peaks_bl{itrial}) * 0.5;
        peaks_ac_sel_trl{itrial}    = peaks_ac{itrial}(peaks_sel_trl);
        locs_ac_sel_trl{itrial}     = locs_ac{itrial}(peaks_sel_trl);
        
    end
    
    % find trials that have no peaks detected (within 500 ms)
    haspeak = ones(1,size(dat_micro.trial,2));
    for itrial = 1 : size(dat_micro.trial,2)
        if isempty(locs_ac_sel_avg{itrial})
            haspeak(itrial) = 0;
            fprintf('Could not find peak in trial number %d. Will remove trial...\n',itrial);
        end
    end
    haspeak = find(haspeak);
    
    % select only trials that had a peak detected
    cfg                     = [];
    cfg.trials              = haspeak;
    dat_micro               = ft_selectdata(cfg,dat_micro);
    dat_macro               = ft_selectdata(cfg,dat_macro);
    dat_micro_filt_chanavg  = ft_selectdata(cfg,dat_micro_filt_chanavg);
    locs_ac_sel_avg         = locs_ac_sel_avg(haspeak);
    t1_ac_indx              = t1_ac_indx(haspeak);
    
    % align data to first peak
    dat_micro_align         = dat_micro;
    dat_macro_align         = dat_macro;
    dat_micro_filt_chanavg_align = dat_micro_filt_chanavg;
    
    for itrial = 1 : size(dat_micro_filt_chanavg.trial,2)
        
        % find relevant peak according to method in MICRO channels
        if strcmp(alignmode(imarker),'nearest')
            [~, ip] = min(abs(dat_micro_filt_chanavg.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1))); % peak-time closes to zero
        elseif strcmp(alignmode(imarker),'max')
            [~, ip] = max(dat_micro_filt_chanavg.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1)); % max peak
        elseif strcmp(alignmode(imarker),'first')
            ip = 1; % first peak after start search toi
        end
        
        % align time axis to relevant peak
        dat_micro_align.time{itrial} = dat_micro.time{itrial} - dat_micro.time{itrial}(locs_ac_sel_avg{itrial}(ip)+t1_ac_indx(itrial)-1);
        dat_macro_align.time{itrial} = dat_macro.time{itrial} - dat_macro.time{itrial}(locs_ac_sel_avg{itrial}(ip)+t1_ac_indx(itrial)-1);
        dat_micro_filt_chanavg_align.time{itrial} = dat_micro_filt_chanavg.time{itrial} - dat_micro_filt_chanavg.time{itrial}(locs_ac_sel_avg{itrial}(ip)+t1_ac_indx(itrial)-1);
        
    end
    clear dat_micro dat_macro
    
    % save to disk
    save(fname);
    % end
    
    % map all trials onto time-independant grid (average over channels)
    mintime                 = 0;
    maxtime                 = 0;
    maxsamples              = 0;
    for itrial = 1 : size(dat_micro_align.trial,2)
        if dat_micro_align.time{itrial}(1) < mintime
            mintime = dat_micro_align.time{itrial}(1);
        end
        if dat_micro_align.time{itrial}(end) > maxtime
            maxtime = dat_micro_align.time{itrial}(end);
        end
        if length(dat_micro_align.time{itrial}) > maxsamples
            maxsamples = length(dat_micro_align.time{itrial});
        end
    end
    
    time      = mintime : 1/resamplefs : maxtime;
    trialgrid_micro = zeros(size(dat_micro_align.trial,2),length(time));
    trialgrid_macro = zeros(size(dat_macro_align.trial,2),length(time));
    
    i = 1;
    for itrial = 1 : size(dat_micro_align.trial,2)
        istart = find(time > dat_micro_align.time{itrial}(1),1,'first');
        l = length(dat_micro_align.trial{itrial});
        n = mean(dat_micro_align.trial{itrial},1);
        n = n - mean(n);
        n = n ./ max(abs(n));
        trialgrid_micro(i,istart:istart+l-1) = n;
        i = i + 1;
    end
    
    i = 1;
    for itrial = 1 : size(dat_macro_align.trial,2)
        istart = find(time > dat_macro_align.time{itrial}(1),1,'first');
        l = length(dat_macro_align.trial{itrial});
        n = mean(dat_macro_align.trial{itrial},1);
        n = n - mean(n);
        n = n ./ max(abs(n));
        trialgrid_macro(i,istart:istart+l-1) = n;
        i = i + 1;
    end
    
    % frequency analysis of whole trials peaks
    cfg                 = [];
    cfg.method          = 'mtmfft';
    cfg.output          = 'pow';
    cfg.taper           = 'hanning';
    cfg.foilim          = [0.5,15];
    cfg.foilim          = [1,15];
    cfg.pad             = 'nextpow2';
    FFT_micro_trials    = ft_freqanalysis(cfg,dat_micro_align);
    FFT_macro_trials    = ft_freqanalysis(cfg,dat_macro_align);
    
    % time frequency analysis around peaks
    cfg                 = [];
    cfg.method          = 'mtmconvol';
    cfg.output          = 'pow';
    cfg.taper           = 'hanning';
    cfg.pad             = 'nextpow2';
    cfg.foi             = 10:1:300;
    cfg.t_ftimwin       = ones(size(cfg.foi))*0.5;
    cfg.toi             = mintime:0.01:maxtime;
    TFR_micro_trials    = ft_freqanalysis(cfg,dat_micro_align);
    TFR_macro_trials    = ft_freqanalysis(cfg,dat_macro_align);
    
    % plot overview
    close all
    fig = figure;
    
    % plot peak detection results
    cfg                 = [];
    cfg.avgoverchan     = 'yes';
    dat_micro_align_chanavg   = ft_selectdata(cfg,dat_micro_align);
    dat_macro_align_chanavg   = ft_selectdata(cfg,dat_macro_align);
    
    n = 10; % number of timecourses
    h = 150;
    
    subplot(5,2,1);
    hold;
    for itrial = 1 : n
        plot(dat_micro_filt_chanavg_align.time{itrial},dat_micro_filt_chanavg_align.trial{itrial} + itrial*h,'linewidth',1,'color','r');
        plot(dat_micro_filt_chanavg_align.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)),dat_micro_filt_chanavg_align.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)) + itrial*h,'.r','markersize',h/7);
        plot(dat_micro_align_chanavg.time{itrial},dat_micro_align_chanavg.trial{itrial} + itrial*h,'color','k');
    end
    title('Peak detection examples (micro)');
    ylabel('Trials');
    xlabel('Time (s)');
    axis tight
    
    subplot(5,2,2);
    hold;
    for itrial = 1 : n
        plot(dat_macro_align_chanavg.time{itrial},dat_macro_align_chanavg.trial{itrial} + itrial*h,'color','k');
    end
    title('Peak detection examples (macro)');
    ylabel('Trials');
    xlabel('Time (s)');
    axis tight
    
    % plot grid of trial voltages (micro)
    subplot(5,2,3);
    image(trialgrid_micro*255);
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
    
    % plot grid of trial voltages (macro)
    subplot(5,2,4);
    image(trialgrid_macro*255);
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
    
    % plot TFR micro
    subplot(5,2,5);
    cfg = [];
    cfg.baseline = 'yes';
    cfg.baselinetype = 'relative';
    % cfg.zlim = [0.9 1.1];
    ft_singleplotTFR(cfg,TFR_micro_trials);
    title('Time-frequency-representation');
    xlabel('Time (s)');
    ylabel('Freq');
    title(' ');
    hcb = colorbar;
    title(hcb,'Rel. change')
    
    % plot TFR macro
    subplot(5,2,6);
    cfg = [];
    cfg.baseline = 'yes';
    cfg.baselinetype = 'relative';
    % cfg.zlim = [0.9 1.1];
    ft_singleplotTFR(cfg,TFR_macro_trials);
    title('Time-frequency-representation');
    xlabel('Time (s)');
    ylabel('Freq');
    title(' ');
    hcb = colorbar;
    title(hcb,'Rel. change')
    
    % plot FFT micro
    subplot(5,2,7); hold;
    plot(FFT_micro_trials.freq,mean(FFT_micro_trials.powspctrm,1),'k');
    [ymax,imax] = max(mean(FFT_micro_trials.powspctrm,1));
    line([FFT_micro_trials.freq(imax),FFT_micro_trials.freq(imax)],[0,ymax],'color','r','linewidth',2);
    txt = ['\leftarrow ', num2str(FFT_micro_trials.freq(imax),3),'Hz'];
    text(FFT_micro_trials.freq(imax),ymax,txt,'color','r');
    title('Spectral analysis');
    xlabel('Hz');
    ylabel('Power');
    
    % plot FFT macro
    subplot(5,2,8); hold;
    plot(FFT_macro_trials.freq,mean(FFT_macro_trials.powspctrm,1),'k');
    [ymax,imax] = max(mean(FFT_macro_trials.powspctrm,1));
    line([FFT_macro_trials.freq(imax),FFT_macro_trials.freq(imax)],[0,ymax],'color','r','linewidth',2);
    txt = ['\leftarrow ', num2str(FFT_macro_trials.freq(imax),3),'Hz'];
    text(FFT_macro_trials.freq(imax),ymax,txt,'color','r');
    title('Spectral analysis');
    xlabel('Hz');
    ylabel('Power');
    
    % plot histogram of durations
    subplot(5,1,5); hold;
    triallength = zeros(size(dat_micro_align.trial,2),1);
    for itrial = 1 : size(dat_micro_align.trial,2)
        triallength(itrial) = dat_micro_align.time{itrial}(end) - poststim(imarker);
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
    clear trialgrid dat_chanavg dat_filt_chanavg temp
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make overview of average of peaks %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fname_data = fullfile(datasavedir,['trialdat_append_',label{imarker},'.mat']);
    %     if exist(fname,'file')
    %         fprintf('Found data in %s, loading...',fname);
    %         load(fname);
    %     else
    
    i = 1;
    for itrial = 1 : size(dat_micro_align.trial,2)
        fprintf('Found %d peaks in trialnr. %d\n',size(locs_ac_sel_avg{itrial},2),itrial);
        
        if size(locs_ac_sel_avg{itrial},2) > 0
            if strcmp(avgmode(imarker),'all')
                for ipeak = 1 : size(locs_ac_sel_avg{itrial},2)
                    cfg = [];
                    cfg.trials          = itrial;
                    center              = find(dat_micro_align.time{itrial} > 0,1,'first');
                    cfg.begsample       = center - resamplefs*prepeak(imarker);
                    cfg.endsample       = center + resamplefs*postpeak(imarker);
                    
                    % don't read before or after data
                    if cfg.begsample > 0 && cfg.endsample <= size(dat_micro_align.trial{itrial},2)
                        trialdat_micro{i}         = ft_redefinetrial(cfg,dat_micro_align);
                        trialdat_macro{i}         = ft_redefinetrial(cfg,dat_macro_align);
                        i = i + 1;
                    end
                end
            end
            if strcmp(avgmode(imarker),'first')
                cfg = [];
                cfg.trials          = itrial;
                center              = find(dat_micro_align.time{itrial} > 0,1,'first');
                cfg.begsample       = center - resamplefs*prepeak(imarker);
                cfg.endsample       = center + resamplefs*postpeak(imarker);
                
                % don't read before or after data
                if cfg.begsample > 0 && cfg.endsample <= size(dat_micro_align.trial{itrial},2)
                    trialdat_micro{i}         = ft_redefinetrial(cfg,dat_micro_align);
                    trialdat_macro{i}         = ft_redefinetrial(cfg,dat_macro_align);
                    i = i + 1;
                end
            end
        else
            fprintf('Did not find a peak in trialnr. %d\n',itrial);
        end
    end
    
    trialdat_micro_append = ft_appenddata([],trialdat_micro{:});
    trialdat_macro_append = ft_appenddata([],trialdat_macro{:});
    
    clear trialdat
    trialdat_micro_append = rmfield(trialdat_micro_append,'cfg');
    trialdat_macro_append = rmfield(trialdat_macro_append,'cfg');
    
    %     end
    
    cfg                 = [];
    cfg.vartrllength    = 2;
    avg_micro           = ft_timelockanalysis(cfg,trialdat_micro_append);
    avg_macro           = ft_timelockanalysis(cfg,trialdat_macro_append);
    avg_micro           = rmfield(avg_micro,'cfg');
    avg_macro           = rmfield(avg_macro,'cfg');
    
    % time frequency analysis around peaks
    %     cfg             = [];
    %     cfg.method      = 'mtmconvol';
    %     cfg.output      = 'pow';
    %     cfg.taper       = 'hanning';
    %     cfg.foi         = 10:1:300;
    %     cfg.t_ftimwin   = ones(size(cfg.foi))*0.5;
    %     cfg.toi         = -prepeak(imarker):0.01:postpeak(imarker);
    %     freq_micro      = ft_freqanalysis(cfg,trialdat_micro_append);
    %     freq_macro      = ft_freqanalysis(cfg,trialdat_macro_append);
    
    % save to disk
    %     save(fname_data);
    
    % plot overview of peaks
    close all
    fig = figure;
    for ichan = 1:size(avg_micro.label,1) % do not bother with artefacted channels or reference channel; 1 : size(P.fdir{idir},1)
        
        subplot(3,5,ichan); hold;
        [~,~,W,P] = findpeaks(avg_micro.avg(ichan,:),avg_micro.time,'MinPeakProminence',4,'Annotate','extents');
        findpeaks(avg_micro.avg(ichan,:),avg_micro.time,'MinPeakProminence',4,'Annotate','extents');
        legend('off');
        
        [~, ci] = max(W);
        title(sprintf('W=%.0fms, P=%.2f',W(ci)*1000,P(ci)));
        axis([-prepeak(imarker),postpeak(imarker), -100 150]);
        
        subplot(3,5,5+ichan); hold;
        patch([avg_micro.time, avg_micro.time(end:-1:1)],[avg_micro.avg(ichan,:) - sqrt(avg_micro.var(ichan,:)), avg_micro.avg(ichan,end:-1:1) + sqrt(avg_micro.var(ichan,end:-1:1))],[0 0 0],'facealpha',0.2,'edgecolor','none');
        plot(avg_micro.time,avg_micro.avg(ichan,:),'linewidth',1,'color',[0 0 0]);
        
        title([avg_micro.label{ichan} ' (n=' num2str(size(trialdat_micro_append.trial,2)) ')']);
        xlabel('ms');
        axis tight
        axis([-prepeak(imarker),postpeak(imarker), -100 150]);
        ax = axis;
        line([0,0],[ax(3) ax(4)],'LineStyle',':','color','k');
        xlabel('Time (s)');
        ylabel('Amplitude (microV)');
        
        %         subplot(3,5,10+ichan);
        %         cfg                 = [];
        %         cfg.channel         = 1;
        %         cfg.baseline        = 'yes';
        %         cfg.baselinetype    = 'relative';
        %         cfg.zlim            = 'maxmin';
        %         cfg.xlim            = [-prepeak(imarker), postpeak(imarker)];
        %         cfg.colormap        = hot;
        %         cfg.title           = ' ';
        %         cfg.colorbar        = 'no';
        %         ft_singleplotTFR(cfg,freq_micro);
        %         xlabel('Time (s)');
        %         ylabel('Frequency (Hz)');
        
        % add scaled LFP line
        hold;
        ax = axis;
        scaled = (avg_micro.avg(ichan,:) + abs(min(avg_micro.avg(ichan,:))));
        scaled = (scaled ./ max(scaled)) * (ax(4) - ax(3)) + ax(3);
        plot(avg_micro.time,scaled,'linewidth',1,'color',[0.1 0.5 1]);
        axis([-prepeak(imarker), postpeak(imarker), ax(3), ax(4)]);
        
        %     end
        
        h=gcf;
        set(h,'PaperOrientation','landscape');
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpdf', fullfile(imagesavedir,['overview_average_',label{imarker},'.pdf']),'-r300');
        
        clear freq trialdat_append* avg FFT dat_align
    end
    close all
end