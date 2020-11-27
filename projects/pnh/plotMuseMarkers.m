function plotMuseMarkers(MuseStruct)

%% append markers over files (converted to clocktime)

addpath /Users/stephen.whitmarsh/WhitmarshEpilepsy/
addpath /Users/stephen.whitmarsh/fieldtrip/ 
ft_defaults

patient_directory       = '/Users/stephen.whitmarsh';
patient_directory       = '/Volumes/epimicro/Donnees-analyses/Stephen/pat_02230_0674/eeg';
directory_searchstring  = '02230_2015-02-25_*';
directory_searchstring  = '02230_2015-02-*';
data_searchstring       = '*m1*.ncs';

MuseStruct = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);


MuseStruct_append = rmfield(MuseStruct{1},'filenames');

% for first directory
mrknames = fieldnames(MuseStruct{1}.markers);
for imarker = 1 : numel(mrknames)
    MuseStruct_append.markers.(mrknames{imarker}).sample = [];
    if isfield(MuseStruct{1}.markers.(mrknames{imarker}).events,'sample')
        for ievent = 1 : length(MuseStruct{1}.markers.(mrknames{imarker}).events)
            MuseStruct_append.markers.(mrknames{imarker}).sample = ...
                [MuseStruct_append.markers.(mrknames{imarker}).sample, ...
                MuseStruct{1}.markers.(mrknames{imarker}).events(ievent).sample];
        end
    end
    
    MuseStruct_append.markers.(mrknames{imarker}).sec = [];
    if isfield(MuseStruct{1}.markers.(mrknames{imarker}).events,'time')
        for ievent = 1 : length(MuseStruct{1}.markers.(mrknames{imarker}).events)
            MuseStruct_append.markers.(mrknames{imarker}).sec = ...
                [MuseStruct_append.markers.(mrknames{imarker}).sec, ...
                MuseStruct{1}.markers.(mrknames{imarker}).events(ievent).time];
        end
    end    
end

% adding next directories
for idir = 2 : size(MuseStruct,2)
    try % some might be empty
        mrknames = fieldnames(MuseStruct{idir}.markers);
        for imarker = 1 : numel(mrknames)
            if ~isfield(MuseStruct_append.markers,(mrknames{imarker}))
                MuseStruct_append.markers.(mrknames{imarker}) = MuseStruct{idir}.markers.(mrknames{imarker});
            end
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'clock')
                MuseStruct_append.markers.(mrknames{imarker}).clock = [];
            end
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}),'clock')
                MuseStruct_append.markers.(mrknames{imarker}).clock = ...
                    [MuseStruct_append.markers.(mrknames{imarker}).clock, ...
                    MuseStruct{idir}.markers.(mrknames{imarker}).clock];
            end
            
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'samples')
                MuseStruct_append.markers.(mrknames{imarker}).sample = [];
            end            
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}).events,'sample')
                for ievent = 1 : length(MuseStruct{idir}.markers.(mrknames{imarker}).events)
                MuseStruct_append.markers.(mrknames{imarker}).sample = ...
                    [MuseStruct_append.markers.(mrknames{imarker}).sample, ...
                    MuseStruct{idir}.markers.(mrknames{imarker}).events(ievent).sample];
                end
            end  
            
            if ~isfield(MuseStruct_append.markers.(mrknames{imarker}),'sec')
                MuseStruct_append.markers.(mrknames{imarker}).sec = [];
            end            
            if isfield(MuseStruct{idir}.markers.(mrknames{imarker}).events,'time')
                for ievent = 1 : length(MuseStruct{idir}.markers.(mrknames{imarker}).events)
                MuseStruct_append.markers.(mrknames{imarker}).sec = ...
                    [MuseStruct_append.markers.(mrknames{imarker}).sec, ...
                    MuseStruct{idir}.markers.(mrknames{imarker}).events(ievent).time];
                end
            end  
            
            
        end
    catch
    end
end

%% Plot overview

h = figure; 

subplot(8,1,1); hold;
fill([MuseStruct_append.markers.StartRecord.clock(1),MuseStruct_append.markers.StopRecord.clock(end),MuseStruct_append.markers.StopRecord.clock(end),MuseStruct_append.markers.StartRecord.clock(1)],[0 0 1 1],[1 0 0],'EdgeColor','none','facealpha',1);
for i = 1 : length(MuseStruct_append.markers.StartRecord.clock)
    fill([MuseStruct_append.markers.StartRecord.clock(i),MuseStruct_append.markers.StopRecord.clock(i),MuseStruct_append.markers.StopRecord.clock(i),MuseStruct_append.markers.StartRecord.clock(i)],[0 0 1 1],[0 1 0],'EdgeColor',[0 0 1],'facealpha',1);
end
title('Data');
axis tight
axx = xlim;
xlim(axx)

subplot(8,1,2); hold;
for i = 1 : length(MuseStruct_append.markers.BAD__START__.clock)
    fill([MuseStruct_append.markers.BAD__START__.clock(i),MuseStruct_append.markers.BAD__END__.clock(i),MuseStruct_append.markers.BAD__END__.clock(i),MuseStruct_append.markers.BAD__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('Artifacts');
xlim(axx)

subplot(8,1,3); hold;
for i = 1 : length(MuseStruct_append.markers.RR__START__.clock)
    fill([MuseStruct_append.markers.RR__START__.clock(i),MuseStruct_append.markers.RR__END__.clock(i),MuseStruct_append.markers.RR__END__.clock(i),MuseStruct_append.markers.RR__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('Fast ripples');
xlim(axx)

subplot(8,1,4); hold;
for i = 1 : length(MuseStruct_append.markers.VF1__START__.clock)
    fill([MuseStruct_append.markers.VF1__START__.clock(i),MuseStruct_append.markers.VF1__END__.clock(i),MuseStruct_append.markers.VF1__END__.clock(i),MuseStruct_append.markers.VF1__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF1');
xlim(axx)

subplot(8,1,5); hold;
for i = 1 : length(MuseStruct_append.markers.VF3__START__.clock)
    fill([MuseStruct_append.markers.VF3__START__.clock(i),MuseStruct_append.markers.VF3__END__.clock(i),MuseStruct_append.markers.VF3__END__.clock(i),MuseStruct_append.markers.VF3__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF3');
xlim(axx)

subplot(8,1,6); hold;
for i = 1 : length(MuseStruct_append.markers.VF4__START__.clock)
    fill([MuseStruct_append.markers.VF4__START__.clock(i),MuseStruct_append.markers.VF4__END__.clock(i),MuseStruct_append.markers.VF4__END__.clock(i),MuseStruct_append.markers.VF4__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('VF4');
xlim(axx)

subplot(8,1,7); hold;
for i = 1 : length(MuseStruct_append.markers.PP__START__.clock)
    fill([MuseStruct_append.markers.PP__START__.clock(i),MuseStruct_append.markers.PP__END__.clock(i),MuseStruct_append.markers.PP__END__.clock(i),MuseStruct_append.markers.PP__START__.clock(i)],[0 0 1 1],[0 0 0],'EdgeColor',[0 0 0],'facealpha',0.5);
end
title('PP');
xlim(axx)


subplot(8,1,8); hold;

[P.minutes, P.edges]    = discretize(MuseStruct_append.markers.P.clock,'minute');
[P.histcount,edges,bin] = histcounts(P.minutes,-0.5:1:P.minutes(end)+0.5);

bar(P.edges,P.histcount,1,'facecolor','k')
ylabel('Spikes per minute');
axis tight;
xlim(axx)

c = smooth(P.histcount);
plot(P.edges,c,'r','linewidth',1)

% 
% % set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)
% set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)
% set(gca,'Xtick',linspace(axx(1),axx(2),12),'XTickLabelRotation',45)

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'overview_visual.pdf','-r600');
% 
% % polar
% figure;
% 
% [P.hours, P.edges]    = discretize(P.clocktime_append(:,1),'hour');
% [P.histcount,edges,bin] = histcounts(P.hours,-0.5:1:P.hours(end)+0.5);
% 
% figure
% h = polarhistogram(P.minutes / 1440*pi*2,linspace(0,pi*2,1440));
% 
% figure
% h = polarhistogram(P.hours / 24*pi*2,linspace(0,pi*2,24));
% 
% h.DisplayStyle = 'stairs';
% 
% addpath ~/Documents/MATLAB/CircStat2012a/
% figure; circ_plot(P.hours / 24*pi*2,'density','r.',true,'linewidth',2,'color','r')
% figure; circ_plot(P.minutes / 1440*pi*2)

%% Text files for overview markers

startend = {'BAD__START__','BAD__END__';...
    'VF1__START__','VF1__END__';...
    'VF3__START__','VF3__END__';...
    'VF4__START__','VF4__END__';...
    'VF5__START__','VF5__END__';...
    'P','P';...
    'PP__START__','PP__END__';...
    'RR__START__','RR__END__'};

for imarker = 1 : length(startend)
    Starttime       = [];
    Endtime         = [];
    Startsample     = [];
    Endsample       = [];
    Duration        = [];
    Directory       = [];
    Event_start     = [];
    Event_end       = [];
    for idir = 1 : size(MuseStruct,2)
        try % some might be empty
            for ievent = 1 : size(MuseStruct{idir}.markers.(startend{imarker,1}).events,2)
                Event_start = [Event_start; startend{imarker,1}];
                Event_end   = [Event_end; startend{imarker,2}];
                Directory   = [Directory;   MuseStruct{idir}.directory];          
                Starttime   = [Starttime;   MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).time];
                Endtime     = [Endtime;     MuseStruct{idir}.markers.(startend{imarker,2}).events(ievent).time];
                Duration    = [Duration;    MuseStruct{idir}.markers.(startend{imarker,2}).events(ievent).time - MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).time];
                Startsample = [Startsample; MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).sample];
                Endsample   = [Endsample;   MuseStruct{idir}.markers.(startend{imarker,1}).events(ievent).sample];
            end
        catch
        end
    end
    overview = table(Event_start,Event_end,Directory,Starttime,Endtime,Duration,Startsample,Endsample);
    writetable(overview,['/Users/stephen.whitmarsh/Documents/VF/',startend{imarker,1} '-to-' startend{imarker,2} '.txt']); 
end

%% Place artefact file for Spyking-Circus in each data directory
for idir = 1 : size(MuseStruct,2)
    Starttime = [];
    Endtime = [];
    
    try
        if size(MuseStruct{idir}.markers.BAD__START__.events,2) > 0
            for ievent = 1 : size(MuseStruct{idir}.markers.BAD__START__.events,2)
                Starttime   = [Starttime;   MuseStruct{idir}.markers.BAD__START__.events(ievent).time];
                Endtime     = [Endtime;     MuseStruct{idir}.markers.BAD__END__.events(ievent).time];
            end
            
            filename = fullfile(MuseStruct{idir}.directory,'artifact.txt');
            fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
            % has to be in ms instead of seconds
            dlmwrite(filename,[Starttime*1000,Endtime*1000],'delimiter',' ','precision','%.4f'); % high precision, i.e sub-millisecond, does not make sense
        else
            fprintf('No artifact markers found in: %s\n',filename);
        end
    catch
        fprintf('No markerfile found in: %s\n',filename);
        
    end
end


%% extract VF1 epochs

startend    = { 'VF1__START__','VF1__END__';...
                'PP__START__','PP__END__'};
prestim     = [1,0.200];    % list of all marker onset and offsets. MIGHT REMOVE
poststim    = [1,0.850];
filenrs     = [1,4,6,7,8];  % do not bother with artefacted channels or reference channel
resamplefs  = 640;          % in Hz
imarker     = 1;

% for imarker = 1 : size(startend,1)
    
    for idir = 1:length(MuseStruct)
        
        hdr = ft_read_header(fullfile(MuseStruct{idir}.directory,MuseStruct{idir}.filenames(1).name));

        if isfield(MuseStruct{idir},'markers')
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
                
                % downsample data
                cfg = [];
                cfg.resamplefs              = resamplefs;
                cfg.demean                  = 'no';
                dirdat{idir}                = ft_resampledata(cfg,dirdat{idir} );
            end
        end
    end
% end

% concatinate data over trials
dat = ft_appenddata([],dirdat{[1:8,12:14]});
clear dirdat

% filtering data for peak detection
dat_filt = dat;
for itrial = 1 : size(dat.trial,2)
    for ichan = 1 : size(dat.label,1)
        dat_filt.trial{itrial}(ichan,:) = bandpassFilter(dat.trial{itrial}(ichan,:),resamplefs,1.0,5.0);
    end
end

% get maximum amplitude in baseline period 
cfg = [];
cfg.avgoverchan     = 'yes';
dat_filt_chanavg    = ft_selectdata(cfg,dat_filt); 
clear dat_filt
dat_chanavg         = ft_selectdata(cfg,dat);

% peak threshold detection based on baseline period 
bl_latency          =[-1, 0];
ac_latency          = [0, 30];
max_peaks_ac_list   = [];

for itrial = 1 : size(dat.trial,2)
    t1_bl_indx(itrial)                  = find(dat_filt_chanavg.time{itrial} > bl_latency(1),1,'first');
    t2_bl_indx(itrial)                  = find(dat_filt_chanavg.time{itrial} < bl_latency(2),1,'last');
    t1_ac_indx(itrial)                  = find(dat_filt_chanavg.time{itrial} > ac_latency(1),1,'first');
    t2_ac_indx(itrial)                  = find(dat_filt_chanavg.time{itrial} < ac_latency(2),1,'last');
    
    [peaks_bl{itrial},~]                = findpeaks(dat_filt_chanavg.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
    max_peaks_bl(itrial)                = max(peaks_bl{itrial});
    [peaks_ac{itrial},locs_ac{itrial}]  = findpeaks(dat_filt_chanavg.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));
    max_peaks_ac(itrial)                = max(peaks_ac{itrial});
    max_peaks_ac_list                   = [max_peaks_ac_list, peaks_ac{itrial}];
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

% plot peak detection results
figure; hold;
n = 50;
h = 150;
for itrial = 1 : n % size(dat.trial,2)
    plot(dat_chanavg.time{itrial},dat_chanavg.trial{itrial} + itrial*h);
    plot(dat_filt_chanavg.time{itrial},dat_filt_chanavg.trial{itrial} + itrial*h,'linewidth',1,'color','k');
    plot(dat_filt_chanavg.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1),dat_filt_chanavg.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1) + itrial*h,'.r','markersize',h/5);
    plot(dat_filt_chanavg.time{itrial}(locs_ac_sel_trl{itrial}+t1_ac_indx(itrial)-1),dat_filt_chanavg.trial{itrial}(locs_ac_sel_trl{itrial}+t1_ac_indx(itrial)-1) + itrial*h,'.g','markersize',h/7);
end
axis tight
clear dat_filt_chanavg

% segment data around peaks
prepeak = 0.5;
postpeak = 0.5;
i = 1;
for itrial = 1 : size(dat.trial,2)
    if size(locs_ac_sel_avg{itrial},2) > 0
        for ipeak = 1 : size(locs_ac_sel_avg{itrial},2)
            cfg = [];
            cfg.trials          = itrial;
            center              = locs_ac_sel_avg{itrial}(ipeak)+t1_ac_indx(itrial)-1;
            cfg.begsample       = center - resamplefs*prepeak;
            cfg.endsample       = center + resamplefs*postpeak;
            
            % don't read before or after data
            if cfg.begsample >= 0 && cfg.endsample <= size(dat.trial{itrial},2)
                trialdat{i}         = ft_redefinetrial(cfg,dat);
                
                % redefine time axis around peak
                trialdat{i}.time{1} = trialdat{i}.time{1} - trialdat{i}.time{1}(ceil(size(trialdat{i}.time{1},2)/2));
                i = i + 1;
            end
        end
    end
end

trialdat_append = ft_appenddata([],trialdat{:});
clear trialdat
trialdat_append = rmfield(trialdat_append,'cfg');


cfg = [];
trialavg = ft_timelockanalysis(cfg,trialdat_append);
trialavg = rmfield(trialavg,'cfg');

figure;
plot(trialavg.time,trialavg.avg);


% time frequency analysis
cfg             = [];
cfg.method      = 'mtmconvol';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = 10:5:200;
cfg.t_ftimwin   = ones(size(cfg.foi))*0.5;
cfg.toi         = -prepeak:0.01:postpeak;
freq            = ft_freqanalysis(cfg,trialdat_append);



















%%

figure; hold;
plot(trialdat{1}.time{1},trialdat{1}.trial{1});
% 
plot(trialdat{2}.time{1},trialdat{2}.trial{1});

% % get peaks in each trial that surpass maximum in trial-based baseline
% latency = 0;
% 
% for itrial = 1 : size(dat_filt.trial,2)
%     t1_indx = find(dat_filt_chanavg.time{itrial} > latency(1),1,'first');
%     t2_indx = size(dat_filt_chanavg.time{itrial},2);
%     [peakamp{itrial},peakindx{itrial}] = findpeaks(dat_filt_chanavg.trial{itrial}(t1_indx:t2_indx));
%     peakindx{itrial} = peakindx{itrial} + t1_indx -1;
% %     max_bl(itrial) = max(PKS);
% end   
% 
% 
% 
% 
% 
% % dat_filt_bl(itrial) = ;
% end
% 
% figure; hist(q)
% 
% figure; hold;
% for itrial = 1 : 15 % size(dat.trial,2)
%     plot(dat.time{itrial}(1,640:end-640),dat.trial{itrial}(1,640:end-640));
% end
% figure;   plot(dat.time{1}(1,640:end-640),dat.trial{1}(1,640:end-640));
% 
% 
% cfg             = [];
% cfg.lpfilter    = 'yes';
% cfg.lpfreq      = 5;
% dat_filt        = ft_preprocessing(cfg,dat); 
% 
% for itrial = 1 : size(dat_filt.trial,2)
%     
% [pks,locs,widths,proms] = findpeaks(dat_filt.trial{itrial},x);
% 
%         
%             % align trial times to max (peak)
%             corrwindow = 0.05; % maximum correction in seconds
%             for itrial = 1 : size(filedat{ifile}.trial,2)
%                 leftside                            = find(filedat{ifile}.time{itrial} > -corrwindow,1,'first');
%                 rightside                           = find(filedat{ifile}.time{itrial} <  corrwindow,1,'last');
%                 [ymax, imax]                        = max(filedat{ifile}.trial{itrial}(leftside:rightside));
%                 imax                                = imax + leftside;
%                 
%                 % shift timing of marker file (Muse) as well, only have to
%                 % do it once on the first file
%                 timediff = filedat{ifile}.time{itrial}(imax);
%                 
%                 if ifile == 1
%                     P.clocktime_corr{idir}(itrial,1) = P.clocktime{idir}(itrial,1) + timediff;
%                     P.clocktime_corr{idir}(itrial,2) = P.clocktime{idir}(itrial,2) + timediff;
%                 end
%                 
%                 filedat{ifile}.time{itrial}         = filedat{ifile}.time{itrial} - filedat{ifile}.time{itrial}(imax);
%                 filedat{ifile}.trialinfo(itrial,3)  = - filedat{ifile}.time{itrial}(imax);
%             end
%             
%             % crop trials around peak-annotations
%             cfg = [];
%             cfg.toilim      = [-prestim+corrwindow,poststim-corrwindow]; % in seconds
%             filedat{ifile}  = ft_redefinetrial(cfg,filedat{ifile});
%             
%             % demean and downsample data
%             cfg = [];
%             cfg.resamplefs      = 640;
%             cfg.demean          = 'yes';
%             cfg.baselinewindow  = [-0.150,-0.050];
%             filedat{ifile}      = ft_resampledata(cfg,filedat{ifile});
%             
%             % bugfix for floating-point comparisons
%             for itrial = 1 : size(filedat{ifile}.trial,2)
%                 filedat{ifile}.time{itrial} = linspace(-prestim+corrwindow,poststim-corrwindow,length(filedat{ifile}.time{itrial}));
%             end
%             
%         end
%     end
% end
