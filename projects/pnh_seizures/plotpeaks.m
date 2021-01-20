% 

patient_directory       = '/Volumes/epimicro/Donnees-analyses/Stephen/2230 files with markers/';
directory_searchstring  = '02230_*';
data_searchstring       = '*m1*.ncs';

eventtype_start         = 'P';
eventtype_end           = 'P';
[P.clocktime,P.clocktime_append,P.sample,P.data,P.ldir,P.fdir] = extractMuseMarkers(patient_directory,directory_searchstring,data_searchstring,eventtype_start,eventtype_end);

hdr = ft_read_header(fullfile(P.fdir{1}(1).folder,P.fdir{1}(1).name));

for idir = 1 : size(P.ldir,1)
    
    % in case of no events in markerfile
    if ~isempty(P.sample{idir})
        
        for ifile = [1,4,6,7,8] % do not bother with artefacted channels or reference channel; 1 : size(P.fdir{idir},1)
            
            % filter on the whole trace to avoid aliasing artefacts
            cfg = [];
            cfg.dataset     = fullfile(P.fdir{idir}(ifile).folder,P.fdir{idir}(ifile).name);
            cfg.hpfilter    = 'no';
            cfg.hpfreq      = 0.1;
            cfg.lpfilter    = 'yes';
            cfg.lpfreq      = 300;
            temp            = ft_preprocessing(cfg);
            
            % segment into trials around peak-annotations, larger to allow
            % re-alignment
            cfg = [];
            prestim         = 0.200; % in seconds, including 50 ms for realignment
            poststim        = 0.850;
            cfg.trl         = P.sample{idir};
            cfg.trl(:,1)    = P.sample{idir}(:,1) - round(hdr.Fs * prestim);
            cfg.trl(:,2)    = P.sample{idir}(:,1) + round(hdr.Fs * poststim);
            cfg.trl(:,3)    = ones(size(P.sample{idir},1),1) * - round(hdr.Fs * prestim);
            cfg.trl(:,4)    = idir;
            cfg.trl(:,5)    = ifile;
            filedat{ifile}  = ft_redefinetrial(cfg,temp);
            clear temp
            
            % truncate label
            filedat{ifile}.label{1} = filedat{ifile}.label{1}(end-6:end);
            
            % align trial times to max (peak)
            corrwindow = 0.05; % maximum correction in seconds
            for itrial = 1 : size(filedat{ifile}.trial,2)
                leftside                            = find(filedat{ifile}.time{itrial} > -corrwindow,1,'first');
                rightside                           = find(filedat{ifile}.time{itrial} <  corrwindow,1,'last');
                [ymax, imax]                        = max(filedat{ifile}.trial{itrial}(leftside:rightside));
                imax                                = imax + leftside;
                
                % shift timing of marker file (Muse) as well, only have to
                % do it once on the first file
                timediff = filedat{ifile}.time{itrial}(imax);
                
                if ifile == 1 
                    P.clocktime_corr{idir}(itrial,1) = P.clocktime{idir}(itrial,1) + timediff;
                    P.clocktime_corr{idir}(itrial,2) = P.clocktime{idir}(itrial,2) + timediff;                  
                end

                filedat{ifile}.time{itrial}         = filedat{ifile}.time{itrial} - filedat{ifile}.time{itrial}(imax);
                filedat{ifile}.trialinfo(itrial,3)  = - filedat{ifile}.time{itrial}(imax);
            end
            
            % crop trials around peak-annotations
            cfg = [];
            cfg.toilim      = [-prestim+corrwindow,poststim-corrwindow]; % in seconds
            filedat{ifile}  = ft_redefinetrial(cfg,filedat{ifile});
            
            % demean and downsample data
            cfg = [];
            cfg.resamplefs      = 640;
            cfg.demean          = 'yes';
            cfg.baselinewindow  = [-0.150,-0.050];
            filedat{ifile}      = ft_resampledata(cfg,filedat{ifile});
            
            % bugfix for floating-point comparisons
            for itrial = 1 : size(filedat{ifile}.trial,2)
                filedat{ifile}.time{itrial} = linspace(-prestim+corrwindow,poststim-corrwindow,length(filedat{ifile}.time{itrial}));
            end
            
        end
        
        cfg = [];
        cfg.appenddim = 'chan';
        dirdat{idir} = ft_appenddata(cfg,filedat{[1,4,6,7,8]});
        clear filedat
    end
%     
%     % add time-corrected markers
%     name_mrk        = fullfile(P.ldir(idir).folder,P.ldir(idir).name,'Events.mrk');
%     field_name      = 'Pc';
%     field_sample    = P.clocktime_corr{idir};
%     field_trial     = zeros(size(field_sample));
%     field_comment   = 'Added by Stephen';
%     field_nrsamples = size(field_sample,1);
%     addMuseMarker(name_mrk,field_name,field_trial,field_sample,field_nrsamples,field_comment)

end

dat             = ft_appenddata([],dirdat{:});
clear dirdat

cfg             = [];
cfg.method      = 'mtmconvol';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = 10:200;
cfg.t_ftimwin   = ones(size(cfg.foi))*0.1;
cfg.toi         = -prestim+corrwindow : 0.01 : poststim-corrwindow;
freq            = ft_freqanalysis(cfg,dat);
avg             = ft_timelockanalysis([],dat);

figure;
for ichan = 1:size(avg.label,1) % do not bother with artefacted channels or reference channel; 1 : size(P.fdir{idir},1)
    subplot(2,5,ichan); hold;
    patch([avg.time, avg.time(end:-1:1)],[avg.avg(ichan,:) - sqrt(avg.var(ichan,:)), avg.avg(ichan,end:-1:1) + sqrt(avg.var(ichan,end:-1:1))],[0 0 0],'facealpha',0.2,'edgecolor','none');
    plot(avg.time,avg.avg(ichan,:),'linewidth',1,'color',[0 0 0]);
    
    title([avg.label{ichan} ' (n=' num2str(size(dat.trial,2)) ')']);
    xlabel('ms');
    axis tight
    ax = axis;
    axis([-0.200, 0.850, -100 300]);
    ax = axis;
    line([0,0],[ax(3) ax(4)],'LineStyle',':','color','k');
    xlabel('Time (s)');
    ylabel('Amplitude (microV)');
    
    subplot(2,5,5+ichan);
    cfg                 = [];
    cfg.channel         = 1;
    cfg.baseline        = [-0.1 -0.05];
    cfg.baselinetype    = 'relative';
    cfg.zlim            = 'maxmin';
    cfg.xlim            = [-0.200, 0.850];
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
    axis([-0.200, 0.850, ax(3), ax(4)]);
    
end

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'overview_peaks_long.pdf','-r600');