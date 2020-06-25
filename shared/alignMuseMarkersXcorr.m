function [MuseStruct] = alignMuseMarkersXcorr(cfg, MuseStruct, force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [MuseStruct] = alignMuseMarkersXcorr(cfg, MuseStruct, force)
%
% Automatically determines timing shift of events read from MUSE, according
% to crosscorrelation between channels. Multiple channels can be used simultaneously.
% Results are saved in a 'timeshift' field for every marker event in the MuseStruct.
% Later functions can use this to align the data (or just read data)
% Put last argument on TRUE if you want to function to re-compute the
% aligment, or FALSE to load from previous results.
%
% Necessary fields (as defined in _setparams function):
% 
% cfg.epoch.toi{1}            = [-0.5  1];
% cfg.epoch.toi{2}            = [-0.5  1];
% cfg.epoch.pad               = {0.5, 0.5, 0.5};
% cfg.align.name              = e.g.: {'Hspike','SpikeDetect'};                               % Name of markers/patterns to align
% cfg.align.channel           = e.g.: {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};    % Channels to use for alignment
% cfg.align.demean            = e.g.: 'yes';
% cfg.align.baselinewindow    = e.g.: [-0.2 -0.1];
% cfg.align.reref             = e.g.: 'yes';
% cfg.align.refmethod         = e.g.: 'bipolar';
% cfg.align.latency           = e.g.: [-0.1, 0.2];                                            % timeperiod to use for crosscorrelation
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_alignedXcorr.mat']);

write   = ft_getopt(cfg.align, 'write', true);
latency = ft_getopt(cfg.align, 'latency', 'all');

if exist(fname,'file') && force == false
    fprintf('*******************************\n');
    fprintf('** Loading results alignment **\n');
    fprintf('*******************************\n\n');
    load(fname,'MuseStruct');
    return
end

fprintf('**************\n');
fprintf('** Aligning **\n');
fprintf('**************\n\n');

cfgtemp                 = rmfield(cfg,'LFP');
cfgtemp.LFP.name        = cfg.align.name;
cfgtemp.LFP.channel     = cfg.align.channel;
cfgtemp.LFP.write       = false;
[LFP]                   = readLFP(cfgtemp, MuseStruct, true);

for ipart = 1:size(cfg.directorylist,2)
    
    for imarker = 1 : size(cfg.align.name,1)
        
        % baseline correct & bipolar rereferencing
        cfgtemp                 = [];
        cfgtemp.demean          = ft_getopt(cfg.align, 'demean', 'no');
        cfgtemp.baselinewindow  = ft_getopt(cfg.align, 'baselinewindow', 'no');
        cfgtemp.reref           = ft_getopt(cfg.align, 'reref', 'no');
        cfgtemp.refmethod       = ft_getopt(cfg.align, 'refmethod', 'bipolar');
        dat                     = ft_preprocessing(cfgtemp, LFP{ipart}{imarker});
        
        % select data
        if strcmp('latency','all')
            latency = [dat.time{1}(1), dat.time{1}(end)];
        end
        
        cfgtemp                 = [];
        cfgtemp.latency         = latency;
        dat_timesel             = ft_selectdata(cfgtemp, dat);
        dat_avg_orig            = ft_timelockanalysis([], dat);

        % put selected timeperiod in single matrix with concatinated channels
        fprintf('Preparing alignment');
        d = size(dat_timesel.trial{1});
        LFP_concatinated = nan(size(dat_timesel.trial,2),d(1)*d(2));
        for itrial = 1 : size(dat.trial,2)
            LFP_concatinated(itrial,:) = reshape(dat_timesel.trial{itrial}',1,d(1)*d(2));
        end
        
        % align data to closest-peak in cross correlation
        [shifted, nshift]   = alignXcorr(LFP_concatinated, 10);
        
        % remove those that moved too much
        rejected            = (nshift < mean(nshift) - std(nshift)*3) | (nshift > mean(nshift) + std(nshift)*3);
        shifted_clean       = shifted(~rejected,:);
        nshift_clean        = nshift(~rejected);
        
        % make time-shifted average
        cfgtemp             = [];
        cfgtemp.trials      = find(~rejected);
        dat_shifted         = ft_selectdata(cfgtemp, dat);
        for itrial = 1 : size(nshift_clean,2)
            dat_shifted.time{itrial} = dat_shifted.time{itrial} + nshift_clean(itrial) / dat_shifted.fsample;
        end
        dat_avg_shifted     = ft_timelockanalysis([], dat_shifted);

        % draw figure
        fig = figure('visible','off');
        fig.Renderer = 'Painters';
        
        subplot(1,5,1);
        imagesc(LFP_concatinated(~rejected,:));
        title(sprintf('Cleaned (%d removed)',sum(rejected)));
        set(gca,'xtick',[]);
        
        subplot(1,5,2);
        imagesc(shifted(~rejected,:));
        title('Aligned');
        set(gca,'xtick',[]);
        
        subplot(1,5,3);
        scatter(nshift(~rejected) / dat.fsample * 1000,sum(~rejected):-1:1,'k.');
        set(gca,'ytick',[]);
        xlabel('ms');
        title('Timeshift');
        axis tight;
        
        subplot(2,5,[4 5]);
        plot(dat_avg_orig.time, dat_avg_orig.avg');
        xlim([dat_avg_orig.time(1), dat_avg_orig.time(end)]);
        ax = axis;
        patch([latency(1), latency(2), latency(2), latency(1)],[ax(3), ax(3), ax(4), ax(4)],'r','facealpha',0.1,'edgecolor','none');        
        title('Original');
        
        subplot(2,5,[9 10]);
        plot(dat_avg_shifted.time, dat_avg_shifted.avg');
        xlim([dat_avg_orig.time(1), dat_avg_orig.time(end)]);
        patch([latency(1), latency(2), latency(2), latency(1)],[ax(3), ax(3), ax(4), ax(4)],'r','facealpha',0.1,'edgecolor','none');        
        title('Aligned');
        
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'alignmentXcorr.png']));
        
        close all
        clear shifted shifted_clear shifted_clean_z
        
        % correct MuseStruct
        i = 1;
        for idir = 1:length(MuseStruct{ipart})
            if ~isfield(MuseStruct{ipart}{idir},'markers')
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers,(cfg.muse.startend{imarker,1}))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                continue
            end
            todelete = [];
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)
                if any(dat.trialinfo(:,1) == ievent & dat.trialinfo(:,3) == idir)
                    timeshift = nshift(i) * 1/dat.fsample;
                    fprintf('Timeshifting event %d in part %d by %d samples (%0.3f seconds) \n',ievent,ipart,nshift(i),timeshift);
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(ievent)      = timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent) + timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(ievent)          = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(ievent) + seconds(timeshift);
                    i = i + 1;
                else
                    fprintf('Removing event %d in part %d\n',ievent,ipart);
                    todelete = [todelete, ievent];
                end
            end
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(todelete)       = [];
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(todelete)          = [];
        end % idir
    end % imarker
end % ipart

% save data
if write
    save(fname,'MuseStruct');
end
