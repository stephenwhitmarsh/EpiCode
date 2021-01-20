function [micromed_markers, stim_micro_aligned] = alignstimdata(cfg,varargin)

if nargin > 0
    if nargin > 1
        micromed_markers = varargin{1};
        markerdat_micro  = varargin{2};
    end        
else
    fprintf('Not enough input arguments!');
end

fname_output = fullfile(cfg.datasavedir,'aligned_stim_markers.mat');
if exist(fname_output,'file') && cfg.force == false
    load(fname_output,'micromed_markers','stim_micro_aligned');
else
    
    cfgtemp = [];
    stim_micro    = ft_appenddata(cfgtemp,markerdat_micro{:});
    clear markerdat_micro
    
    cfgtemp = [];
    cfgtemp.hpfilter    = 'yes';
    cfgtemp.hpfreq      = 100;
    cfgtemp.channel     = cfg.channel_micro;
    stim_micro          = ft_preprocessing(cfg,stim_micro);
    
    stim_micro_aligned  = stim_micro;
    for istim = 1 : size(stim_micro.trial,2)
        
        for ichan = 1 : size(stim_micro.label,1)
            v(ichan) = max(abs(stim_micro.trial{istim}(ichan,:)));
        end
        [~, chanindx] = max(v);
        
        if micromed_markers.Frequency(istim) == 50
            
            range(1) = floor(1/micromed_markers.Frequency(istim) * stim_micro.fsample) - 3;
            range(2) = ceil(1/micromed_markers.Frequency(istim) * stim_micro.fsample) + 3;
            
            t1 = find(stim_micro.time{istim} > 3.5,1,'first');
            t2 = find(stim_micro.time{istim} > 4,1,'first');
            
            [PKSp,LOCSp] = findpeaks(stim_micro.trial{istim}(chanindx,:),'annotate','peaks');
            [PKSn,LOCSn] = findpeaks(-stim_micro.trial{istim}(chanindx,:),'annotate','peaks');
            
            if length(PKSp) <= length(PKSn)
                LOCS = LOCSp;
                PKS = PKSp;
                m = max(stim_micro.trial{istim}(chanindx,t1:t2)) * 0.5;
                
                [PKSr,LOCSr] = findpeaks(stim_micro.trial{istim}(chanindx,t1:t2),'Npeaks',40,'MinPeakHeight',m,'annotate','peaks','SortStr','descend');
                prange(1) = mean(PKSr) - 4 * std(PKSr);
                prange(2) = mean(PKSr) + 4 * std(PKSr);
                
            else
                LOCS = LOCSn;
                PKS = PKSn;
                m = max(-stim_micro.trial{istim}(chanindx,t1:t2)) * 0.5;
                [PKSr,LOCSr] = findpeaks(-stim_micro.trial{istim}(chanindx,t1:t2),'Npeaks',40,'MinPeakHeight',m,'annotate','peaks','SortStr','descend');
                prange(1) = mean(PKSr) - 4 * std(PKSr);
                prange(2) = mean(PKSr) + 4 * std(PKSr);
            end
            
            LOCS = LOCS(PKS >= prange(1) & PKS <= prange(2));
            LOCS = LOCS(diff(LOCS) >= range(1) & diff(LOCS) <= range(2));
            
        else
            [PKSp,LOCSp] = findpeaks(stim_micro.trial{istim}(chanindx,:),'NPeaks',5,'SortStr','descend','annotate','peaks');
            [PKSn,LOCSn] = findpeaks(-stim_micro.trial{istim}(chanindx,:),'NPeaks',5,'SortStr','descend','annotate','peaks');
            LOCSp = sort(LOCSp);
            LOCSn = sort(LOCSn);
            
            if var(diff(LOCSp)) <= var(diff(LOCSn))
                LOCS = LOCSp;
            else
                LOCS = LOCSn;
            end
        end
        
        try
            stim_micro_aligned.time{istim}              = stim_micro_aligned.time{istim} - stim_micro.time{istim}(LOCS(1));
            micromed_markers.Sample_NL_aligned(istim)   = micromed_markers.Sample_NL(istim) + stim_micro.time{istim}(LOCS(1)) * stim_micro.fsample; % or +, never know...
            micromed_markers.Seconds_NL_aligned(istim)  = micromed_markers.Seconds_NL(istim) + stim_micro.time{istim}(LOCS(1)); % or +, never know...
            micromed_markers.Seconds_NL_shift(istim)    = stim_micro.time{istim}(LOCS(1));
            
        catch
            micromed_markers.Sample_NL_aligned(istim)   = micromed_markers.Sample_NL(istim); % or +, never know...
            micromed_markers.Seconds_NL_aligned(istim)  = micromed_markers.Seconds_NL(istim); % or +, never know...
            micromed_markers.Seconds_NL_shift(istim)    = 0;
            
        end
        %
        fig = figure; hold;
        plot(stim_micro_aligned.time{istim},stim_micro_aligned.trial{istim}(chanindx,:),'k');
        plot(stim_micro_aligned.time{istim}(LOCS),stim_micro_aligned.trial{istim}(chanindx,LOCS),'r.');
        axis tight
        
        titlestr1 = [stim_micro_aligned.label{chanindx},' ',num2str(micromed_markers.Frequency(istim)),'Hz ',num2str(micromed_markers.Current(istim,:)),'mA ',cell2mat(micromed_markers.Label(istim,:))];
        title(titlestr1);
        fig.Renderer='Painters';
        set(fig,'PaperOrientation','portrait');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        set(gca,'yticklabel','');
        fname = fullfile(cfg.imagesavedir,[titlestr1,'.pdf']);
        fprintf('\n%s\n',fname);
        print(fig, '-dpdf', fname,'-r600');
        close all

        %
    end
    save(fname_output,'micromed_markers','stim_micro_aligned','-v7.3');
    
end
