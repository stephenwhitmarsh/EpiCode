function [MuseStruct] = alignMuseMarkersXcorr(cfg, MuseStruct, force)

% ALIGNMUSEMARKERSXCORR determines timing shift of MUSE markers according
% to crosscorrelation with average
%
% use as
%    [MuseStruct] = alignMuseMarkersXcorr(cfg, MuseStruct, force)
%
% Multiple channels can be used simultaneously.
% Results are saved in a 'timeshift' field for every marker event in the MuseStruct.
% Later functions can use this to align the data (or just read data)
% Put last argument on TRUE if you want to function to re-compute the
% aligment, or FALSE to load from previous results.
%
% Necessary fields (as defined in _setparams function):
%
% config{1}.align.name        = e.g.: {'PSW','FA','ES'}; % names of patterns to align
% config{1}.align.channel     = e.g.: {'m1pNs_1', 'm1pNs_2', 'm1pNs_4', 'm1pNs_6', 'm1pNs_7', 'm1pNs_8'}; % channels to use for alignment
% config{1}.align.reref       = e.g.: 'no';
% config{1}.align.refmethod   = e.g.: 'bipolar';
% config{1}.align.latency.PSW = e.g.: [-0.1 2]; % time period to use for alignment
% config{1}.align.latency.FA  = e.g.: [-0.1 0.4];
% config{1}.align.latency.ES  = e.g.: [-0.1 0.4];
% config{1}.align.reject      = e.g.: 'BAD_cnt'; %fieldname in trialinfo used
%                               to reject trials, see readLFP and cfg.LFP.overlap.
% 
% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

% default settings
write                = ft_getopt(cfg.align, 'write', true);
cfg.visible          = ft_getopt(cfg, 'visible', 'on');
cfg.plotpages        = ft_getopt(cfg, 'plotpages', false);
cfg.align.refchannel = ft_getopt(cfg.align, 'refchannel', []);

fname       = fullfile(cfg.datasavedir,[cfg.prefix, 'MuseStruct_alignedXcorr.mat']);

% load results if requested (force = false) and result file exist

if nargin == 1
    if exist(fname, 'file')
        fprintf('Loading results xcorr alignment: %s\n', fname);
        
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                load(fname,'MuseStruct');
            catch ME
                err_count = err_count + 1;
                disp('Something went wrong loading the file. Trying again...')
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        return
    end
end

fname = fullfile(cfg.datasavedir,[cfg.prefix, 'MuseStruct_alignedXcorr.mat']);
if exist(fname,'file') && force == false
    fprintf('Loading results xcorr alignment: %s\n', fname);
    load(fname,'MuseStruct');
    return
end

fprintf('Aligning with xcorr\n');

% read LFP - based on epoch latency for plotting later
cfgtemp                 = rmfield(cfg, 'LFP');
cfgtemp.LFP.name        = cfg.align.name;
cfgtemp.LFP.channel     = [cfg.align.channel, cfg.align.refchannel];
cfgtemp.LFP.write       = false;
cfgtemp.LFP.lpfilter    = ft_getopt(cfg.align, 'lpfilter', 'no');
cfgtemp.LFP.lpfreq      = ft_getopt(cfg.align, 'lpfreq', '');

for markername = string(cfg.align.name)
    cfgtemp.epoch.pad.(markername) = 0;
end
LFP = readLFP(cfgtemp, MuseStruct, true);

for ipart = 1 : size(cfg.directorylist,2)
    
    for markername = string(cfg.align.name)
                
        if isempty(LFP{ipart}.(markername))
            continue
        end
        
        sel = ft_getopt(cfg.align, 'reject', []);
        cfgtemp = [];
        if ~isempty(sel)
            cfgtemp.trials = ~LFP{ipart}.(markername).trialinfo.(sel);
            LFP{ipart}.(markername) = ft_selectdata(cfgtemp, LFP{ipart}.(markername));          
        end
        
        % baseline correct & bipolar rereferencing
        cfgtemp                 = [];
        cfgtemp.demean          = ft_getopt(cfg.align, 'demean', 'yes');
        cfgtemp.reref           = ft_getopt(cfg.align, 'reref', 'no');
        cfgtemp.refmethod       = ft_getopt(cfg.align, 'refmethod', 'bipolar');
        cfgtemp.refchannel      = cfg.align.refchannel;
        dat                     = ft_preprocessing(cfgtemp, LFP{ipart}.(markername));
        
        %remove reference channels
        if ~isempty(cfgtemp.refchannel)
            cfgtemp         = [];
            cfgtemp.channel = cfg.align.channel;
            dat             = ft_selectdata(cfgtemp, dat);
        end
        
        % clear LFP for memory
        LFP{ipart}.(markername) = [];
        
        % select data
        if strcmp(cfg.align.latency.(markername), 'all')
            latency = [dat.time{1}(1), dat.time{1}(end)];
        else
            latency = cfg.align.latency.(markername);
        end
        
        cfgtemp                 = [];
        cfgtemp.latency         = latency;
        dat_timesel             = ft_selectdata(cfgtemp, dat);
        
        cfgtemp                 = [];
        dat_avg_orig            = ft_timelockanalysis(cfgtemp, dat);
        
        % put selected timeperiod in single matrix with concatinated channels
        fprintf('Preparing alignment...\n');
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
        
        % remove rejected trials of non-shifted for plotting later
        cfgtemp             = [];
        cfgtemp.trials      = find(~rejected);
        dat                 = ft_selectdata(cfgtemp, dat);

        %% Plot alignment
        fig = figure('visible', cfg.visible);
        fig.Renderer = 'Painters';
        
        subplot(1, 5, 1);
        imagesc(LFP_concatinated(~rejected,:));
        title(sprintf('Cleaned (%d-%d)', sum(~rejected), sum(rejected)));
        set(gca,'xtick', []);
        
        subplot(1, 5, 2);
        imagesc(shifted(~rejected, :));
        title('Aligned');
        set(gca,'xtick', []);
        
        subplot(1, 5, 3);
        scatter(nshift(~rejected) / dat.fsample * 1000, sum(~rejected):-1:1,'k.');
        set(gca,'ytick', []);
        xlabel('ms');
        title('Timeshift');
        axis tight;
        
        subplot(2, 5, [4 5]); hold on;
        plot(dat_avg_orig.time, dat_avg_orig.avg');
        axis tight;
        ax = axis;
        patch([latency(1), latency(2), latency(2), latency(1)], [ax(3), ax(3), ax(4), ax(4)], 'r', 'facealpha', 0.1, 'edgecolor', 'none');
        title('Original');
        xlim([dat.time{1}(1), dat.time{1}(end)]);
        
        subplot(2,5,[9 10]); hold on;
        plot(dat_avg_shifted.time, dat_avg_shifted.avg');
        axis tight
        ax = axis;
        patch([latency(1), latency(2), latency(2), latency(1)], [ax(3), ax(3), ax(4), ax(4)], 'r', 'facealpha', 0.1, 'edgecolor', 'none');
        title('Aligned');
        xlim([dat.time{1}(1), dat.time{1}(end)]);
        
        set(fig,'Renderer','Painters');
        set(fig,'PaperOrientation','landscape');
        set(fig,'PaperUnits','normalized');
        set(fig,'PaperPosition', [0 0 1 1]);
        
        % if it doesnt exist, create directory
        if ~exist(cfg.imagesavedir, 'dir')
            sprintf('Creating directory %s', cfg.imagesavedir);
            mkdir(cfg.imagesavedir);
        end
%         if ~exist(fullfile(cfg.imagesavedir, cfg.prefix(1:end-1)), 'dir')
%             sprintf('Creating directory %s', fullfile(cfg.imagesavedir, cfg.prefix(1:end-1)));
%             mkdir(fullfile(cfg.imagesavedir, cfg.prefix(1:end-1)));
%         end
        
        fname_fig = fullfile(cfg.imagesavedir, 'alignment_xcorr', strcat(cfg.prefix, 'p', num2str(ipart), '_', markername, '_alignmentXcorr.png'));
        isdir_or_mkdir(fileparts(fname_fig));
        exportgraphics(fig, fname_fig);
        
        %% Plot alignment
        if cfg.plotpages
            
            nrchan  = size(dat_shifted.label, 1);
            nrtrial = size(dat_shifted.trial,2);
            nrpage  = ceil(nrtrial / 100);
            
            for ipage = 1 : nrpage
                
                % get scaling parameter for all channels
                maxrange = [];
                for  ichan = 1 : size(dat_shifted.label, 1)
                    for itrial = 1 : size(dat_shifted.trial, 2)
                        y           = dat_shifted.trial{itrial}(ichan,:);
                        maxrange    = max(abs([maxrange, y]));
                    end
                end
                maxrange = maxrange / 4;
                
                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                for ichan = 1 : size(dat_shifted.label, 1)
                    
                    subplot(4, nrchan, [ichan ichan+nrchan ichan+nrchan*2 ichan+nrchan*3]); hold;
                    
                    i = 0;
                    for itrial = 1+(ipage-1)*100 : min(nrtrial, 100+(ipage-1)*100)
                        
                        plot(dat.time{itrial}, dat.trial{itrial}(ichan, :) + maxrange * i, 'color', [0.5, 0.5, 0.5]);
                        plot(dat_shifted.time{itrial}, dat_shifted.trial{itrial}(ichan, :) + maxrange * i, 'k');
                        i = i + 1;
                    end
                    
                    axis tight
                    ax = axis;
                    patch([latency(1), latency(2), latency(2), latency(1)], [ax(3), ax(3), ax(4), ax(4)], 'r', 'facealpha', 0.1, 'edgecolor', 'none');
                    plot([0, 0], [0, maxrange * 100], 'b:');
                    set(gca,'children',flipud(get(gca,'children')));
                    set(gca,'Yticklabels','')
                    xlabel('time (s)');
                    if ichan == 1
                        ylabel('trials');
                    end
                    title(dat_shifted.label{ichan}, 'interpreter', 'none');
                    
                end
                set(findall(fig, '-property', 'fontsize'), 'fontsize', 8);
                set(fig,'Renderer','Painters');
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                
                fname_fig = fullfile(cfg.imagesavedir, 'alignment_Xcorr', strcat(cfg.prefix, 'p', num2str(ipart), '_', markername, '_alignmentXcorr_page', num2str(ipage), '.png'));
                isdir_or_mkdir(fileparts(fname_fig));
                exportgraphics(fig, fname_fig);
                
            end
        end
        
        %% correct MuseStruct
        i = 0;
        for idir = 1 : length(MuseStruct{ipart})
            if ~isfield(MuseStruct{ipart}{idir},'markers')
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)), 'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime)
                continue
            end
            todelete = [];
            
            for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2)
                if any(dat.trialinfo.trialnr == ievent & dat.trialinfo.idir == idir)
                    i = i + 1;
                    timeshift = nshift_clean(i) / dat.fsample;
                    fprintf('Timeshifting %s #%d in part %d by %d samples (%0.3f seconds) \n', markername, ievent, ipart, nshift_clean(i),timeshift);
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).timeshift(ievent) = timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent)  = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) - timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent)     = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent) - seconds(timeshift); % direction should be tripple-checked
                else
                    fprintf('Removing event %d in part %d\n', ievent, ipart);
                    todelete = [todelete, ievent];
                end
            end
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(todelete) = [];
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(todelete)    = [];
            
        end % idir
        clear shifted shifted_clear shifted_clean_z dat dat_shifted
        close all
        
    end % imarker
    
end % ipart

% save data
if write
    save(fname,'MuseStruct');
end
