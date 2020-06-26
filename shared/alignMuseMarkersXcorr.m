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

% Paul (removeme) : see how to plot data

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_alignedXcorr.mat']);

write   = ft_getopt(cfg.align, 'write', true);
latency = ft_getopt(cfg.align, 'latency', 'all');

if exist(fname,'file') && force == false
    fprintf('*******************************\n');
    fprintf('** Loading results alignment **\n');
    fprintf('*******************************\n\n');
    load(fname,'MuseStruct');
else
    
    if force == true
        fprintf('*********************************\n');
        fprintf('** Forced redoing of alignment **\n');
        fprintf('*********************************\n\n');
    else
        fprintf('**************\n');
        fprintf('** Aligning **\n');
        fprintf('**************\n\n');
    end
    
    cfgtemp                 = rmfield(cfg,'LFP');
    cfgtemp                 = rmfield(cfgtemp,'epoch');
    %convert align cfg to lfp cfg
    for imarker = 1 : size(cfg.align.name,2)
        %correct latency if 'all'
        try if strcmp(cfg.align.latency{imarker}, 'all'), cfg.align.latency{imarker} = cfg.epoch.toi{imarker}; end, end
        %if find a latency for this marker, select it as a period to load
        %LFP. Else use the LFP parameter.
        try 
            cfgtemp.epoch.toi{imarker} = cfg.align.latency{imarker};
        catch
            cfgtemp.epoch.toi{imarker} = cfg.epoch.toi{imarker};
        end
        %no need of padding for alignment
        cfgtemp.epoch.pad{imarker} = 0;
    end
    cfgtemp.LFP.name        = cfg.align.name;
    cfgtemp.LFP.channel     = cfg.align.channel;
    [LFP]                   = readLFP(cfgtemp, MuseStruct, true);
    
    for ipart = 1:size(cfg.directorylist,2)
        
        for imarker = 1 : size(cfg.align.name,2)
            
            % baseline correct & bipolar rereferencing
            cfgtemp                 = [];
            cfgtemp.demean          = ft_getopt(cfg.align, 'demean', 'no');
            cfgtemp.baselinewindow  = ft_getopt(cfg.align, 'baselinewindow', 'no');
            cfgtemp.reref           = ft_getopt(cfg.align, 'reref', 'no'); 
            cfgtemp.refmethod       = ft_getopt(cfg.align, 'refmethod', 'bipolar');
            dat                     = ft_preprocessing(cfgtemp,LFP{ipart}{imarker});
            
%             % select data Not needed anymore because load only data of
%             interest
%             cfgtemp                 = [];
%             cfgtemp.latency         = ft_getopt(cfg.align, 'latency', 'all');
%             dat                     = ft_selectdata(cfgtemp,dat);
            
            % put all data in single matrix with concatinated channels
            fprintf('Preparing alignment');
            d = size(dat.trial{1});
            LFP_concatinated = nan(size(dat.trial,2),d(1)*d(2));
            for itrial = 1 : size(dat.trial,2)
                LFP_concatinated(itrial,:) = reshape(dat.trial{itrial}',1,d(1)*d(2));
            end
            
            % align data to closest-peak in cross correlation
            [shifted, nshift]       = alignXcorr(LFP_concatinated,5);
            
            % remove those that moved too much
            rejected                = (nshift < mean(nshift) - std(nshift)*3) | (nshift > mean(nshift) + std(nshift)*3);
            shifted_clean           = shifted(~rejected,:);
            
            % z-normalize data
            shifted_clean_z = nanznorm(shifted_clean);
            shifted_clean_z(isnan(shifted_clean_z)) = 0;
           
            try
            % draw figure
            fig = figure(1);
            fig.Renderer = 'Painters';
            
            subplot(2,5,1);
            imagesc(LFP_concatinated(~rejected,:));
            title('Cleaned');
            set(gca,'xtick',[]);
            
            subplot(2,5,6);
            plot(cumsum(repmat(dat.time{1},1,d(1))),LFP_concatinated(~rejected,:)');
%             plot(dat.time{1},LFP_concatinated(~rejected,:)');
            axis tight
            
            subplot(2,5,2);
            imagesc(LFP_concatinated(rejected,:));
            title('Rejected');
            set(gca,'xtick',[]);
            
            subplot(2,5,7);
            plot(cumsum(repmat(dat.time{1},1,d(1))),LFP_concatinated(rejected,:)');
%             plot(dat.time{1},LFP_concatinated(rejected,:)');
            axis tight
            
            subplot(2,5,3);
            imagesc(shifted(~rejected,:));
            title('Cleaned & Aligned');
            set(gca,'xtick',[]);
            
            subplot(2,5,8);
            plot(cumsum(repmat(dat.time{1},1,d(1))),shifted(~rejected,:)');
%             plot(dat.time{1},shifted(~rejected,:)');
            axis tight
            
            subplot(2,5,4);
            imagesc(shifted(rejected,:));
            title('Misaligned');
            set(gca,'xtick',[]);
            subplot(2,5,9);
            plot(cumsum(repmat(dat.time{1},1,d(1))),shifted(rejected,:)');
%             plot(dat.time{1},shifted(rejected,:)');
            axis tight
            
            subplot(2,5,5);
            imagesc(shifted_clean_z);
            title('Cleaned & Aligned (z)')
            set(gca,'xtick',[]);
            
            subplot(2,5,10);
            plot(cumsum(repmat(dat.time{1},1,d(1))),shifted_clean_z');
%             plot(dat.time{1},shifted_clean_z');
            axis tight
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpng', fullfile(cfg.imagesavedir,[cfg.prefix,'p',num2str(ipart),'alignmentXcorr.png']));
            
            catch 
                fprintf('Error with x-axis\n');
            end %REMOVEME
            
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
