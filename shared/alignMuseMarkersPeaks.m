function [MuseStruct] = alignMuseMarkersPeaks(cfg, varargin)

% ALIGNMUSEMARKERS aligns MUSE markers according to morphological features
%
% Use as
%    [MuseStruct] = alignMuseMarkers(cfg, MuseStruct, force)
%
% Results are saved in a 'timeshift' field for every marker event in MuseStruct
% Later functions can use this to align the data.
% Put last argument on TRUE if you want to function to re-compute the
% aligment, or FALSE to load from previous results.
%
% Necessary fields (as defined in _setparams function):
%
% cfg.name                      = e.g.: {'Hspike','SpikeDetect'};                                % Name by which to call the events, can be different from name in MUSE
% cfg.prefix                    = e.g.: 'P1-';                                                   % unique identifier for subject, prefixed to data and image files
% cfg.muse.startend             = e.g.: {'Hspike','Hspike'; 'SpikeDetect','SpikeDetect'};        % start and end Muse marker
% cfg.patientdir                = e.g.: fullfile(rootpath_data,       'pat_02711_1193');         % patient directory
% cfg.rawdir                    = e.g.: fullfile(rootpath_data,       'pat_02711_1193', 'eeg');  % raw data directory
% cfg.datasavedir               = e.g.: fullfile(rootpath_analysis,   'data',   'hspike');       % where to write data
% cfg.imagesavedir              = e.g.: fullfile(rootpath_analysis,   'images', 'hspike');       % where to print images
% cfg.alignpeak.name                = e.g.: {'Hspike','SpikeDetect'};                                % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% cfg.alignpeak.channel             = e.g.: {'_HaT2_1','_HaT2_1'};                                   % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% cfg.alignpeak.abs                 = e.g.: {'no','no'};                                             % whether to determine peaks on absolute (both max and min)
% cfg.alignpeak.method              = e.g.: {'crawlback','min'};                                     % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% cfg.alignpeak.maxtimeshift.(markername)        = e.g.: {0.05,0.05};                                % reject the trial as artefact if timeshift is bigger than this value (in seconds)
% cfg.alignpeak.filter              = e.g.: {'bp','bp'};                                             % what filter to use (bp,hp,lp)
% cfg.alignpeak.freq                = e.g.: {[1, 10],[1, 40]};                                       % lowpass filter freq to smooth peak detection (Hz)
% cfg.alignpeak.demean.(markername)              = e.g.: {'yes', 'yes'}                              % whether to apply baseline correction (usefull for detecting begining of event)
% cfg.alignpeak.hilbert             = e.g.: {'no','no'};                                             % whether to first create a hilbert transform (to determine peak of frequency response rather than timecourse
% cfg.alignpeak.thresh.method       = e.g.: {'trial','trial'}                                        % which peak baseline value is selected for thresholding : 'trial' 'medianbl','both'
% cfg.alignpeak.thresh.value        = e.g.: [1, 0];                                                  % only peaks that are more than x percent above threshold, 0 to disable
% cfg.alignpeak.toiplot{1}          = e.g.: [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.alignpeak.toiactive{1}        = e.g.: [-0.05,  0.150];                                         % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% cfg.alignpeak.toibaseline{1}      = e.g.: [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.alignpeak.toiplot{2}          = e.g.: [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.alignpeak.toiactive{2}        = e.g.: [-0.05,  0.05];                                          % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% cfg.alignpeak.toibaseline{2}      = e.g.: [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
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

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned.mat']);

% check input arguments, and default to false;
if length(varargin) > 3
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 3 optional inputs');
end
if length(varargin) == 2
    MuseStruct = varargin{1};
    force = varargin{2};
end
if length(varargin) == 1
    force = false;
end
if ~exist(fname,'file') && length(varargin) == 1
    force = true;
end

if exist(fname,'file') && force == false
    fprintf('Loading results alignment\n');
    load(fname, 'MuseStruct');
    return
end

if force == true
    fprintf('Forced redoing of peak-alignment\n');
else
    fprintf('Peak-aligning\n');
end

%get format to adapt script for each format specificities :
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% Go through different parts
for ipart = 1 : size(cfg.directorylist,2)

    for markername = string(cfg.alignpeak.name)

        %optional cfg fields
        cfg.alignpeak.findbegin                 = ft_getopt(cfg.alignpeak, 'findbegin', []);
        cfg.alignpeak.findbegin.(markername)    = ft_getopt(cfg.alignpeak.findbegin, convertStringsToChars(markername), 'no');
        cfg.alignpeak.hilbert                   = ft_getopt(cfg.alignpeak,'hilbert',[]);
        cfg.alignpeak.hilbert.(markername)      = ft_getopt(cfg.alignpeak.hilbert, convertStringsToChars(markername), 'no');
        cfg.alignpeak.demean                    = ft_getopt(cfg.alignpeak,'demean',[]);
        cfg.alignpeak.demean.(markername)       = ft_getopt(cfg.alignpeak.demean, convertStringsToChars(markername), 'no');
        cfg.alignpeak.reref                     = ft_getopt(cfg.alignpeak,'reref','no');
        cfg.alignpeak.notch                     = ft_getopt(cfg.alignpeak,'notch','no');

        % find data directories that have the required event
        markerindx = [];
        for idir = 1 : size(MuseStruct{ipart},2)
            if isfield(MuseStruct{ipart}{idir},'markers')
                if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startmarker.(markername))
                    if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)),'synctime')
                        if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime)
                            fprintf('Found %d event timings for %s \n',size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock,2),cfg.muse.startmarker.(markername));
                            markerindx = [markerindx idir];
                        end
                    end
                end
            end
        end
        if isempty(markerindx)
            fprintf('Did not find any events for marker %s\n',markername);
        end

        % loop over all data directories that have the marker
        for idir = markerindx

            % find datafilename corresponding to requested electrode
            if isNeuralynx
                filelist = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
                channelnr = 0;
                for ifile = 1 : size(filelist,1)
                    if strfind(filelist(ifile).name,cfg.alignpeak.channel.(markername))
                        fprintf('Found channel with pattern "%s" in %s\n',cfg.alignpeak.channel.(markername),filelist(ifile).name);
                        dataset = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(ifile).name);
                        channelnr = channelnr + 1;
                    end
                end
                if channelnr == 0
                    fprintf('Did not find of channel pattern %s!\n',cfg.alignpeak.channel.(markername));
                end
                if channelnr > 1
                    fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignpeak.channel.(markername));
                end

            elseif isMicromed
                dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);

            elseif isBrainvision
                dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);

            end

            %load data and header
            hdr = ft_read_header(dataset);

            cfgtemp             = [];

            if strcmp(cfg.alignpeak.reref,'yes')
                cfgtemp.reref       = cfg.alignpeak.reref;
                cfgtemp.rerefmethod = cfg.alignpeak.rerefmethod;
                cfgtemp.refchannel  = cfg.alignpeak.refchannel;
            end

            if strcmp(cfg.alignpeak.notch,'yes')
                cfgtemp.bsfilter       = 'yes';
                cfgtemp.bsfreq         = [49 51];
            end

            cfgtemp.dataset     = dataset;
            dat_sel             = ft_preprocessing(cfgtemp);

            if isMicromed || isBrainvision %choose channel
                cfgtemp         = [];
                cfgtemp.channel = cfg.alignpeak.channel.(markername);
                dat_sel         = ft_selectdata(cfgtemp,dat_sel);
                if size(dat_sel.trial{1},1) == 0
                    fprintf('Did not find of channel pattern %s!\n',cfg.alignpeak.channel.(markername));
                end
                if size(dat_sel.trial{1},1) > 1
                    fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.alignpeak.channel.(markername));
                end
            end

            % flip data if required
            if strcmp(cfg.alignpeak.method.(markername),'min') || strcmp(cfg.alignpeak.method.(markername),'firstmin') || strcmp(cfg.alignpeak.method.(markername),'lastmin') || strcmp(cfg.alignpeak.method.(markername),'nearestmin')
                dat_sel.trial{1} = -dat_sel.trial{1};
            end

            % filtering data of selected channel for peak detection
            switch cfg.alignpeak.filter.(markername)

                case 'lp'
                    fprintf('Lowpass Filtering data for alignment. This can take a while \n');
                    cfgtemp                     = [];
                    cfgtemp.lpfilter            = 'yes';
                    cfgtemp.lpfreq              = cfg.alignpeak.freq.(markername);
                    cfgtemp.lpfilttype          = 'fir';
                    dat_filt                    = ft_preprocessing(cfgtemp,dat_sel);

                case 'bp'
                    fprintf('Bandpass Filtering the data for alignment. This can take a while \n');
                    dat_filt                    = dat_sel;
                    dat_filt.trial{1}           = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.alignpeak.freq.(markername)(1),cfg.alignpeak.freq.(markername)(2));

                    if strcmp(cfg.alignpeak.method.(markername), 'crawlback')
                        dat_filt40              = dat_sel;
                        dat_filt40.trial{1}     = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.alignpeak.freq.(markername)(1),40);
                    end

                otherwise
                    fprintf('Not applying any filtering. \n');
                    dat_filt                    = dat_sel;
            end

            if strcmp(cfg.alignpeak.hilbert.(markername), 'yes')
                fprintf('Converting to Hilbert envelope.\n');
                dat_filt.trial{1}   = envelope(dat_filt.trial{1});
            end

            % segment in trials for alignment
            startsample         = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime*dat_sel.fsample + cfg.alignpeak.toiplot.(markername)(1)*dat_sel.fsample + hdr.nSamplesPre)';
            endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime*dat_sel.fsample + cfg.alignpeak.toiplot.(markername)(2)*dat_sel.fsample + hdr.nSamplesPre)';
            offset              = round(ones(size(startsample))*+cfg.alignpeak.toiplot.(markername)(1)*dat_sel.fsample);
            cfgtemp             = [];
            cfgtemp.trl         = [startsample, endsample, offset];
            cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
            dat_sel_trl         = ft_redefinetrial(cfgtemp,dat_sel);
            clear dat_sel
            dat_filt_trl        = ft_redefinetrial(cfgtemp,dat_filt);
            clear dat_filt

            if strcmp(cfg.alignpeak.demean.(markername), 'yes')
                cfgtemp = [];
                cfgtemp.demean = 'yes';
                cfgtemp.baselinewindow = cfg.alignpeak.toibaseline.(markername);
                dat_sel_trl = ft_preprocessing(cfgtemp, dat_sel_trl);
                dat_filt_trl = ft_preprocessing(cfgtemp, dat_filt_trl);
            end

            if strcmp(cfg.alignpeak.method.(markername), 'crawlback')
                dat_filt40_trl  = ft_redefinetrial(cfgtemp,dat_filt40);
            end

            % peak detection based on baseline and active period
            for itrial = 1 : size(dat_filt_trl.trial,2)

                % time indexes for baseline and active period
                t1_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.alignpeak.toibaseline.(markername)(1), 1, 'first');
                t2_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.alignpeak.toibaseline.(markername)(2), 1, 'last');
                t1_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.alignpeak.toiactive.(markername)(1), 1, 'first');
                t2_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.alignpeak.toiactive.(markername)(2), 1, 'last');

                % find peaks in active period
                [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));

                % find peaks in baseline period
                [peaks_bl{itrial},~]               = findpeaks(dat_filt_trl.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
                if isempty(peaks_bl{itrial})%if no peak, take avg bl value
                    peaks_bl{itrial} = nanmean(dat_filt_trl.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
                end

                % select maximum peak in baseline
                max_peaks_bl(itrial)               = max(peaks_bl{itrial});

            end

            % threshold peaks with baseline, if asked (i.e. if not 0)
            haspeak         = true(1,size(dat_filt_trl.trial,2));
            hasartefact     = false(1,size(dat_filt_trl.trial,2));

            for itrial = 1 : size(dat_filt_trl.trial,2)

                if strcmp(cfg.alignpeak.thresh.method.(markername), 'medianbl')
                    % threshold based on median of max peaks in all baseline periods
                    peaks_sel               = peaks_ac{itrial} > nanmedian(max_peaks_bl) * cfg.alignpeak.thresh.value.(markername);
                    peaks_ac_sel{itrial}    = peaks_ac{itrial}(peaks_sel);
                    locs_ac_sel{itrial}     = locs_ac{itrial}(peaks_sel);
                elseif strcmp(cfg.alignpeak.thresh.method.(markername), 'trial')
                    % threshold based on max peak in trial-by-trial baseline
                    peaks_sel               = peaks_ac{itrial} > max(peaks_bl{itrial}) * cfg.alignpeak.thresh.value.(markername);
                    peaks_ac_sel{itrial}    = peaks_ac{itrial}(peaks_sel);
                    locs_ac_sel{itrial}     = locs_ac{itrial}(peaks_sel);
                elseif strcmp(cfg.alignpeak.thresh.method.(markername), 'both')
                    peaks_sel               = (peaks_ac{itrial} > nanmedian(max_peaks_bl) * cfg.alignpeak.thresh.value.(markername) && peaks_ac{itrial} > max(peaks_bl{itrial}) * cfg.alignpeak.thresh.value.(markername));
                    peaks_ac_sel{itrial}    = peaks_ac{itrial}(peaks_sel);
                    locs_ac_sel{itrial}     = locs_ac{itrial}(peaks_sel);
                else
                    error('wrong cfg.alignpeak.thresh.method parameter');
                end

                if isempty(locs_ac_sel{itrial})
                    haspeak(itrial) = false;
                    fprintf('Could not find peak in trial number %d \n',itrial);
                end
            end

            dat_sel_aligned     = dat_sel_trl;
            dat_filt_aligned    = dat_filt_trl; %for plot
            n_haspeak           = 0;

            if strcmp(cfg.alignpeak.findbegin.(markername),'yes')

                %prepare dummy variables for begin alignment
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    mean_bl(itrial)         = nanmean(dat_filt_trl.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
                end
                med_all_bl = nanmedian(mean_bl);
                %zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %find time cross zero
            end

            for itrial = 1 : size(dat_filt_trl.trial,2)
                if haspeak(itrial)
                    n_haspeak = n_haspeak+1;

                    % find relevant peak according to method
                    if strcmp(cfg.alignpeak.method.(markername),'nearestmax') ||strcmp(cfg.alignpeak.method.(markername),'nearestmin') % Min same as 'max' because timecourse was flipped
                        [~, ip(itrial)] = min(abs(dat_filt_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1))); % peak-time closest to zero

                    elseif strcmp(cfg.alignpeak.method.(markername),'max')
                        [~, ip(itrial)] = max(dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1));

                    elseif strcmp(cfg.alignpeak.method.(markername),'min') % Same as 'max' because timecourse was flipped
                        [~, ip(itrial)] = max(dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1));

                    elseif strcmp(cfg.alignpeak.method.(markername),'firstmax')
                        ip(itrial)      = 1; % first peak within search toi

                    elseif strcmp(cfg.alignpeak.method.(markername),'firstmin')
                        ip(itrial)      = 1; % first peak within search toi

                    elseif strcmp(cfg.alignpeak.method.(markername),'lastmax')
                        [~, indx]       = max(dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1));
                        ip(itrial)      = indx(end); % last peak within search toi

                    elseif strcmp(cfg.alignpeak.method.(markername),'lastmin')
                        [~, indx]       = max(dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1));
                        ip(itrial)      = indx(end); % last peak within search toi

                    elseif strcmp(cfg.alignpeak.method.(markername),'crawlback')
                        % start at max in search window
                        [~, indx1]          = max(dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx(itrial)-1));
                        [~, indx2]          = findpeaks(-dat_filt40_trl.trial{itrial}(1:locs_ac_sel{itrial}(indx1)+t1_ac_indx(itrial)-50));
                        indx2               = indx2(end);
                        ip(itrial)          = 1;
                        locs_ac_sel{itrial} = indx2 - t1_ac_indx(itrial);
                    end
                    timeshift               = dat_filt_trl.time{itrial}(locs_ac_sel{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);

                    %align to the begin of the event if asked
                    if strcmp(cfg.alignpeak.findbegin.(markername),'yes')
                        %find begin according to peak found before
                        peak_loc(itrial)        = locs_ac_sel{itrial}(ip(itrial))+t1_ac_indx(itrial)-1;
                        thresh(itrial)          = med_all_bl + (peaks_ac_sel{itrial}(ip(itrial)) - med_all_bl) * cfg.alignpeak.beginthresh.(markername);

                        x1 = dat_filt_trl.time{itrial}(t1_ac_indx(itrial):peak_loc(itrial));
                        y1 = dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):peak_loc(itrial));
                        x2 = dat_filt_trl.time{itrial}(t1_ac_indx(itrial):peak_loc(itrial));
                        y2 = ones(size(dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):peak_loc(itrial)))) .* thresh(itrial);
                        [x_intersect{itrial}, y_intersect{itrial}] = intersections(x1,y1,x2,y2,true);

                        if ~isempty(x_intersect{itrial})
                            x_intersect{itrial} = x_intersect{itrial}(end);
                            timeshift           = x_intersect{itrial};

                            %remove bad detection (assuming that marker is put after wave begining)
                            if timeshift > 0
                                hasartefact(itrial) = true;
                            end
                        else
                            haspeak(itrial) = false;
                        end
                    end

                    % align time axis to relevant (peak) index
                    dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                    dat_filt_aligned.time{itrial}   = dat_filt_trl.time{itrial} - timeshift; %for plot

                    if abs(timeshift) > cfg.alignpeak.maxtimeshift.(markername)
                        hasartefact(itrial) = true;
                    end

                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).timeshift(itrial) = timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(itrial)  = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(itrial) + timeshift;
                    MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(itrial)     = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(itrial) + seconds(timeshift);
                end
            end

            %remove supplemental trials in case some are ignored
            %because of ~haspeak or hasartefact
            last_envent =  size(dat_filt_trl.trial,2) - sum(~haspeak) - sum(hasartefact);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(1:last_envent);
            MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(1:last_envent);

            %% Plot alignment
            fig             = figure('visible','off');
            fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap

            i_h_temp = 0;
            for ipeak = 1:length(peaks_ac_sel)
                if ~isempty((peaks_ac_sel{ipeak}))
                    i_h_temp = i_h_temp+1;
                    h_temp(i_h_temp) = max(peaks_ac_sel{ipeak});
                end
            end

            h =0;%overdraw of all data

            %Find position of the line to plot
            for itrial = 1 : size(dat_filt_trl.trial,2)
                if haspeak(itrial)
                    line_idx(itrial) = locs_ac_sel{itrial}(ip(itrial))+t1_ac_indx(itrial)-1;
                    if strcmp(cfg.alignpeak.findbegin.(markername), 'yes')
                        line_idx(itrial) = nearest(dat_filt_trl.time{itrial},x_intersect{itrial});
                    end
                end
            end

            subplot(2,2,1);
            hold;
            for itrial = 1 : size(dat_filt_trl.trial,2)
                if haspeak(itrial)
                    color = 'k';
                    t     = dat_filt_trl.time{itrial}(line_idx(itrial));
                    line_position = dat_filt_trl.trial{itrial}(line_idx(itrial));
                    scatter(t,line_position+itrial*h,'bx');
                    %line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                else
                    color = 'r';
                end
                if hasartefact(itrial)
                    color = 'c';
%                     else
%                         color = 'k';
                end
                p = plot(dat_filt_trl.time{itrial},dat_filt_trl.trial{itrial} + itrial*h,'color',color);
                p.ZData = ones(size(p.YData)) .* -1;
            end
            title('Detection');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            ax = axis;
            fill([cfg.alignpeak.toiactive.(markername)(1), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
            fill([cfg.alignpeak.toibaseline.(markername)(1), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
            xlim(cfg.alignpeak.toiplot.(markername));

            subplot(2,2,2);
            hold;
            for itrial = 1 : size(dat_sel_trl.trial,2)
                if haspeak(itrial) && ~hasartefact(itrial)
                    color = 'k';
                    t       = 0;
                    line_position = dat_filt_trl.trial{itrial}(line_idx(itrial));
                    scatter(t,line_position+itrial*h,'bx');
                    %line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
%
%                     else
%                         color = 'r';
%                     end
%                     if hasartefact(itrial)
%                         color = 'c';
%                     else
%                         color = 'k';
%                     end
               p = plot(dat_filt_aligned.time{itrial},dat_filt_aligned.trial{itrial}+ itrial*h,'color',color);
               p.ZData = ones(size(p.YData)) .* -1;
                end
            end
            title('Alignment');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            ax = axis;
            fill([cfg.alignpeak.toiactive.(markername)(1), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
            fill([cfg.alignpeak.toibaseline.(markername)(1), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
            xlim(cfg.alignpeak.toiplot.(markername));

            subplot(2,2,3);
            hold;
            for itrial = 1 : size(dat_sel_trl.trial,2)
                if haspeak(itrial)
                    color = 'k';
                    t       = dat_sel_trl.time{itrial}(line_idx(itrial));
                    line_position = dat_sel_trl.trial{itrial}(line_idx(itrial));
                    scatter(t,line_position+itrial*h,'bx');
                    %line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                else
                    color = 'r';
                end
                if hasartefact(itrial)
                    color = 'c';
%                     else
%                         color = 'k';
                end
                p = plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
                p.ZData = ones(size(p.YData)) .* -1;
            end
            title('Original data');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            ax = axis;
            %fill([cfg.alignpeak.toiactive.(markername)(1), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(2), cfg.alignpeak.toiactive.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
            %fill([cfg.alignpeak.toibaseline.(markername)(1), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(2), cfg.alignpeak.toibaseline.(markername)(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
            xlim(cfg.alignpeak.toiplot.(markername));

            subplot(2,2,4);
            hold;
            for itrial = 1 : size(dat_sel_aligned.trial,2)
                if haspeak(itrial) && ~hasartefact(itrial)
                    color = 'k';
                    t       = 0;
                    line_position = dat_sel_trl.trial{itrial}(line_idx(itrial));
                    scatter(t,line_position+itrial*h,'bx');
                    %line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
%                     else
%                         color = 'r';
%                     end
%                     if hasartefact(itrial)
%                         color = 'c';
%                     else
%                         color = 'k';
%                     end
                p = plot(dat_sel_aligned.time{itrial},dat_sel_aligned.trial{itrial} + itrial*h,'color',color);
                p.ZData = ones(size(p.YData)) .* -1;
                end

            end
            title('Aligned data');
            set(gca, 'YTickLabel', '');
            xlabel('Time (s)');
            axis tight
            xlim(cfg.alignpeak.toiplot.(markername));

            % print to file


            % check if images directory exists, if not create
            if ~isfolder(cfg.imagesavedir)
                ft_notice('creating directory %s', cfg.imagesavedir);
                mkdir(cfg.imagesavedir);
            end

            % check if aligment subdirectory exists, if not create
            isdir_or_mkdir(fullfile(cfg.imagesavedir, 'alignment_peak'));

            % print to file
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
%             print(fig, '-dpdf', fullfile(cfg.imagesavedir, 'alignment_peak',[cfg.prefix,'p',num2str(ipart),'_',convertStringsToChars(markername),'_',dat_sel_aligned.label{1},'_',num2str(idir),'.pdf']));
            print(fig, '-dpng', fullfile(cfg.imagesavedir, 'alignment_peak',[cfg.prefix,'p',num2str(ipart),'_',convertStringsToChars(markername),'_',dat_sel_aligned.label{1},'_',num2str(idir),'.png']),'-r600');


            close all
        end % idir
    end % imarker
end % ipart

if isempty(markername)
    fprintf('cfg.alignpeak.name is empty, no alignment is done\n');
else
    % save data
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned']),'MuseStruct');
end