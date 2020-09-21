function [MuseStruct] = alignMuseMarkers(cfg, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supposer que le pic est vers le milieu de l'événement
% function [MuseStruct] = alignMuseMarkers(cfg, MuseStruct, force)
%
% Automatically determines timing shift of events read from MUSE, according
% to different alignment methods. Results are saved in a 'timeshift' field 
% for every marker event in the MuseStruct. 
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
% cfg.align.name                = e.g.: {'Hspike','SpikeDetect'};                                % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% cfg.align.channel             = e.g.: {'_HaT2_1','_HaT2_1'};                                   % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% cfg.align.abs                 = e.g.: {'no','no'};                                             % whether to determine peaks on absolute (both max and min)         
% cfg.align.method              = e.g.: {'crawlback','min'};                                     % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% cfg.align.filter              = e.g.: {'bp','bp'};                                             % what filter to use (bp,hp,lp)
% cfg.align.freq                = e.g.: {[1, 10],[1, 40]};                                       % lowpass filter freq to smooth peak detection (Hz)
% cfg.align.hilbert             = e.g.: {'no','no'};                                             % whether to first create a hilbert transform (to determine peak of frequency response rather than timecourse
% cfg.align.thresh              = e.g.: [1, 0];                                                  % only peaks that are more than x percent above threshold, 0 to disable
% cfg.align.toiplot{1}          = e.g.: [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.align.toiactive{1}        = e.g.: [-0.05,  0.150];                                         % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% cfg.align.toibaseline{1}      = e.g.: [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.align.toiplot{2}          = e.g.: [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% cfg.align.toiactive{2}        = e.g.: [-0.05,  0.05];                                          % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% cfg.align.toibaseline{2}      = e.g.: [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if results exist
fname = fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned_begin.mat']);

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
    
    %get format to adapt script for each format
    %specificities :
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);
    

    % select those markers to align
%     markerlist = [];
%     for i = 1 : size(cfg.align.name,2)
%         if ismember(cfg.align.name{i},cfg.name)
%             markerlist = [markerlist, i];
%         end
%     end
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        

        for imarker = 1 : size(cfg.align.name,2)%markerlist

            % find data directories that have the required event
            markerindx = [];
            for idir = 1 : size(MuseStruct{ipart},2)
                if isfield(MuseStruct{ipart}{idir},'markers')
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{imarker,1})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')                        
                            if ~isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                                fprintf('Found %d event timings for %s \n',size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock,2),cfg.muse.startend{imarker,1});
                                markerindx = [markerindx idir];
                            end
                        end
                    end
                end
            end
            if isempty(markerindx)
                fprintf('Did not find any events for marker %s\n',cfg.name{imarker});
            end
            
            % loop over all data directories that have the marker
            for idir = markerindx
                
                % find datafilename corresponding to requested electrode 
                if isNeuralynx                    
                    filelist = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},'*.ncs'));
                    channelnr = 0;
                    for ifile = 1 : size(filelist,1)
                        if strfind(filelist(ifile).name,cfg.align.channel{imarker})
                            fprintf('Found channel with pattern "%s" in %s\n',cfg.align.channel{imarker},filelist(ifile).name);
                            dataset = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},filelist(ifile).name);
                            channelnr = channelnr + 1;
                        end
                    end
                    if channelnr == 0
                        fprintf('Did not find of channel pattern %s!\n',cfg.align.channel{imarker});
                    end
                    if channelnr > 1
                        fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.align.channel{imarker});
                    end
                    
                elseif isMicromed
                    dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
                    
                elseif isBrainvision
                    dataset = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
                    
                end
                
                
                           
                %load data and header
                hdr = ft_read_header(dataset);
                
                
                cfgtemp             = [];
                
                if isfield(cfg.align, 'reref')
                    if strcmp(cfg.align.reref,'yes')
                        cfgtemp.reref       = 'yes';
                        cfgtemp.rerefmethod = 'avg';
                        cfgtemp.refchannel  = cfg.labels.macro';
                    end
                end
                
                if isfield(cfg.align, 'notch')
                    if strcmp(cfg.align.notch,'yes')
                        cfgtemp.bsfilter       = 'yes';
                        cfgtemp.bsfreq         = [49 51];
                    end
                end
                
                cfgtemp.dataset     = dataset;
                dat_sel             = ft_preprocessing(cfgtemp);
                
                if isMicromed || isBrainvision %choose channel
                    cfgtemp         = [];
                    cfgtemp.channel = cfg.align.channel{imarker};
                    dat_sel         = ft_selectdata(cfgtemp,dat_sel);
                    if size(dat_sel.trial{1},1) == 0
                        fprintf('Did not find of channel pattern %s!\n',cfg.align.channel{imarker});
                    end
                    if size(dat_sel.trial{1},1) > 1
                        fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.align.channel{imarker});
                    end
                end
                
                %FLIP DATA
                if strcmp(cfg.align.flip, 'yes')
                    dat_sel.trial{1} = -dat_sel.trial{1};
                end
                
%                 % flip data if required
%                 if strcmp(cfg.align.method{imarker},'min') || strcmp(cfg.align.method{imarker},'firstmin') || strcmp(cfg.align.method{imarker},'lastmin') || strcmp(cfg.align.method{imarker},'nearestmin')
%                     dat_sel.trial{1} = -dat_sel.trial{1};
%                 end
                
%                 % filtering data of selected channel for peak detection
%                 switch cfg.align.filter{imarker}
%                     
%                     case 'lp'
                        fprintf('Lowpass Filtering data for alignment. This can take a while \n');
                        cfgtemp                     = [];
                        cfgtemp.lpfilter            = 'yes';
                        cfgtemp.lpfreq              = 10;%cfg.align.freq{imarker};
                        cfgtemp.lpfilttype          = 'fir';
                        dat_filt                    = ft_preprocessing(cfgtemp,dat_sel);
%                         
%                     case 'bp'
%                         fprintf('Bandpass Filtering the data for alignment. This can take a while \n');
%                         dat_filt                    = dat_sel;
%                         dat_filt.trial{1}           = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.align.freq{imarker}(1),cfg.align.freq{imarker}(2));
%                         
%                         if strcmp(cfg.align.method{imarker}, 'crawlback')
%                             dat_filt40              = dat_sel;
%                             dat_filt40.trial{1}     = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.align.freq{imarker}(1),40);
%                         end
%                         
%                     otherwise
%                         fprintf('Not applying any filtering. \n');
%                         dat_filt                    = dat_sel;
%                 end
%                 if strcmp(cfg.align.hilbert{imarker},'yes')
%                     fprintf('Converting to Hilbert envelope.\n');
%                     dat_filt.trial{1}   = envelope(dat_filt.trial{1});
%                 end
                
                % segment in trials for alignment
                startsample         = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime*dat_sel.fsample + cfg.align.toiplot{imarker}(1)*dat_sel.fsample + hdr.nSamplesPre)';          
                endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime*dat_sel.fsample + cfg.align.toiplot{imarker}(2)*dat_sel.fsample + hdr.nSamplesPre)'; 
                offset              = round(ones(size(startsample))*+cfg.align.toiplot{imarker}(1)*dat_sel.fsample);
                cfgtemp             = [];
                cfgtemp.trl         = [startsample, endsample, offset];
                cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                dat_sel_trl         = ft_redefinetrial(cfgtemp,dat_sel);
                   
%                 dat_filt_trl = dat_sel_trl;
                endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime*dat_sel.fsample + hdr.nSamplesPre)'; 
                cfgtemp             = [];
                cfgtemp.trl         = [startsample, endsample, offset];
                cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                dat_filt_trl        = ft_redefinetrial(cfgtemp,dat_filt);
                clear dat_sel dat_filt
%                 
%                 if strcmp(cfg.align.method{imarker}, 'crawlback')
%                     dat_filt40_trl  = ft_redefinetrial(cfgtemp,dat_filt40);
%                 end     
                
                % peak detection based on baseline and active period
                
%                          
%                 test = dat_filt_trl;
%                 dat_filt_trl = test;
%                 dat_filt_trl = dat_sel_trl; 
%                 test = dat_sel_trl; 
%                 plot(dat_filt_trl.time{itrial}, dat_filt_trl.trial{itrial});
%                 plot(test.time{itrial}, test.trial{itrial});
                
                %set minimum as zero so the integral is only increasing values
                for itrial =  1: size(dat_filt_trl.trial,2)
                    dat_filt_trl.trial{itrial} = dat_filt_trl.trial{itrial} - min(dat_filt_trl.trial{itrial});
                end
%                 
%                 cfgtemp = [];
%                 cfgtemp.lpfilter = 'yes';
%                 cfgtemp.lpfreq = 10;
%                 dat_filt_trl = ft_preprocessing(cfgtemp,dat_filt_trl);
% %                    baseline correction (so integral is ~ zero on baseline)
%                 cfgtemp = [];
%                 cfgtemp.demean = 'yes';
%                 cfgtemp.baselinewindow = cfg.align.toibaseline{imarker};
%                 dat_filt_trl = ft_preprocessing(cfgtemp, dat_filt_trl);
%                 
                % compute integral
                for itrial = 1: size(dat_filt_trl.trial,2)
                    dat_filt_trl.trial{itrial} = cumtrapz(dat_filt_trl.trial{itrial});
                end
%                 
                %remove linear trend
                cfgtemp = [];
                cfgtemp.detrend = 'yes';
                dat_filt_trl = ft_preprocessing(cfgtemp,dat_filt_trl);
                
                %flip for peak detection
                for itrial = 1: size(dat_filt_trl.trial,2)
                    dat_filt_trl.trial{itrial} =-(dat_filt_trl.trial{itrial});
                end
              
                 %initialize variables for alignment
                dat_sel_aligned     = dat_sel_trl;
                dat_filt_aligned    = dat_filt_trl; %for plot
                haspeak         = true(1,size(dat_filt_trl.trial,2));
                hasartefact     = false(1,size(dat_filt_trl.trial,2));
                keeptrial = 0;
                
                % time indexes for baseline and active period
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    t1_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.align.toibaseline{imarker}(1),1,'first');
                    t2_bl_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.align.toibaseline{imarker}(2),1,'last');
                    t1_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} > cfg.align.toiactive{imarker}(1),1,'first');
                    t2_ac_indx{itrial}                 = find(dat_filt_trl.time{itrial} < cfg.align.toiactive{imarker}(2),1,'last');
                end
                
                %find max peak of dat_filt
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    
                    % find peaks in active period
                    [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}:t2_ac_indx{itrial}));
                    
                    % find peaks in baseline period
                    [peaks_bl{itrial},~]               = findpeaks(dat_filt_trl.trial{itrial}(t1_bl_indx{itrial}:t2_bl_indx{itrial}));
                    if isempty(peaks_bl{itrial})%if no peak, take avg bl value
                        peaks_bl{itrial} = nanmean(dat_filt_trl.trial{itrial}(t1_bl_indx{itrial}:t2_bl_indx{itrial}));
                    end
                    
                    % select maximum peak in baseline
                    max_peaks_bl{itrial}           = max(peaks_bl{itrial});
                end
                
                %select peak
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    %remove peaks lower than baseline
                    peaks_sel               = peaks_ac{itrial} > max(peaks_bl{itrial});
                    peaks_ac_sel{itrial}    = peaks_ac{itrial}(peaks_sel);
                    locs_ac_sel{itrial}     = locs_ac{itrial}(peaks_sel);
                    
%                     
                    %choose maximum of selected peak
                    if ~isempty(peaks_ac_sel{itrial})
                        peaks_ac_sel{itrial}    = peaks_ac{itrial}(end); %last peak superior of bl
                        locs_ac_sel{itrial}     = locs_ac{itrial}(end);
                        [peaks_ac_sel{itrial}, loc_idx]    = max(peaks_ac_sel{itrial});
                        locs_ac_sel{itrial}                = locs_ac_sel{itrial}(loc_idx);
%                         peaks_ac_sel{itrial}               = peaks_ac_sel{itrial};(1)
%                         locs_ac_sel{itrial}                = locs_ac_sel{itrial}(1);
                        timeshift               = dat_filt_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    else
                        locs_ac_sel{itrial}     = [];
                        haspeak(itrial) = false;
                        fprintf('Could not find peak in trial number %d \n',itrial);
                        timeshift = 0;
                    end
                    
                    % align data
                    dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                    dat_filt_aligned.time{itrial}   = dat_filt_trl.time{itrial} - timeshift; %for plot
                    
                    % if no detection in a certain window, set as artefact
                    % because marker put in peak of event : timeshift has
                    % to be negative (begin of event before peak)
                    if abs(timeshift) > cfg.align.maxtimeshift || timeshift > -0.03
                        hasartefact(itrial) = true;
                    end
                    
                    %keep only trials with good detection
                    if ~hasartefact(itrial) && haspeak(itrial)
                        keeptrial = keeptrial +1;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial) + timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial)          = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial) + seconds(timeshift);
                    end
                end
                
                %remove supplemental trials in case some are ignored
                %because of hasartefact or ~haspeak
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(1:keeptrial);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(1:keeptrial);
                
                %% Plot alignement
                
                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    h_temp(itrial) = max(dat_filt_trl.trial{itrial}(t1_ac_indx{itrial}: t2_ac_indx{itrial})); %amplitude of peak vs baseline
                end
                h               = mean(h_temp)/2;
                
                
                subplot(2,2,1);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    color = 'k';
                    t     = dat_filt_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line_position =   dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    if hasartefact(itrial)
                        color = 'c';
                    end
                    if ~haspeak(itrial)
                        color = 'r';
                    end
                    plot(dat_filt_trl.time{itrial},dat_filt_trl.trial{itrial} + itrial*h,'color',color);
                    xlim([-0.5 0.5])
                end
                title('Detection : integral of trial + detrend + flip (find max peak)');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.align.toiplot{imarker});
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 0 1],'edgecolor','none','facealpha',0.1);
                
                subplot(2,2,2);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    color = 'k';
                    t       = 0;
                    line_position = dat_filt_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    if isempty(line_position); t=[]; end
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    if ~hasartefact(itrial) && haspeak(itrial)
                        plot(dat_filt_aligned.time{itrial},dat_filt_aligned.trial{itrial}+ itrial*h,'color',color);
                    end
                end
                title('Alignment');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.align.toiplot{imarker});
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 0 1],'edgecolor','none','facealpha',0.1);
                
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    h_temp(itrial) = max(dat_sel_trl.trial{itrial}(t1_ac_indx{itrial}: t2_ac_indx{itrial})); %amplitude of peak vs baseline
                end
                h               = mean(h_temp)/2;
                
                subplot(2,2,3);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    color = 'k';
                    t       = dat_sel_trl.time{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line_position = dat_sel_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    
                    if hasartefact(itrial)
                        color = 'c';
                    end
                    if ~haspeak(itrial)
                        color = 'r';
                    end
                    plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
                end
                title('Original data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.align.toiplot{imarker});
                
                subplot(2,2,4);
                hold;
                for itrial = 1 : size(dat_sel_aligned.trial,2)
                    color = 'k';
                    t       = 0;
                    line_position = dat_sel_trl.trial{itrial}(locs_ac_sel{itrial}+t1_ac_indx{itrial}-1);
                    if isempty(line_position); t=[]; end
                    line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    if ~hasartefact(itrial) && haspeak(itrial)
                        plot(dat_sel_aligned.time{itrial},dat_sel_aligned.trial{itrial} + itrial*h,'color',color);
                    end
                    
                end
                title('Aligned data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                xlim(cfg.align.toiplot{imarker});
                
                % print to file
                    
                   
                % check if images directory exists, if not create
                if ~isfolder(cfg.imagesavedir)
                    ft_notice('creating directory %s', cfg.imagesavedir);
                    mkdir(cfg.imagesavedir);
                end
                
                % check if aligment subdirectory exists, if not create
                if ~isfolder(fullfile(cfg.imagesavedir,'alignment_begin'))
                    ft_notice('creating directory %s', fullfile(cfg.imagesavedir,'alignment_begin'));
                    mkdir(fullfile(cfg.imagesavedir,'alignment_begin'));
                end
                
                % print to file    
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'alignment_begin',[cfg.prefix,'p',num2str(ipart),'_',cfg.name{imarker},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.pdf']));
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'alignment_begin',[cfg.prefix,'p',num2str(ipart),'_',cfg.name{imarker},'_',dat_sel_aligned.label{1},'_',num2str(idir),'.png']),'-r600');


                close all
            end % idir
        end % imarker
    end % ipart
    
    % save data
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned']),'MuseStruct');
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned']),'MuseStruct');
end

