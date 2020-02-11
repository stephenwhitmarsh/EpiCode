function [MuseStruct] = alignMuseMarkers(cfg, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%Paul : I did'nt understand this part : What if length(varargin) == 3 ? And
%if length(varargin) == 1, where is called MuseStruct ?

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
    markerlist = [];
    for i = 1 : size(cfg.align.name,2)
        if ismember(cfg.align.name{i},cfg.name)
            markerlist = [markerlist, i];
        end
    end
    
    % Go through different parts
    for ipart = 1 : size(cfg.directorylist,2)
        
        for imarker = markerlist
            
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
                % and load header and data
                
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
                    
                    hdr = ft_read_header(dataset);
                    cfgtemp             = [];
                    cfgtemp.dataset     = dataset;
                    dat_sel             = ft_preprocessing(cfgtemp);
                
                elseif isMicromed   
                    hdr = ft_read_header(fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC'])); 
                    dat_sel = preprocessing_eeg_emg(cfg,ipart,idir,true); %notch filter, otherwise the figure at the end is less meaningfull
                    channelnr = 0;
                    channelindex = [];
                    %Find channel of interest
                    for ichannel = 1 : length(dat_sel.label)
                        if strfind(dat_sel.label{ichannel},cfg.align.channel{imarker})
                            fprintf('Found channel with pattern "%s" in %s\n',cfg.align.channel{imarker},MuseStruct{ipart}{idir}.directory);
                            channelnr = channelnr + 1;
                            channelindex = [channelindex ichannel];
                        end
                    end
                    dat_sel.trial{1} = dat_sel.trial{1}(channelindex,:);
                    dat_sel.label = dat_sel.label(channelindex);
                    
                end
                
                
                if channelnr == 0
                    fprintf('Did not find of channel pattern %s!\n',cfg.align.channel{imarker});
                end
                
                if channelnr > 1
                    fprintf('Found more than 1 occurance of channel pattern %s!\n',cfg.align.channel{imarker});
                end
                       

                
                % flip data if required
                if strcmp(cfg.align.method{imarker},'min') || strcmp(cfg.align.method{imarker},'firstmin') || strcmp(cfg.align.method{imarker},'lastmin')
                    dat_sel.trial{1} = -dat_sel.trial{1};
                end
                
                % filtering data of selected channel for peak detection
                switch cfg.align.filter{imarker}
                    
                    case 'lp'
                        fprintf('Lowpass Filtering data for alignment. This can take a while \n');
                        cfgtemp                     = [];
                        cfgtemp.lpfilter            = 'yes';
                        cfgtemp.lpfreq              = cfg.align.freq{imarker};
                        cfgtemp.lpfilttype          = 'fir';
                        dat_filt                    = ft_preprocessing(cfgtemp,dat_sel);
                        
                    case 'bp'
                        fprintf('Bandpass Filtering the data for alignment. This can take a while \n');
                        dat_filt                    = dat_sel;
                        dat_filt.trial{1}           = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.align.freq{imarker}(1),cfg.align.freq{imarker}(2));
                        
                        if strcmp(cfg.align.method{imarker}, 'crawlback')
                            dat_filt40              = dat_sel;
                            dat_filt40.trial{1}     = bandpassFilter(dat_sel.trial{1},dat_sel.fsample,cfg.align.freq{imarker}(1),40);
                        end
                        
                    otherwise
                        fprintf('Not applying any filtering. \n');
                        dat_filt                    = dat_sel;
                end
                if strcmp(cfg.align.hilbert{imarker},'yes')
                    fprintf('Converting to Hilbert envelope.\n');
                    dat_filt.trial{1}   = envelope(dat_filt.trial{1});
                end
                
                % segment in trials for alignment
                startsample         = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime*dat_sel.fsample + cfg.align.toiplot{imarker}(1)*dat_sel.fsample + hdr.nSamplesPre)';          
                endsample           = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime*dat_sel.fsample + cfg.align.toiplot{imarker}(2)*dat_sel.fsample + hdr.nSamplesPre)'; 
                offset              = round(ones(size(startsample))*+cfg.align.toiplot{imarker}(1)*dat_sel.fsample);
                cfgtemp             = [];
                cfgtemp.trl         = [startsample, endsample, offset];
                cfgtemp.trl(:,4)    = 1:size(cfgtemp.trl,1); % try to find trials that are missing aftewards
                dat_sel_trl         = ft_redefinetrial(cfgtemp,dat_sel);
                dat_filt_trl        = ft_redefinetrial(cfgtemp,dat_filt);
                
                if strcmp(cfg.align.method{imarker}, 'crawlback')
                    dat_filt40_trl  = ft_redefinetrial(cfgtemp,dat_filt40);
                end     
                
                % peak detection based on baseline and active period
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    
                    % time indexes for baseline and active period
                    t1_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.align.toibaseline{imarker}(1),1,'first');
                    t2_bl_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.align.toibaseline{imarker}(2),1,'last');
                    t1_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} > cfg.align.toiactive{imarker}(1),1,'first');
                    t2_ac_indx(itrial)                 = find(dat_filt_trl.time{itrial} < cfg.align.toiactive{imarker}(2),1,'last');
                    
                    % find peaks in active period
                    [peaks_ac{itrial},locs_ac{itrial}] = findpeaks(dat_filt_trl.trial{itrial}(t1_ac_indx(itrial):t2_ac_indx(itrial)));
                    
                    % find peaks in baseline period
                    [peaks_bl{itrial},~]               = findpeaks(dat_filt_trl.trial{itrial}(t1_bl_indx(itrial):t2_bl_indx(itrial)));
                    if isempty(peaks_bl{itrial})
                        peaks_bl{itrial} = 0;
                    end
                    
                    % select maximum peak in baseline
                    if ~isempty(peaks_bl{itrial})
                        max_peaks_bl(itrial)           = max(peaks_bl{itrial});
                    else
                        max_peaks_bl(itrial)           = nan;
                    end
                end
                
                % threshold peaks with baseline, if asked (i.e. if not 0)
                haspeak         = true(1,size(dat_filt_trl.trial,2));
                hasartefact     = false(1,size(dat_filt_trl.trial,2));
                
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    
                    % threshold based on median of max peaks in all baseline periods
                    peaks_sel_avg               = peaks_ac{itrial} > nanmean(max_peaks_bl) * cfg.align.thresh(imarker);
                    peaks_ac_sel_avg{itrial}    = peaks_ac{itrial}(peaks_sel_avg);
                    locs_ac_sel_avg{itrial}     = locs_ac{itrial}(peaks_sel_avg);
                    
                    % threshold based on median of max peaks in trial-by-trial baseline
                    peaks_sel_trl               = peaks_ac{itrial} > max(peaks_bl{itrial}) * cfg.align.thresh(imarker);
                    peaks_ac_sel_trl{itrial}    = peaks_ac{itrial}(peaks_sel_trl);
                    locs_ac_sel_trl{itrial}     = locs_ac{itrial}(peaks_sel_trl);
                    
                    if isempty(locs_ac_sel_avg{itrial})
                        haspeak(itrial) = false;
                        fprintf('Could not find peak in trial number %d \n',itrial);
                    end
                end
                
                dat_sel_aligned = dat_sel_trl;
                dat_filt_aligned = dat_filt_trl; %for plot
                
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    if haspeak(itrial)
                        
                        % find relevant peak according to method
                        if strcmp(cfg.align.method{imarker},'nearest')
                            [~, ip(itrial)] = min(abs(dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1))); % peak-time closest to zero
                            
                        elseif strcmp(cfg.align.method{imarker},'max')
                            [~, ip(itrial)] = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1)); 
                            
                        elseif strcmp(cfg.align.method{imarker},'min') % Same as 'max' because timecourse was flipped
                            [~, ip(itrial)] = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1));
                            
                        elseif strcmp(cfg.align.method{imarker},'firstmax')
                            ip(itrial)      = 1; % first peak within search toi
                            
                        elseif strcmp(cfg.align.method{imarker},'firstmin')
                            ip(itrial)      = 1; % first peak within search toi     
                            
                        elseif strcmp(cfg.align.method{imarker},'lastmax')
                            [~, indx]       = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1));   
                            ip(itrial)      = indx(end); % last peak within search toi   
                            
                        elseif strcmp(cfg.align.method{imarker},'lastmin')
                            [~, indx]       = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1));                             
                            ip(itrial)      = indx(end); % last peak within search toi  
                            
                        elseif strcmp(cfg.align.method{imarker},'crawlback')
                            % start at max in search window
                            [~, indx1]              = max(dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}+t1_ac_indx(itrial)-1));       
                            [~, indx2]              = findpeaks(-dat_filt40_trl.trial{itrial}(1:locs_ac_sel_avg{itrial}(indx1)+t1_ac_indx(itrial)-50));
                            indx2                   = indx2(end);
                            ip(itrial)              = 1;
                            locs_ac_sel_avg{itrial} = indx2 - t1_ac_indx(itrial);      
                        end
                        
                        % align time axis to relevant (peak) index 
                        timeshift                       = dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        dat_sel_aligned.time{itrial}    = dat_sel_trl.time{itrial} - timeshift;
                        dat_filt_aligned.time{itrial}   = dat_filt_trl.time{itrial} - timeshift; %for plot
%                         if abs(timeshift) > 0.050 
%                             hasartefact(itrial) = true;
%                         end 
%Paul : this is not adapted to my data, it detects lot of artefacts whereas it is
%not. Could we put the artefact detection of the last version of alignMuseMarkers, which is
%done according to Muse BAD__START__ and BAD__END__ ?

                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).timeshift(itrial)      = timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(itrial) + timeshift;
                        MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial)          = MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).clock(itrial) + seconds(timeshift);
                    end
                end

                %% Plot alignement

                fig             = figure;
                fig.Renderer    = 'Painters'; % Else pdf is saved to bitmap
                
                for itemp = 1:length(peaks_ac_sel_trl)
                    peaks_ac_sel_trl_max(itemp) = max(peaks_ac_sel_trl{itemp});
                end
                h               = mean(peaks_ac_sel_trl_max)/10;%1200; %adapt h because it has to be different between my 3 formats of files. Paul
           
                
                subplot(2,2,1);
                hold;
                for itrial = 1 : size(dat_filt_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t     = dat_filt_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line_position = dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
                    else
                        color = 'k';
                    end
                    plot(dat_filt_trl.time{itrial},dat_filt_trl.trial{itrial} + itrial*h,'color',color);
                end
                title('Peak detection');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                
                subplot(2,2,2);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = 0;
                        line_position = dat_filt_trl.trial{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r'); 
                        
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
                    else
                        color = 'k';
                    end
                   plot(dat_filt_aligned.time{itrial},dat_filt_aligned.trial{itrial}+ itrial*h,'color',color); 
                    
                end
                title('Alignment');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                
               
                subplot(2,2,3);
                hold;
                for itrial = 1 : size(dat_sel_trl.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = dat_sel_trl.time{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line_position = dat_sel_trl.trial{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
                    else
                        color = 'k';
                    end
                    plot(dat_sel_trl.time{itrial},dat_sel_trl.trial{itrial} + itrial*h,'color',color);
                end
                title('Original data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                ax = axis;
                fill([cfg.align.toiactive{imarker}(1), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(2), cfg.align.toiactive{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 0],'edgecolor','none','facealpha',0.1);
                fill([cfg.align.toibaseline{imarker}(1), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(2), cfg.align.toibaseline{imarker}(1)],[ax(3), ax(3),  ax(4),  ax(4)],[0 1 1],'edgecolor','none','facealpha',0.1);
                
                subplot(2,2,4);
                hold;
                for itrial = 1 : size(dat_sel_aligned.trial,2)
                    if haspeak(itrial)
                        color = 'k';
                        t       = 0;
                        line_position = dat_sel_trl.trial{itrial}(locs_ac_sel_avg{itrial}(ip(itrial))+t1_ac_indx(itrial)-1);
                        line([t,t],[(line_position-h)+itrial*h, (line_position+h)+itrial*h],'color','r');
                    else
                        color = 'r';
                    end
                    if hasartefact(itrial)
                        color = 'c';
                    else
                        color = 'k';
                    end
                    plot(dat_sel_aligned.time{itrial},dat_sel_aligned.trial{itrial} + itrial*h,'color',color);
          
                end
                title('Aligned data');
                set(gca, 'YTickLabel', '');
                xlabel('Time (s)');
                axis tight
                
                % print to file
                set(fig,'PaperOrientation','landscape');
                set(fig,'PaperUnits','normalized');
                set(fig,'PaperPosition', [0 0 1 1]);
                print(fig, '-dpdf', fullfile(cfg.imagesavedir,'alignement',[cfg.prefix,'p',num2str(ipart),'alignment_',cfg.name{imarker},'_',dat_sel.label{1},'_',num2str(idir),'.pdf']),'-r600');
                print(fig, '-dpng', fullfile(cfg.imagesavedir,'alignement',[cfg.prefix,'p',num2str(ipart),'alignment_',cfg.name{imarker},'_',dat_sel.label{1},'_',num2str(idir),'.png']),'-r600');
                close all
            end % idir
        end % imarker
    end % ipart
    
    % save data
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned']),'MuseStruct');
    save(fullfile(cfg.datasavedir,[cfg.prefix,'MuseStruct_aligned']),'MuseStruct');
end

