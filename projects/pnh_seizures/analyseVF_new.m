%% set up paths
% requires bandpassFilter.m
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/

% addpath /Volumes/iss01.charpier/stephen.whitmarsh/data/02230_2015-02-25_14-36/m1pNs
% addpath /Volumes/iss01.charpier/stephen.whitmarsh/fieldtrip/

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%% load or create MuseStruct
%% MAKE FUNCTION
if exist('MuseStruct_micro.mat','file')
    load('MuseStruct_micro','MuseStruct_micro');
    load('MuseStruct_macro','MuseStruct_macro');
else
    patient_directory       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/';
    directory_searchstring  = '02230_2015-02-*';
    
    data_searchstring       = '*m1pNs*.ncs';
    MuseStruct_micro        = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);
    save('MuseStruct_micro','MuseStruct_micro');

    data_searchstring       = '*_1pNs*.ncs';
    MuseStruct_macro        = readMuseMarkers(patient_directory,directory_searchstring,data_searchstring);
    save('MuseStruct_macro','MuseStruct_macro');
end

hpfilter                = 'no';
hpfreq                  = 1;
label                   = { 'VF1','VF2','P','PP'};                                                                      % marker name
startend                = { 'VF1__START__','VF1__END__';'RR__START__','RR__END__';'P','P';'PP__START__','PP__END__'};   % start and end Muse marker
prestim                 = [0.5, 0.25,   0.2,    0.25];                                                                                    % list of onset timing with respect to start-marker (s)
poststim                = [2,   1,      0.2,    0.5];  % list of offset timing with respect to end-marker (s)
pad                     = 0.25;
slidestep               = [0.01,0.001,0.001,0.001];
lpfreq                  = [4,4,40,40];                                                                                  % lowpass filter freq to smooth peak detection (Hz)
toiactive               = [-0.1, 0.8;  -0.1,  0.3;  -0.1, 0.1;   0.0,  0.3];                                             % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
toibaseline             = [-1.0,-0.1;  -1.0, -0.1;  -1.0 -0.1;  -1.0,  0.3];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
resamplefs              = 640;                                                                                          % resample data to speed up and reduce memory (Hz), and align micro and macro
channel                 = 'm1pNs_4';                                                                                    % pattern to identify channel on which to based peak detection
alignmethod             = {'first','first','max','first'};                                                                % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
pkthresh                = [0,0.25,0.25,0.25];                                                                           % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period 
datasavedir             = '/network/lustre/iss01/charpier/stephen.whitmarsh/analysis/data';                             % where to write data
imagesavedir            = '/network/lustre/iss01/charpier/stephen.whitmarsh/analysis/images';                           % where to print images

for imarker = 1 : length(label)
    
    %% align Muse markers according to peaks and detect whether they contain artefacts 
    cfg = [];
    cfg.force                = false; % whether to force recalculations
    cfg.markerstart          = startend{imarker,1};     
    cfg.markerend            = startend{imarker,2};         
    cfg.channel              = channel;                
    cfg.lpfreq               = lpfreq(imarker);                       
    cfg.method               = alignmethod{imarker};                    
    cfg.prestim              = prestim(imarker);                        
    cfg.poststim             = poststim(imarker);                        
    cfg.toiactive            = toiactive(imarker,:);              
    cfg.toibaseline          = toibaseline(imarker,:);
    cfg.datasavedir          = datasavedir;      
    cfg.imagesavedir         = imagesavedir;       
    cfg.pkthresh             = pkthresh(imarker);
    [MuseStruct_micro_aligned,MuseStruct_macro_aligned] = alignMuseMarkers(cfg,MuseStruct_micro,MuseStruct_macro);
    
    %% write deadfile for Spyking Circus
 
%     writeSpykingCircusDeadFile(MuseStruct_micro_aligned)
%     writeSpykingCircusDeadFile_concatinated(MuseStruct_micro_aligned)
   
    %% read spike data from Spyking Circus
    
%     cfg                     = [];
%     cfg.force               = false; % whether to force recalculations
%     cfg.startend            = startend(imarker,:);                            
%     cfg.prestim             = prestim(imarker);    
%     cfg.poststim            = poststim(imarker); 
%     cfg.datasavedir         = datasavedir;                     
%     [SpikeRaw, SpikeTrials] = readSpykingCircus(cfg,MuseStruct_micro_aligned);
%     
    %% read EEG data
    
    cfg                     = [];
    cfg.force               = true; % whether to force recalculations
    cfg.startend            = startend(imarker,:);                    
    cfg.label               = label{imarker};     
    cfg.prestim             = prestim(imarker)+pad;                                    
    cfg.poststim            = poststim(imarker)+pad;                     % in seconds
    cfg.hpfilter            = hpfilter;
    cfg.hpfreq              = hpfreq;
    cfg.baseline            = 'no';
    cfg.micro_labels        = {'m1pNs_1','m1pNs_4','m1pNs_6','m1pNs_7','m1pNs_8'};
    cfg.macro_labels        = {'_1pNs_1','_1pNs_2','_1pNs_3','_1pNs_4','_1pNs_5','_1pNs_6','_1pNs_7','_1pNs_8'};
    cfg.resamplefs          = resamplefs;
    cfg.datasavedir         = datasavedir;                     
    [dat_micro, dat_macro]  = readEEG(cfg,MuseStruct_micro_aligned,MuseStruct_macro_aligned);

    %% remove artefacts based on RMS
    chanindx = 0;
    for ichan = 1 : size(dat_micro.label,1)
        if strcmp(channel,dat_micro.label{ichan})
            chanindx = ichan;
        end
    end
    if chanindx == 0
        fprintf('Can not find channel %s! \n',channel);
    end
    datrms = zeros(size(dat_micro.trial,2),1);
    for itrial = 1 : size(dat_micro.trial,2)
        datrms(itrial) = rms(dat_micro.trial{itrial}(chanindx,:));
    end
    fprintf('removing %d artefacted trials out of %d (> 3 STD(RMS)) \n',sum((datrms-mean(datrms)) >= std(datrms)*3),length(datrms));
    
    cfg = [];
    cfg.trials = find(datrms-mean(datrms) < std(datrms)*3);
    dat_micro_clean = ft_selectdata(cfg,dat_micro);
    dat_macro_clean = ft_selectdata(cfg,dat_macro);
    
    cfg = [];
    cfg.vartrllength = 2;
    avg = ft_timelockanalysis(cfg,dat_micro_clean);
    clear c
    for itrial = 1 : size(dat_micro_clean.trial,2)
        fprintf('correlating trial %d out of %d \n',itrial,size(dat_micro_clean.trial,2));
        s = min(size(avg.avg(chanindx,:),2),size(dat_micro_clean.trial{itrial}(chanindx,:),2));
        s = find(avg.time > prestim(imarker),1,'first');
        c(itrial) = corr(avg.avg(chanindx,1:s)',dat_micro_clean.trial{itrial}(chanindx,1:s)');
    end
    
    cfg = [];
    cfg.trials = find(c > 0);
    fprintf('removing %d trials out of %d (c <= 0) \n',sum(c<=0),size(c,2));
    dat_micro_clean = ft_selectdata(cfg,dat_micro_clean);
    dat_macro_clean = ft_selectdata(cfg,dat_macro_clean);
    c = c(c>0);
    
    %% plot EEG data
    
    cfg                     = [];
    cfg.force               = true; % whether to force recalculations
    cfg.label               = label{imarker};
    cfg.prestim             = prestim(imarker); 
    cfg.poststim            = poststim(imarker);
    cfg.slidestep           = slidestep(imarker);
    [Y, I]                  = sort(c,'descend');
    cfg.representtrials     = I(1:20);
    cfg.resamplefs          = resamplefs;
    cfg.channel             = channel;
    cfg.binsize             = 0.1;
    cfg.datasavedir         = datasavedir;                     
    cfg.imagesavedir        = imagesavedir;
%     plotdata(cfg,dat_micro_clean,dat_macro_clean,SpikeTrials)
    plotdata_corr(cfg,dat_micro_clean,dat_macro_clean,[])

end



