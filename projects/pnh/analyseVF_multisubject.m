%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/mlib6/

%
% addpath /Volumes/iss01.charpier/stephen.whitmarsh/data/02230_2015-02-25_14-36/m1pNs
% addpath /Volumes/iss01.charpier/stephen.whitmarsh/fieldtrip/

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
maxNumCompThreads(8)

% Setting parameters

% Patient 1, perivetricular heterotopia #1
hpfilter{1}                = {'no','no','no','no'};
hpfreq{1}                  = {1,1,1,1};
bpfilter{1}                = {'no','no','no','no'};
bpfreq{1}                  = {1,1,1,1};
hilbertenv{1}              = {'no','no','no','no'};
label{1}                   = { 'VF1','RR','P','PP'};                                                                      % marker name
startend{1}                = { 'VF1__START__','VF1__END__';'RR__START__','RR__END__';'P','P';'PP__START__','PP__END__'};   % start and end Muse marker
prestim{1}                 = [0.50,  0.50,  0.2,  0.25];                                                                  % list of onset timing with respect to start-marker (s)
poststim{1}                = [2.00,  0.50,  0.2,  0.50];                                                                   % list of offset timing with respect to end-marker (s)
pad{1}                     = 0.25;
slidestep{1}               = [0.01,0.001,0.001,0.001];
lpfreq{1}                  = [4,4,40,40];                                                                                  % lowpass filter freq to smooth peak detection (Hz)
toiactive{1}               = [-0.1, 0.8;  -0.1,  0.3;  -0.1, 0.1;   0.0,  0.3];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
toibaseline{1}             = [-1.0,-0.1;  -1.0, -0.1;  -1.0 -0.1;  -1.0,  0.3];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
resamplefs{1}              = 640;                                                                                          % resample data to speed up and reduce memory (Hz), and align micro and macro
channel{1}                 = {'m1pNs_4','m1pNs_4','m1pNs_4','m1pNs_4'};                                                                                    % pattern to identify channel on which to based peak detection
alignmethod{1}             = {'first','first','max','first'};                                                              % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
pkthresh{1}                = [0,0.25,0.25,0.25];                                                                           % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
datasavedir{1}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/data';         % where to write data
imagesavedir{1}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/images';       % where to print images
patient_directory{1}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg';
directory_searchstring{1}  = '02230_2015-02-*';
micro_searchstring{1}{1}   = '*m1pNs*.ncs';
macro_searchstring{1}{1}   = '*_1pNs*.ncs';
micro_labels{1,1}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,1}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,2}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,2}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,3}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,3}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,4}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,4}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
timwin{1,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,2}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,3}                = [-0.0025 0.0025];  % for plotting spikerate
timwin{1,4}                = [-0.01 0.01];      % for plotting spikerate

% Patient 2, microgyria
hpfilter{2}                = {'no','no','no','no'};
hpfreq{2}                  = {1,1,1,1};
bpfilter{2}                = {'no','no','no','no'};
bpfreq{2}                  = {1,1,1,1};
hilbertenv{2}              = {'no','no','no','no'};
label{2}                   = { 'RR','other'};
startend{2}                = { 'RR__START__','RR__END__';'other__START__','other__END__'};
prestim{2}                 = [0.5,0.5];
poststim{2}                = [2,2];
pad{2}                     = 0.25;
slidestep{2}               = [0.01,0.01];
lpfreq{2}                  = [4,4];
toiactive{2}               = [-0.1, 0.8; -0.1, 0.8];
toibaseline{2}             = [-1.0,-0.1; -1.0,-0.1];
resamplefs{2}              = 640;
channel{2}                 = {'mPlan_2','mPlan_2','mPlan_2','mPlan_2'};
alignmethod{2}             = {'max','max'};
pkthresh{2}                = [0.25, 0.25];
datasavedir{2}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/analysis/data';
imagesavedir{2}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/analysis/images';
patient_directory{2}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/eeg';
directory_searchstring{2}  = '02476_2017-04-*';
micro_searchstring{2}{1}   = '*mPlan*.ncs';
macro_searchstring{2}{1}   = '*_Plan*.ncs';
micro_labels{2,1}          = ["mPlan_1","mPlan_2","mPlan_4","mPlan_5","mPlan_6","mPlan_7","mPlan_8"];
macro_labels{2,1}          = ["_Plan_1","_Plan_2","_Plan_3","_Plan_4","_Plan_5","_Plan_6","_Plan_7","_Plan_8"];
micro_labels{2,2}          = ["mPlan_1","mPlan_2","mPlan_4","mPlan_5","mPlan_6","mPlan_7","mPlan_8"];
macro_labels{2,2}          = ["_Plan_1","_Plan_2","_Plan_3","_Plan_4","_Plan_5","_Plan_6","_Plan_7","_Plan_8"];
timwin{2,1}                = [-0.01 0.01]; % for plotting spikerate
timwin{2,2}                = [-0.01 0.01]; % for plotting spikerate

% Patient 3A, periventricular heterotopia #2
hpfilter{3}                = {'no','yes','no','no'};
hpfreq{3}                  = {1,1,1,1};
bpfilter{3}                = {'no','no','no','no'};
bpfreq{3}                  = {1,1,1,1};
hilbertenv{3}              = {'no','no','no','no'};
label{3}                   = { 'H02614_01','H02614_02','H02614_03','H02614_04'};
startend{3}                = { 'H02614_01__START__','H02614_01__END__';'H02614_02','H02614_02';'H02614_03__START__','H02614_03__END__';'H02614_04__START__','H02614_04__END__'};
prestim{3}                 = [2, 0.2, 0.5, 0.5];
poststim{3}                = [2, 0.2, 5,   2.5];
pad{3}                     = 0.25;
slidestep{3}               = [0.01,0.001,0.01,0.01];
lpfreq{3}                  = [4,40,4,40];
toiactive{3}               = [-0.1, 0.8; -0.1, 0.1; -0.1, 0.8; -0.1, 0.2];
toibaseline{3}             = [-1.0,-0.1; -1.0,-0.1; -1.0,-0.1; -1.0,-0.2];
resamplefs{3}              = 640;
channel{3}                 = {'mCasd_2','mCasd_2','mCasd_2','mCasd_2'};
alignmethod{3}             = {'first','max','first','first'};
pkthresh{3}                = [0.25, 0.25, 0.25, 0.25];
datasavedir{3}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';
imagesavedir{3}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';
patient_directory{3}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg';
directory_searchstring{3}  = '02614_2018-06-*';
micro_searchstring{3}{1}   = '*_m*.ncs'; %'*mCasd*.ncs';
macro_searchstring{3}{1}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*_m*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*_m*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*_m*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_labels{3,1}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7"];
macro_labels{3,1}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,2}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7"];
macro_labels{3,2}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,3}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7"];
macro_labels{3,3}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,4}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7"];
macro_labels{3,4}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];

% micro_labels{3,1}          = ["mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
% macro_labels{3,1}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
% micro_labels{3,2}          = ["mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8","mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
% macro_labels{3,2}          = ["_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
% micro_labels{3,3}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8","mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
% macro_labels{3,3}          = ["_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
% micro_labels{3,4}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8","mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
% macro_labels{3,4}          = ["_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];

timwin{3,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{3,2}                = [-0.01 0.01];      % for plotting spikerate
timwin{3,3}                = [-0.1 0.1];  % for plotting spikerate
timwin{3,4}                = [-0.1 0.1];      % for plotting spikerate


% Patient 3B, periventricular heterotopia #2
hpfilter{4}                = {'no','no','no'};
hpfreq{4}                  = {1,1,1,1};
bpfilter{4}                = {'yes','yes','no'};
bpfreq{4}                  = {[50, 150],[5 15],1,1};
hilbertenv{4}              = {'yes','no','no'};
label{4}                   = { 'H02614_07','H02614_08','H02614_09'};
startend{4}                = { 'H02614_07__START__','H02614_07__END__';'H02614_08__START__','H02614_08__END__';'H02614_09','H02614_09'}; % 9=spike
prestim{4}                 = [0.5, 1, .2];
poststim{4}                = [0.5, 2, .2];
pad{4}                     = 0.25;
slidestep{4}               = [0.01,0.01,0.001];
lpfreq{4}                  = [4,4,40];
toiactive{4}               = [-0.1, 0.8; -0.1, 0.5; -0.1, 0.1];
toibaseline{4}             = [-1.0,-0.1; -1.0,-0.1; -1.0,-0.1];
resamplefs{4}              = 640;
channel{4}                 = {'mTNmi_3','mTNmi_2','mTNmi_3'};
alignmethod{4}             = {'max','first','max'};
pkthresh{4}                = [0.25, 0.25, 0.25];
datasavedir{4}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';
imagesavedir{4}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';
patient_directory{4}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg';
directory_searchstring{4}  = '02614_2018-06-*';
micro_searchstring{4}{1}   = '*mTNmi*.ncs';
macro_searchstring{4}{1}   = '*_TNmi*.ncs';
micro_searchstring{4}{2}   = '*mTNmi*.ncs';
macro_searchstring{4}{2}   = '*_TNmi*.ncs';
micro_searchstring{4}{3}   = '*mTNmi*.ncs';
macro_searchstring{4}{3}   = '*_TNmi*.ncs';
micro_labels{4,1}          = ["mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
macro_labels{4,1}          = ["_TNmi_1","_TNmi_2","_TNmi_3","_TNmi_4","_TNmi_5","_TNmi_6","_TNmi_7","_TNmi_8"];
micro_labels{4,2}          = ["mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
macro_labels{4,2}          = ["_TNmi_1","_TNmi_2","_TNmi_3","_TNmi_4","_TNmi_5","_TNmi_6","_TNmi_7","_TNmi_8"];
micro_labels{4,3}          = ["mTNmi_1","mTNmi_2","mTNmi_3","mTNmi_4","mTNmi_5","mTNmi_6","mTNmi_7","mTNmi_8"];
macro_labels{4,3}          = ["_TNmi_1","_TNmi_2","_TNmi_3","_TNmi_4","_TNmi_5","_TNmi_6","_TNmi_7","_TNmi_8"];
timwin{4,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{4,2}                = [-0.1 0.1];      % for plotting spikerate
timwin{4,3}                = [-0.01 0.01];  % for plotting spikerate



ipatient = 1;
imarker = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis project Valerio Frazzini %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ipatient = 1 : 1
    
    % read muse markers
    cfg = [];
    cfg.force                   = true;
    cfg.patient_directory       = patient_directory{ipatient};
    cfg.directory_searchstring  = directory_searchstring{ipatient};
    cfg.channel_searchstring    = micro_searchstring{ipatient}{1};
    cfg.datasavedir             = datasavedir{ipatient};
    cfg.suffix                  = [num2str(ipatient),'_micro'];
    MuseStruct_micro            = readMuseMarkers(cfg);
    
    cfg.suffix                  = [num2str(ipatient),'_macro'];
    cfg.channel_searchstring    = macro_searchstring{ipatient}{1};
    MuseStruct_macro            = readMuseMarkers(cfg);
    
    % loop over all markers. Last saves all markers in musestruct
    for imarker = 1 : length(label{ipatient})
        
        % align Muse markers according to peaks and detect whether they contain artefacts
        cfg = [];
        cfg.prefix              = [num2str(ipatient),'-'];
        cfg.force               = false; % whether to force recalculations
        cfg.markerstart         = startend{ipatient}{imarker,1};
        cfg.markerend           = startend{ipatient}{imarker,2};
        cfg.channel             = channel{ipatient}{imarker};
        cfg.lpfreq              = lpfreq{ipatient}(imarker);
        cfg.method              = alignmethod{ipatient}{imarker};
        cfg.prestim             = prestim{ipatient}(imarker);
        cfg.poststim            = poststim{ipatient}(imarker);
        cfg.toiactive           = toiactive{ipatient}(imarker,:);
        cfg.toibaseline         = toibaseline{ipatient}(imarker,:);
        cfg.datasavedir         = datasavedir{ipatient};
        cfg.imagesavedir        = imagesavedir{ipatient};
        cfg.pkthresh            = pkthresh{ipatient}(imarker);
        cfg.bpfilter            = bpfilter{ipatient}{imarker};
        cfg.bpfreq              = bpfreq{ipatient}{imarker};
        cfg.hilbert             = hilbertenv{ipatient}{imarker};
        [MuseStruct_micro,MuseStruct_macro] = alignMuseMarkers(cfg,MuseStruct_micro,MuseStruct_macro);
        
    end
    
    % write all data as one big file for SC
    cfg                     = [];
    cfg.prefix              = [num2str(ipatient),'-'];   
    cfg.force               = false;
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.channel             = micro_labels{ipatient,1}; % inconsistent
    [sampleinfo,fnames_ncs,deadfile_ms,deadfile_samples] = writemicroforspykingcircus_alldata(cfg, MuseStruct_micro);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOW RUN SPYKING-CIRCUS ON *-all_data_*.ncs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % read and plot all data as one big file from SC
    cfg                     = [];
    cfg.prefix              = [num2str(ipatient),'-'];      
    cfg.force               = false;
    cfg.startend            = startend{ipatient};
    cfg.prestim             = prestim{ipatient}; % +pad{ipatient};
    cfg.poststim            = poststim{ipatient}; % +pad{ipatient};
    cfg.imagesavedir        = imagesavedir{ipatient};
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.label               = label{ipatient};
    cfg.channel             = channel{ipatient}{imarker};
    cfg.resamplefs          = resamplefs{ipatient};
    cfg.timwin              = timwin(ipatient,:);
    cfg.suffix              = '-1';
    cfg.sampleinfo          = sampleinfo;
    [SpikeRaw, SpikeTrials] = readSpykingCircus_allmarkers_alldata(cfg,MuseStruct_micro);
    
    
    for imarker = 1 : length(label{ipatient})
        
        %% read LFP data
        cfg                     = [];
        cfg.force               = true; % whether to force recalculations
        cfg.startend            = startend{ipatient}(imarker,:);
        cfg.label               = label{ipatient}{imarker};
        cfg.prestim             = prestim{ipatient}(imarker)  + pad{ipatient};
        cfg.poststim            = poststim{ipatient}(imarker) + pad{ipatient};
        cfg.hpfilter            = hpfilter{ipatient}{imarker};
        cfg.hpfreq              = hpfreq{ipatient}{imarker};
        cfg.baseline            = 'no';
        cfg.micro_labels        = micro_labels{ipatient,imarker};
        cfg.macro_labels        = macro_labels{ipatient,imarker};
        cfg.resamplefs          = resamplefs{ipatient};
        cfg.datasavedir         = datasavedir{ipatient};
        [dat_micro, dat_macro]  = readEEG(cfg,MuseStruct_micro,MuseStruct_macro); % rename to readLFP
        
        %% remove artefacts based on RMS - make function
        cfg = [];
        cfg.prestim             = prestim{ipatient}(imarker); % +pad{ipatient};
        cfg.channel             = channel{ipatient}{imarker};
        [dat_micro,dat_macro,c] = removeartefacts_corr(cfg,dat_micro, dat_macro); % remove c, an d place in plotLFP function
        [~,~,c] = removeartefacts_corr(cfg,dat_micro, dat_macro); % remove c, an d place in plotLFP function
        
        %% plot EEG data
        cfg                     = [];
        cfg.force               = true;
        cfg.label               = label{ipatient}{imarker};
        cfg.prestim             = prestim{ipatient}(imarker);
        cfg.poststim            = poststim{ipatient}(imarker);
        cfg.slidestep           = slidestep{ipatient}(imarker);
        [Y, I]                  = sort(c,'descend');
        cfg.representtrials     = I(1:20); % do within function
        cfg.resamplefs          = resamplefs{ipatient};
        cfg.channel             = channel{ipatient}{imarker};
        cfg.binsize             = 0.1;
        cfg.datasavedir         = datasavedir{ipatient};
        cfg.imagesavedir        = imagesavedir{ipatient};
        plotdata_corr(cfg,dat_micro,dat_macro);
        
        %% write  data for spyking-circus analysis for single patterns
        cfg                     = [];
        cfg.force               = true;
        cfg.forcereload         = true;
        cfg.prefix              = [num2str(ipatient),'-'];
        cfg.hpfilter            = 'yes';
        cfg.hpfreq              = 100;
        cfg.label               = label{ipatient}{imarker};
        cfg.startend            = startend{ipatient}(imarker,:);
        cfg.prestim             = prestim{ipatient}(imarker);
        cfg.poststim            = poststim{ipatient}(imarker);
        cfg.channel             = micro_labels{ipatient,imarker};
        cfg.datasavedir         = datasavedir{ipatient};
        [trialinfo_single{imarker},circuspath{imarker}]  = writemicroforspykingcircus_singlepattern(cfg, MuseStruct_micro);
        
        %         % read and plot spike data from Spyking Circus
        %         cfg                     = [];
        %         cfg.force               = false; % not yet implemented
        %         cfg.startend            = startend{ipatient}(imarker,:);
        %         cfg.prestim             = prestim{ipatient}(imarker);
        %         cfg.poststim            = poststim{ipatient}(imarker);
        %         cfg.imagesavedir        = imagesavedir{ipatient};
        %         cfg.datasavedir         = datasavedir{ipatient};
        %         cfg.label               = label{ipatient}{imarker};
        %         cfg.channel             = channel{ipatient}{1};
        %         cfg.resamplefs          = resamplefs{ipatient};
        %         cfg.trialinfo           = trialinfo_single{imarker};
        %         cfg.suffix              = '-1';
        %         [SpikeRaw, SpikeTrials] = readSpykingCircus_selected(cfg,MuseStruct_micro_aligned);
        
        close all
    end % imarker
    
    %% write concatinated data for spyking-circus -  concatinating all markers
    cfg                     = [];
    cfg.force               = true;
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.circuspath          = circuspath;
    cfg.channel             = micro_labels{ipatient,imarker};
    cfg.prefix              = [num2str(ipatient),'-'];
    trialinfo_allmarkers    = writemicroforspykingcircus_allmarkers(cfg);
    
    %%
    fname_trialinfos = fullfile(datasavedir{ipatient},[cfg.prefix,'trialinfos-',num2str(ipatient),'.mat']);
%     load(fname_trialinfos,'trialinfo_single','trialinfo_allmarkers');
    save(fname_trialinfos,'trialinfo_single','trialinfo_allmarkers');
    
    %% read and plot spike data from Spyking Circus - reads for all concatinated markers
    cfg                     = [];
    cfg.force               = true;
    cfg.startend            = startend{ipatient};
    cfg.prestim             = prestim{ipatient}; % +pad{ipatient};
    cfg.poststim            = poststim{ipatient}; % +pad{ipatient};
    cfg.imagesavedir        = imagesavedir{ipatient};
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.label               = label{ipatient};
    cfg.channel             = channel{ipatient}{imarker};
    cfg.resamplefs          = resamplefs{ipatient};
    cfg.trialinfo_all       = trialinfo_allmarkers;
    cfg.trialinfo_single    = trialinfo_single;
    cfg.timwin              = timwin(ipatient,:);
    cfg.suffix              = '-1';
    cfg.prefix              = [num2str(ipatient),'-'];
    [SpikeRaw, SpikeTrials] = readSpykingCircus_allmarkers(cfg,MuseStruct_micro);
    
    %% analyse patterns further separately
    
    for imarker = 1 : length(label{ipatient})
        
        %         %% read LFP data
        %         cfg                     = [];
        %         cfg.force               = false; % whether to force recalculations - does not make sense as one markers are not saved separately
        %         cfg.startend            = startend{ipatient}(imarker,:);
        %         cfg.label               = label{ipatient}{imarker};
        %         cfg.prestim             = prestim{ipatient}(imarker); %+pad{ipatient};
        %         cfg.poststim            = poststim{ipatient}(imarker); %+pad{ipatient};
        %         cfg.hpfilter            = hpfilter{ipatient}{imarker};
        %         cfg.hpfreq              = hpfreq{ipatient}{imarker};
        %         cfg.baseline            = 'no';
        %         cfg.micro_labels        = micro_labels{ipatient,imarker};
        %         cfg.macro_labels        = macro_labels{ipatient,imarker};
        %         cfg.resamplefs          = resamplefs{ipatient};
        %         cfg.datasavedir         = datasavedir{ipatient};
        %         [dat_micro, dat_macro]  = readEEG(cfg,MuseStruct_micro,MuseStruct_macro); % rename to readLFP
        %
        %         %% cross correlation
        %         cfg = [];
        %         cfg.binsize = 0.001;
        %         cfg.maxlag = 0.010;
        %
        %         cfg.debias = 'yes';
        %         cfg.method = 'xcorr';
        %         stat_x = ft_spike_xcorr(cfg,SpikeTrials{imarker});
        %         %
        %         %         cfg.method = 'shiftpredictor';
        %         %         stat_s      = ft_spike_xcorr(cfg,SpikeTrials{imarker});
        %         %         stat_s.xcorr = stat_s.shiftpredictor;
        %         %
        %         %         stat_diff = stat_x;
        %         %         stat_diff.xcorr = stat_x.xcorr - stat_s.shiftpredictor;
        %
        %         % plot crosscorrelation
        %         cfg = [];
        %         cfg.imagesavedir        = imagesavedir{ipatient};
        %         cfg.prefix = label{ipatient}{imarker};
        %         plotxcorr(cfg,stat_x);
        %
        %         figure;
        %         bar(stat_x.time,squeeze(stat_x.xcorr(4,7,:)),1,'facecolor',[0 0 0],'edgecolor',c);
        
        %% read data at original high samplerate
        cfgtemp                     = [];
        cfgtemp.force               = true;
        cfgtemp.startend            = startend{ipatient}(imarker,:);
        cfgtemp.prestim             = prestim{ipatient}(imarker); % +pad{ipatient};
        cfgtemp.poststim            = poststim{ipatient}(imarker); % +pad{ipatient};
        cfgtemp.datasavedir         = datasavedir{ipatient};
        cfgtemp.label               = label{ipatient}{imarker};
        cfgtemp.channel             = micro_labels{ipatient,imarker};
        cfgtemp.hpfilter            = 'no';
        cfgtemp.hpfreq              = 100;
        dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);
        
        
        %         % correlate between channels over trials
        %         for itrial = 1 : size(dat_microFs.trial,2)
        %             disp(num2str(itrial));
        %             c(itrial,:,:) = corr(dat_microFs.trial{itrial}');
        %         end
        %
        %         c_avg = squeeze(nanmean(c,1));
        %         c_std = squeeze(nanstd(c,1));
        %         c_3std = c_std * 3;
        %
        %         fig = figure;
        %         im = image(c_avg*255);
        %         colormap(jet(255));
        %         im.CDataMapping = 'scaled';
        %         colormap jet
        %         [H,P,CI,STATS] = ttest(c-0.5);
        %
        %         % print to file
        %         set(fig,'PaperOrientation','landscape');
        %         set(fig,'PaperUnits','normalized');
        %         set(fig,'PaperPosition', [0 0 1 1]);
        %         print(fig, '-dpdf', fullfile(imagesavedir{ipatient},'correlation.pdf'),'-r600');
        %         set(fig,'PaperOrientation','portrait');
        %         print(fig, '-dpng', fullfile(imagesavedir{ipatient},'correlation.png'),'-r600');
        %
        %         %% connectivity analysis
        %         cfg = [];
        %         cfg.method = 'mtmfft';
        %         cfg.taper = 'hanning';
        %         cfg.output = 'fourier';
        %         cfg.foi = 1:100;
        %         cfg.pad = 'nextpow2';
        %         freq = ft_freqanalysis(cfg,dat_microFs);
        %
        %         cfg = [];
        %         cfg.method = 'coh';
        %         coh = ft_connectivityanalysis(cfg,freq);
        %
        %         fig = figure;
        %         cfg = [];
        %         cfg.parameter = 'cohspctrm';
        %         cfg.zlim = [0 1];
        %         cfg.xlim = [1 40];
        %         ft_connectivityplot(cfg, coh);
        %
        %
        %         % print to file
        %         set(fig,'PaperOrientation','landscape');
        %         set(fig,'PaperUnits','normalized');
        %         set(fig,'PaperPosition', [0 0 1 1]);
        %         print(fig, '-dpdf', fullfile(imagesavedir{ipatient},'coherence.pdf'),'-r600');
        %         print(fig, '-dpng', fullfile(imagesavedir{ipatient},'coherence.png'),'-r600');
        
        %% combining LFP with spikes
        
        fname_trialinfos = fullfile(datasavedir{ipatient},[num2str(ipatient),'-trialinfos-',num2str(ipatient),'.mat']);
        load(fname_trialinfos,'trialinfo_single','trialinfo_allmarkers');
        prefix                      = [num2str(ipatient),'-'];
        
        % add header - better to do it in readMicroFs
        temp                        = dir(fullfile(datasavedir{ipatient},[prefix,'all_concatinated_',channel{ipatient}{imarker}(1:end-2),'_*.ncs']));
        hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
        dat_microFs.hdr             = hdr;
        
        % recreate trl
        trl                         = trialinfo_single{imarker} + trialinfo_allmarkers(imarker,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
        trl(:,3)                    = -ones(size(trl,1),1) * prestim{ipatient}(imarker) * hdr.Fs;
        
        %         trl = [];
        %         trl(:,1) = cumsum((dat_microFs.sampleinfo(:,2) - dat_microFs.sampleinfo(:,1)));
        %         trl(:,2) = circshift(trl(:,1),-1)-1;
        %         trl = circshift(trl,1);
        %         trl(1,1) = 1;
        %
        dat_microFs.cfg.trl         = trl;
        
        % combine spike and LFP
        data_all = ft_appendspike([],dat_microFs,SpikeTrials{imarker});
        %
        %         % plot spikes together with lfp
        %         cfg                 = [];
        %         cfg.viewmode        = 'vertical';
        %         cfg.continuous      = 'yes';
        %         cfg.mychan          = ft_channelselection('temp*', data_all, 'all');
        %         cfg.mychanscale     = repmat(max(max(data_all.trial{1}))*2,size(cfg.mychan,1));
        %         cfg.channelcolormap = [0,0,0;1,0,0];
        %         cfg.colorgroups     = [1,1,1,1,1,1,1,1,2,2,2,2,2,2];
        %         cfg.colorgroups     = [ones(1,size(data_all.label,1)-size(cfg.mychan,1)) ones(1,size(cfg.mychan,1))*2];
        %
        %         cfg.linewidth       = 1;
        %         ft_databrowser(cfg,data_all);
        %
        %         %% plot trials with LFP and spikes
        %
        %         for itrial = 1:3
        %             h = 150;
        %             n = 1;
        %
        %             figure; hold;
        %             for ichannel = 1 : size(data_all.label,1)
        %                 if isempty(strfind(data_all.label{ichannel},'temp'))
        %                     i1 = find(data_all.time{itrial} >= -prestim{ipatient}(imarker),1,'first');
        %                     i2 = find(data_all.time{itrial} >=  prestim{ipatient}(imarker),1,'first');
        %                     plot(data_all.time{itrial}(i1:i2),data_all.trial{itrial}(ichannel,i1:i2) + n*h,'color','k');
        %                     n = n + 1;
        %                 end
        %             end
        %             ax = axis;
        %             ic = 1;
        %             cmap = jet(size(cfg.mychan,1));
        %             for ichannel = 1 : size(data_all.label,1)
        %                 if ~isempty(strfind(data_all.label{ichannel},'temp'))
        %                     i1 = find(data_all.time{itrial} >= -prestim{ipatient}(imarker),1,'first');
        %                     i2 = find(data_all.time{itrial} >=  prestim{ipatient}(imarker),1,'first');
        %                     plot(data_all.time{itrial}(i1:i2),data_all.trial{itrial}(ichannel,i1:i2) * n*h,'color',cmap(ic,:));
        %                     ic = ic + 1;
        %                 end
        %             end
        %             legend
        %         end
        %         ylabel('Trials');
        %         xlabel('Time (s)');
        %         axis tight
        
        
        template_labels = ft_channelselection('temp*', data_all, 'all');
        
        for itemp = 1 : size(template_labels,1)
            
            %% interpolate LFP around spikes
            cfg                         = [];
            cfg.method                  = 'nan';
            cfg.timwin                  = [-0.002 0.002];
            cfg.spikechannel            = template_labels{itemp};
            cfg.channel                 = cellstr(micro_labels{ipatient,imarker});
            data_i                      = ft_spiketriggeredinterpolation(cfg,data_all);
            
            
            cfg                         = [];
            cfg.timwin                  = [-1 1];
            cfg.spikechannel            = template_labels{itemp};
            cfg.channel                 = cellstr(micro_labels{ipatient,imarker});
            cfg.latency                 = [-prestim{ipatient}(imarker) poststim{ipatient}(imarker)];
            spikeavg                    = ft_spiketriggeredaverage(cfg, data_i);
            
            fig = figure;
            plot(spikeavg.time,spikeavg.avg);
            legend
            title(['Spike-triggered average ',label{ipatient}{imarker}]);
            
            % print to file
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(imagesavedir{ipatient},[num2str(ipatient),'-spiketriggered_average_',label{ipatient}{imarker},'_template',num2str(itemp),'.pdf']),'-r600');
            set(fig,'PaperOrientation','portrait');
            print(fig, '-dpng', fullfile(imagesavedir{ipatient},[num2str(ipatient),'-spiketriggered_average_',label{ipatient}{imarker},'_template',num2str(itemp),'.png']),'-r600');
        end % itemplate
        
    end % imarker
        
    
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.foi = 20:10:100;
    cfg.t_ftimwin = 5./cfg.foi;
    cfg.taper = 'hanning';
    stsConcol = ft_spiketriggeredspectrum(cfg,data_all);
    
    
    
    
end
close all

end % ipatient




%% test correlations between nodules in patient #2

cfgtemp                     = [];
cfgtemp.force               = true;
cfgtemp.datasavedir         = datasavedir{ipatient};
cfgtemp.label               = label{ipatient}{imarker};
cfgtemp.startend            = startend{ipatient};
cfgtemp.channel             = ["_TNmi_1","_TNmi_2","_TNmi_3","_TNmi_4","_TNmi_5","_TNmi_6","_TNmi_7","_TNmi_8", "mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);


















% ilabel = 1;
%
% cfgtemp                         = [];
% cfgtemp.trl                     = cfg.trialinfo_single{ilabel} + cfg.trialinfo_all(ilabel,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
% cfgtemp.trl(:,3)                = -ones(size(cfgtemp.trl,1),1) * cfg.prestim(ilabel) * hdr.Fs;
% cfgtemp.trlunit                 = 'samples';
% cfgtemp.hdr                     = hdr;
% SpikeTrials                     = ft_spike_maketrials(cfgtemp,SpikeRaw);


cfgtemp                     = [];
cfgtemp.force               = false;
cfgtemp.datasavedir         = datasavedir{ipatient};
cfgtemp.label               = label{ipatient}{ilabel};
dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);

prefix                      = [num2str(ipatient),'-'];

% add header - better to do it in readMicroFs
temp                        = dir(fullfile(datasavedir{ipatient},[prefix,'all_concatinated_',channel{ipatient}{imarker}(1:end-2),'_*.ncs']));

hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data

dat_microFs.hdr = hdr;
data_all = ft_appendspike([],dat_microFs,SpikeTrials);


cfgtemp                             = [];
cfgtemp.timwin                      = [-0.25 0.25];
cfgtemp.spikechannel                = channel{ipatient}{ilabel}; % first unit
cfgtemp.channel                     = micro_labels{ipatient,ilabel}; % first four chans
cfgtemp.latency                     = [-prestim{ipatient}(ilabel) poststim{ipatient}(ilabel)];
staPost                             = ft_spiketriggeredaverage(cfg, SpikeTrials);

dat_micro
%% write deadfile for Spyking Circus
% writeSpykingCircusDeadFile(MuseStruct_micro_aligned)
% writeSpykingCircusDeadFile_concatinated(MuseStruct_micro_aligned)


