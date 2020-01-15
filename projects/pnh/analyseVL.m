%% Analysis script for Virgini Lambreque's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario Chavez
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/

addpath /Volumes/iss01.charpier/stephen.whitmarsh/data/02230_2015-02-25_14-36/m1pNs
addpath /Volumes/iss01.charpier/stephen.whitmarsh/fieldtrip/

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
maxNumCompThreads(8)

% Setting parameters

% Patient 1, perivetricular heterotopia #1
hpfilter{1}                = {'no','no','no','no'};
hpfreq{1}                  = {1,1,1,1};
label{1}                   = { 'inRD','RDsS','P','PP'};                                                                      % marker name
startend{1}                = { 'inRD__START__','inRD__END__';'RDsS__START__','RDsS__END__';'P','P';'PP__START__','PP__END__'};   % start and end Muse marker
prestim{1}                 = [0, 0, 0, 0];                                                                  % list of onset timing with respect to start-marker (s)
poststim{1}                = [0, 0, 0, 0];                                                                   % list of offset timing with respect to end-marker (s)
pad{1}                     = 0.25;
slidestep{1}               = [0.01,0.001,0.001,0.001];
lpfreq{1}                  = [4,4,40,40];                                                                                  % lowpass filter freq to smooth peak detection (Hz)
toiactive{1}               = [-0.1, 0.8;  -0.1,  0.3;  -0.1, 0.1;   0.0,  0.3];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
toibaseline{1}             = [-1.0,-0.1;  -1.0, -0.1;  -1.0 -0.1;  -1.0,  0.3];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
resamplefs{1}              = 640;                                                                                          % resample data to speed up and reduce memory (Hz), and align micro and macro
channel{1}                 = {'m1pNs_4','m1pNs_4','m1pNs_4','m1pNs_4'};                                                                                    % pattern to identify channel on which to based peak detection
alignmethod{1}             = {'first','first','max','first'};                                                              % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
pkthresh{1}                = [0,0.25,0.25,0.25];                                                                           % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
datasavedir{1}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02316_0761/analysis/data';         % where to write data
imagesavedir{1}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02316_0761/analysis/images';       % where to print images
patient_directory{1}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02316_0761/eeg';
directory_searchstring{1}  = '02316*';
micro_searchstring{1}{1}   = '*mHa2d*.ncs';
macro_searchstring{1}{1}   = '*_Ha2d*.ncs';
micro_labels{1,1}          = ["mHa2d_1","mHa2d_2","mHa2d_3","mHa2d_4","mHa2d_5","mHa2d_6","mHa2d_7","mHa2d_8"];
macro_labels{1,1}          = ["_Ha2d_1","_Ha2d_2","_Ha2d_3","_mHa2d_4","_Ha2d_5","_Ha2d_6","_Ha2d_7","_Ha2d_8"];
micro_labels{1,2}          = ["mHa2d_1","mHa2d_2","mHa2d_3","mHa2d_4","mHa2d_5","mHa2d_6","mHa2d_7","mHa2d_8"];
macro_labels{1,2}          = ["_Ha2d_1","_Ha2d_2","_Ha2d_3","_mHa2d_4","_Ha2d_5","_Ha2d_6","_Ha2d_7","_Ha2d_8"];
% micro_labels{1,3}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
% macro_labels{1,3}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
% micro_labels{1,4}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
% macro_labels{1,4}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
timwin{1,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,2}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,3}                = [-0.0025 0.0025];  % for plotting spikerate
timwin{1,4}                = [-0.01 0.01];      % for plotting spikerate


ipatient = 1;
imarker = 1;


%% loop through patient analyses

for ipatient = 1 : 1
    
    for imarker = 1 : 2 %length(label{ipatient})
        
%         read muse markers
        MuseStruct_micro = readMuseMarkers(patient_directory{ipatient},directory_searchstring{ipatient},micro_searchstring{ipatient}{1},datasavedir{ipatient},'micro',false);
        MuseStruct_macro = readMuseMarkers(patient_directory{ipatient},directory_searchstring{ipatient},macro_searchstring{ipatient}{1},datasavedir{ipatient},'macro',false);
        
        % read LFP data
        cfg                     = [];
        cfg.force               = false; % whether to force recalculations
        cfg.startend            = startend{ipatient}(imarker,:);
        cfg.label               = label{ipatient}{imarker};
        cfg.prestim             = prestim{ipatient}(imarker); %+pad{ipatient};
        cfg.poststim            = poststim{ipatient}(imarker); %+pad{ipatient};
        cfg.hpfilter            = hpfilter{ipatient}{imarker};
        cfg.hpfreq              = hpfreq{ipatient}{imarker};
        cfg.baseline            = 'no';
        cfg.micro_labels        = micro_labels{ipatient,imarker};
        cfg.macro_labels        = macro_labels{ipatient,imarker};
        cfg.resamplefs          = resamplefs{ipatient};
        cfg.datasavedir         = datasavedir{ipatient};
        [dat_micro, dat_macro]  = readEEG(cfg,MuseStruct_micro,MuseStruct_macro); % rename to readLFP
        
%         write concatinated data for spyking-circus
        cfg                     = [];
        cfg.force               = true;
        cfg.forcereload         = true;
        cfg.hpfilter            = 'yes';
        cfg.hpfreq              = 100;
        cfg.label               = label{ipatient}{imarker};
        cfg.startend            = startend{ipatient}(imarker,:);
        cfg.prestim             = prestim{ipatient}(imarker);
        cfg.poststim            = poststim{ipatient}(imarker);
        cfg.channel             = micro_labels{ipatient,imarker};
        cfg.datasavedir         = datasavedir{ipatient};
        [trialinfo_single{imarker},circuspath{imarker}]  = writemicroforspykingcircus(cfg, MuseStruct_micro);
        
    end
  
    % write concatinated data for spyking-circus -  concatinating all markers
    cfg                     = [];
    cfg.force               = true;
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.circuspath          = circuspath;
    cfg.channel             = micro_labels{ipatient,imarker};
    trialinfo_allmarkers    = writemicroforspykingcircus_allmarkers(cfg);
    
    fname_trialinfos = fullfile(datasavedir{ipatient},'trialinfos_VL.mat');
    save(fname_trialinfos,'trialinfo_single','trialinfo_allmarkers');
    
    fname_ncs = fullfile(datasavedir{ipatient},'all_concatinated_mHa2d_1.ncs');   
    dat = ft_read_data(fname_ncs);
    figure; plot(dat)
    
    %% read and plot spike data from Spyking Circus
    cfg                     = [];
    cfg.force               = false; % not yet implemented
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
    [SpikeRaw, SpikeTrials] = readSpykingCircus_allmarkers(cfg,MuseStruct_micro_aligned);
%     [SpikeRaw, SpikeTrials] = readSpykingCircus_allmarkers(cfg,MuseStruct_micro);
%     
    close all
    
end % ipatient







%% write deadfile for Spyking Circus
% writeSpykingCircusDeadFile(MuseStruct_micro_aligned)
% writeSpykingCircusDeadFile_concatinated(MuseStruct_micro_aligned)


