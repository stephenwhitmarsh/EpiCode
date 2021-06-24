%% Stimulation micro analysis script
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
% addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/scottclowe-matlab-schemer-f8115af/ %https://github.com/scottclowe/matlab-schemer/blob/master/schemes/README.md#matlab-schemes
% schemer_import('/network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/scottclowe-matlab-schemer-f8115af/schemes/darksteel.prf')

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/

% addpath /Volumes/iss01.charpier/stephen.whitmarsh/data/02230_2015-02-25_14-36/m1pNs
% addpath /Volumes/iss01.charpier/stephen.whitmarsh/fieldtrip/

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

patientdir{1}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/';
imagesavedir{1}     = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';
datasavedir{1}      = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';
prestim{1}          = 10;
poststim{1}         = 15;
micro_labels{1}     = {'mHimg_2','mHimg_3','mHimg_4','mHimg_5','mHimg_6','mHimg_7','mHimg_8', ...
    'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_6','mTNmi_7','mTNmi_1',}; % mCasd too noisy
% micro_labels{1}     = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_6','mTNmi_7'}; % mCasd too noisy
timecode{1}         = '*SYNC*';
timecodeori{1}      = -1;
slidestep{1}        = 0.05;


patientdir{2}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02619_1078/';
imagesavedir{2}     = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02619_1078/analysis/images';
datasavedir{2}      = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02619_1078/analysis/data';
prestim{2}          = 2;
poststim{2}         = 7;
micro_labels{2}     = {'mAmT2_1','mAmT2_2','mAmT2_3','mAmT2_4','mAmT2_5','mAmT2_6','mAmT2_7', 'mAmT2_8', ...
    'mHaTB_1','mHaTB_2','mHaTB_3','mHaTB_4','mHaTB_5','mHaTB_6','mHaTB_7', 'mHaTB_8'};

timecode{2}         = '*SLI*';
timecodeoriMM{2}    = -1;
timecodeoriNL{2}    = 1;


slidestep{2}        = 0.05;

ipatient = 1;

maxNumCompThreads(4)

%% synchronize between MicroMed notes and NeuraLynx data.
cfg                 = [];
cfg.force           = false;
cfg.patientdir      = patientdir{ipatient};
cfg.imagesavedir    = imagesavedir{ipatient};
cfg.datasavedir     = datasavedir{ipatient};
cfg.timecode        = timecode{ipatient};
cfg.timecodeoriMM   = timecodeoriMM{ipatient};
cfg.timecodeoriNL   = timecodeoriNL{ipatient};
[micromed_markers]  = synchronize_micromed_neuralynx(cfg);

%% read micro data segment for all stimulations
cfg                 = [];
cfg.force           = false;
cfg.patientdir      = patientdir{ipatient};
cfg.imagesavedir    = imagesavedir{ipatient};
cfg.datasavedir     = datasavedir{ipatient};
cfg.prestim         = prestim{ipatient}; % in seconds
cfg.poststim        = poststim{ipatient}; % in seconds
cfg.samplecolumnNL  = 'Sample_NL';
cfg.samplecolumnMM  = 'Sample';
cfg.postfix         = '';
[markerdat_micro,~] = readstimdata(cfg,micromed_markers);

%% align according to artefacts
cfg = [];
cfg.force           = true; % if false, only supply cfg argument
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 100;
cfg.channel_micro   = {'mHimg_3','mHimg_4','mHimg_5','mHimg_6','mHimg_7','mHimg_8'};
cfg.plotfigures     = false;
cfg.imagesavedir    = imagesavedir{ipatient};
cfg.datasavedir     = datasavedir{ipatient};
[micromed_markers_aligned,markerdat_micro_aligned] = alignstimdata(cfg,micromed_markers,markerdat_micro);
% [micromed_markers_aligned] = alignstimdata(cfg);

% read micro data segment for all stimulations - now with alignement 
cfg                 = [];
cfg.force           = true;
cfg.patientdir      = patientdir{ipatient};
cfg.imagesavedir    = imagesavedir{ipatient};
cfg.datasavedir     = datasavedir{ipatient};
cfg.prestim         = prestim{ipatient}; % in seconds
cfg.poststim        = poststim{ipatient}; % in seconds
cfg.samplecolumnNL  = 'Sample_NL_aligned';
cfg.samplecolumnMM  = 'Sample';
cfg.postfix         = 'aligned';
[markerdat_micro_aligned,~] = readstimdata(cfg,micromed_markers_aligned);

%% plot EEG data
cfg                     = [];
cfg.force               = true;
cfg.prestim             = prestim{ipatient};
cfg.poststim            = poststim{ipatient};
cfg.slidestep           = slidestep{ipatient};
cfg.binsize             = 0.1;
cfg.datasavedir         = datasavedir{ipatient};
cfg.imagesavedir        = imagesavedir{ipatient};
cfg.latency             = [-prestim{ipatient} poststim{ipatient}];
cfg.frequency           = 50;
cfg.current             = 'all';
cfg.channel_micro       = {'mHimg_3','mHimg_4','mHimg_5','mHimg_6','mHimg_7','mHimg_8'};
cfg.channel_micro       = {'mHimg_3'};
cfg.channel_macro       = {'Himg3','Himg4','Himg5','Himg6','Himg7','Himg8'};
cfg.contact1            = 'Himg_3';
cfg.contact2            = 'Himg_4';
cfg.channel_micro       = {'mCasd_1'};
cfg.channel_micro       = {'mTNmi_4'};
cfg.contact1            = 'Caid_6';
cfg.contact2            = 'Caid_7';
cfg.current             = 0.5;
cfg.markers             = micromed_markers_aligned;
cfg.hpfilter            = 'no';
cfg.hpfreq              = 100;
plotdata_stim_trial(cfg,markerdat_micro_aligned, []);

%% write data for Spyking-Circus
cfg                     = [];
cfg.force               = true;
cfg.patientdir          = patientdir{ipatient};
cfg.datasavedir         = datasavedir{ipatient};
cfg.channel             = micro_labels{ipatient};
cfg.start_deadtime      = prestim{ipatient} - 1;
cfg.end_deadtime        = prestim{ipatient} + 6;

cfg.poststim            = poststim{ipatient};
[trialinfo,dat]         = writemicrostimforspykingcircus(cfg, markerdat_micro_aligned, micromed_markers_aligned);

%% run Spyking-Circus
        
%% read and plot spike data from Spyking Circus
cfg                     = [];
cfg.force               = true; % not yet implemented
cfg.prestim             = prestim{ipatient};
cfg.poststim            = poststim{ipatient};
cfg.imagesavedir        = imagesavedir{ipatient};
cfg.datasavedir         = datasavedir{ipatient};
cfg.channel             = micro_labels{ipatient};
cfg.suffix              = '-1';
cfg.trialinfo           = trialinfo;
cfg.markers             = micromed_markers;
[SpikeRaw, SpikeTrials] = readSpykingCircus_selected_stim(cfg);