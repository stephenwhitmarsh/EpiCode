function [config] = hspike_setparams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [config] = hspike_setparams
%
% This function outputs all the settings of a study, to be defined below
%
% Note the consideration of the operating system, because the pointers to
% the file server is dealt with differently. This could be different for
% you.
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
else
    error('Platform not supported')
end

%% Patient 1

config{1}.name                      = {'Hspike'};
config{1}.prefix                    = '2711-';
config{1}.rawdir                    = fullfile(rootpath_data,     'pat_02711_1193', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data',   'hspike');        
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'hspike');     

config{1}.muse.startend             = {'Hspike','Hspike'}; 
config{1}.muse.backupdir            = fullfile(rootpath_analysis, 'markerbackup');

config{1}.hyp.imagesavedir          = fullfile(rootpath_analysis, 'images', 'hspike');
config{1}.hyp.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{1}.hyp.markerdir             = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.hyp.micromedchannel       = 'F3p6';
config{1}.hyp.contains              = {'NO_SCORE','AWAKE','PHASE_1','PHASE_2','PHASE_3','REM'}; 
config{1}.hyp.markers               = {'Hspike','SpikeDetect'};
config{1}.hyp.overwrite             = 'append';
config{1}.hyp.spikewindow           = 60;
config{1}.hyp.spikewindowoverlap    = 0.5;

config{1}.epoch.toi{1}              = [-0.5  1];
config{1}.epoch.pad{1}              = 0.5;

config{1}.LFP.name                  = {'Hspike'};
config{1}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 250;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow{1}     = [-0.5, -0.2];
config{1}.LFP.baselinewindow{2}     = [-0.5, -0.2];

config{1}.align.name                = {'Hspike'};
config{1}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.align.reref               = 'yes';
config{1}.align.refmethod           = 'bipolar';
config{1}.align.latency             = [-0.05 0.05];

config{1}.cluster.name              = {'HSpike'};
config{1}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.cluster.reref             = 'no';
config{1}.cluster.refmethod         = 'bipolar';
config{1}.cluster.resamplefs        = [];
config{1}.cluster.baselinewindow    = [-0.5 -0.2];
config{1}.cluster.latency           = [-0.5 1];
config{1}.cluster.dbscan            = 'no';
config{1}.cluster.kmeans            = 'no';
config{1}.cluster.kmedoids          = 'yes';
config{1}.cluster.N                 = [2, 6];

%% Patient 2

config{2}                           = config{1};
config{2}.name                      = {'SpikeHaT1_1'};
config{2}.prefix                    = '2718-'; % edit in code
config{2}.rawdir                    = fullfile(rootpath_data, 'pat_02718_1201', 'eeg');

config{2}.muse.startend             = {'SpikeHaT1_1','SpikeHaT1_1'};   % start and end Muse marker
config{2}.hyp.micromedchannel       = 'HaT11';
config{2}.hyp.markers               = {'SpikeHaT1_1'};

config{2}.LFP.name                  = {'SpikeHaT1_1'};
config{2}.LFP.channel               = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HaT1_5'};

config{2}.align.name                = {'SpikeHaT1_1'};
config{2}.align.channel             = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HaT1_5'};
config{2}.align.reref               = 'yes';
config{2}.align.refmethod           = 'bipolar';

config{2}.cluster.name              = {'SpikeHaT1_1'};
config{2}.cluster.channel           = {'_HaT1_1','_HaT1_2','_HaT1_3','_HaT1_4','_HaT1_5'};


%% Patient 3

config{3}                           = config{1};
config{3}.name                      = {'Hspike'};
config{3}.prefix                    = '2660-';
config{3}.rawdir                    = fullfile(rootpath_data, 'pat_02660_1136', 'eeg');

config{3}.muse.startend             = {'Hspike','Hspike'};   % start and end Muse marker

config{3}.hyp.micromedchannel       = 'Ha2g1';
config{3}.hyp.markers               = {'HSpike'};
                                    
config{3}.LFP.name                  = {'Hspike'};
config{3}.LFP.channel               = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Ha2g_5','_Ha2g_6','_Ha2g_7','_Ha2g_8','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4','_Hm2g_5','_Hm2g_6'};

config{3}.align.name                = {'Hspike'};  
config{3}.align.channel             = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Ha2g_5','_Ha2g_6','_Ha2g_7','_Ha2g_8','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4','_Hm2g_5','_Hm2g_6'};

config{3}.cluster.name              = {'HSpike'};
config{3}.cluster.channel           = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4','_Hm2g_1','_Hm2g_2','_Hm2g_3','_Hm2g_4'};

%% Patient 4

config{4}                           = config{1};
config{4}.name                      = {'Hspike'};
config{4}.prefix                    = '2680-'; % edit in code
config{4}.rawdir                    = fullfile(rootpath_data,       'pat_02680_1158', 'eeg');

config{4}.muse.startend             = {'Hspike','Hspike'};   % start and end Muse marker

config{4}.hyp.micromedchannel       = 'HaT21';
config{4}.hyp.markers               = {'HSpike'};
                                    
config{4}.LFP.name                  = {'Hspike'};
config{4}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT2_1','_HmT2_2','_HmT2_3'};

config{4}.align.name                = {'Hspike'};  
config{4}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT2_1','_HmT2_2','_HmT2_3'};

config{4}.cluster.name              = {'HSpike'};
config{4}.cluster.channel           = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HmT2_1','_HmT2_2','_HmT2_3'};

%% DATA   

config{1}.directorylist             = [];
config{1}.directorylist{1}          =  {'02711_2019-04-17_12-29',...
                                        '02711_2019-04-17_14-29',...
                                        '02711_2019-04-17_16-29',...
                                        '02711_2019-04-17_18-29',...
                                        '02711_2019-04-17_20-29',...
                                        '02711_2019-04-17_22-29',...
                                        '02711_2019-04-18_00-29',...
                                        '02711_2019-04-18_02-29',...
                                        '02711_2019-04-18_04-29',...
                                        '02711_2019-04-18_06-29',...
                                        '02711_2019-04-18_08-29',...
                                        '02711_2019-04-18_10-29',...
                                        '02711_2019-04-18_11-04'};
config{1}.directorylist{2}          =  {'02711_2019-04-18_13-04',...
                                        '02711_2019-04-18_15-04',...
                                        '02711_2019-04-18_17-04',...
                                        '02711_2019-04-18_19-04',...
                                        '02711_2019-04-18_21-04',...
                                        '02711_2019-04-18_23-04',...
                                        '02711_2019-04-19_01-04',...
                                        '02711_2019-04-19_03-04',...
                                        '02711_2019-04-19_05-04',...
                                        '02711_2019-04-19_07-04',...
                                        '02711_2019-04-19_09-04',...
                                        '02711_2019-04-19_10-00',...
                                        '02711_2019-04-19_12-00'};
config{1}.directorylist{3}          =  {'02711_2019-04-19_14-00',...
                                        '02711_2019-04-19_16-00',...
                                        '02711_2019-04-19_18-00',...
                                        '02711_2019-04-19_20-00',...
                                        '02711_2019-04-19_22-00',...
                                        '02711_2019-04-20_00-00',...
                                        '02711_2019-04-20_02-00',...
                                        '02711_2019-04-20_04-00',...
                                        '02711_2019-04-20_06-00',...
                                        '02711_2019-04-20_08-00',...
                                        '02711_2019-04-20_10-00',...
                                        '02711_2019-04-20_12-00'};

config{2}.directorylist             = [];
config{2}.directorylist{1}          =  {'02718_2019-05-14_12-31',...
                                        '02718_2019-05-14_14-31',...
                                        '02718_2019-05-14_16-31',...
                                        '02718_2019-05-14_18-31',...
                                        '02718_2019-05-14_20-31',...
                                        '02718_2019-05-14_22-31',...
                                        '02718_2019-05-15_00-31',...
                                        '02718_2019-05-15_02-31',...
                                        '02718_2019-05-15_04-31',...
                                        '02718_2019-05-15_06-31',...
                                        '02718_2019-05-15_08-31',...
                                        '02718_2019-05-15_10-31',...
                                        '02718_2019-05-15_10-44',...
                                        '02718_2019-05-15_11-52'};
config{2}.directorylist{2}          =  {'02718_2019-05-15_13-52',...
                                        '02718_2019-05-15_15-52',...
                                        '02718_2019-05-15_17-52',... % VERY ARTEFACTED
                                        '02718_2019-05-15_19-52',...
                                        '02718_2019-05-15_21-52',...
                                        '02718_2019-05-15_23-52',...
                                        '02718_2019-05-16_01-52',...
                                        '02718_2019-05-16_03-52',...
                                        '02718_2019-05-16_05-52',...
                                        '02718_2019-05-16_07-52',...
                                        '02718_2019-05-16_09-52',...
                                        '02718_2019-05-16_10-16',...
                                        '02718_2019-05-16_11-15'};   % very artefacted for most of the recording
config{2}.directorylist{3}          =  {'02718_2019-05-16_13-15',...
                                        '02718_2019-05-16_15-15',...
                                        '02718_2019-05-16_17-15',... % very artefacted
                                        '02718_2019-05-16_19-15',... % very artefacted
                                        '02718_2019-05-16_21-15',... % very artefacted for most of the recording
                                        '02718_2019-05-16_23-15',...
                                        '02718_2019-05-17_01-15',...
                                        '02718_2019-05-17_03-15',...
                                        '02718_2019-05-17_05-15',...
                                        '02718_2019-05-17_07-15',...
                                        '02718_2019-05-17_09-15',...
                                        '02718_2019-05-17_10-12'};                                                                  
  
config{3}.directorylist             = [];
config{3}.directorylist{1}          = { '02660_2018-11-13_16-08',...
                                        '02660_2018-11-13_18-08',...
                                        '02660_2018-11-13_20-08',...
                                        '02660_2018-11-13_22-08',...
                                        '02660_2018-11-14_00-08',...
                                        '02660_2018-11-14_02-08',...
                                        '02660_2018-11-14_04-08',...
                                        '02660_2018-11-14_06-08',...
                                        '02660_2018-11-14_08-08',...
                                        '02660_2018-11-14_09-29',...
                                        '02660_2018-11-14_10-56',...
                                        '02660_2018-11-14_11-51',...
                                        '02660_2018-11-14_13-51',...
                                        '02660_2018-11-14_15-51'};
                                    
config{4}.directorylist             = [];                             
config{4}.directorylist{1}          = { '02680_2019-01-15_12-45'...
                                        '02680_2019-01-15_14-45'...
                                        '02680_2019-01-15_15-31'...
                                        '02680_2019-01-15_17-31'...
                                        '02680_2019-01-15_19-31'...
                                        '02680_2019-01-15_21-31'...
                                        '02680_2019-01-15_23-31'...
                                        '02680_2019-01-16_01-31'...
                                        '02680_2019-01-16_03-31'...
                                        '02680_2019-01-16_05-31'...
                                        '02680_2019-01-16_07-31'...
                                        '02680_2019-01-16_09-31'...
                                        '02680_2019-01-16_09-52'...
                                        '02680_2019-01-16_10-58'...
                                        '02680_2019-01-16_11-32'};
                            
%%                                   
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
% config{1}.hyp.notcontains         = {"ADStartLoss","ADEndLoss","TTL","StartRecord","StopRecord","NLXEvent","BAD"};

% config{1}.align.name                = {'Hspike','SpikeDetect'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.channel             = {'_HaT2_1','_HaT2_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.abs                 = {'no','no'};
% config{1}.align.method              = {'crawlback','min'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{1}.align.filter              = {'bp','bp'};
% config{1}.align.freq                = {[1, 10],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{1}.align.hilbert             = {'no','no'};
% config{1}.align.thresh              = [1, 0];
% config{1}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{1}        = [-0.05,  0.150];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{2}        = [-0.05,  0.05];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% 
% config{1}.circus.channel            = {'mHaT2_1','mHaT2_3','mHaT2_4','mHaT2_6','mHaT2_8'};
% config{1}.circus.reref              = 'yes';
% config{1}.circus.refchan            = 'mHaT2_2';
% config{1}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'hspike', 'SpykingCircus');
% config{1}.circus.hpfilter           = 'no'; 
% config{1}.circus.hpfreq             = 0;    
% config{1}.circus.postfix            = '-1';
% 
% config{1}.TFR.channel               = {'_HaT2_1','_HaT2_6'};
% 
% config{1}.spike.slidestep           = 0.001;
% config{1}.spike.toispikerate{1}     = [-2, 1];         
% config{1}.spike.resamplefs          = 1000;
% config{1}.spike.width               = 15;
% config{1}.spike.ISIbins             = 0 : 0.05 : 0.150;
% 
% config{1}.stats.bltoi{1}            = [-0.5, -0.1];  
% config{1}.stats.actoi{1}            = [0,    1];  
% config{1}.stats.bltoi{2}            = [-0.5, -0.1]; 
% config{1}.stats.actoi{2}            = [0,    1]; 
% config{1}.stats.alpha               = 0.025;

% config{1}.spikedetect.LS            = 5;    % Left half-wave slope; default: 7
% config{1}.spikedetect.RS            = 5;    % Right half-wave slope; default: 7
% config{1}.spikedetect.TAMP          = 400;  % Total amplitude; default: 600
% config{1}.spikedetect.LD            = 1;    % Left half-wave duration; default = 10
% config{1}.spikedetect.RD            = 1;    % Right half-wave duration; default = 10
% config{1}.spikedetect.STDCoeff      = 4;    % Chebyshev inequality coefficient (distance from centre point or mean); default 4
% config{1}.spikedetect.SCALE         = 70;   % Scaling parameter
% config{1}.spikedetect.BlockSize     = 1;    % Data processing block size in minutes
% config{1}.spikedetect.TroughSearch  = 40;   % distance in ms to search for a trough on each side of a detected peak
% config{1}.spikedetect.FilterSpec    = [20; 50; 1; 35;];
% 
% config{2}.align.name                = {'SpikeHaT1_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.align.channel             = {'_HaT1_3','_HaT1_3'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.align.abs                 = {'no','no'};
% config{2}.align.method              = {'nearest','max'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{2}.align.filter              = {'bp','bp'};
% config{2}.align.freq                = {[1, 40],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{2}.align.hilbert             = {'no','no'};
% config{2}.align.thresh              = [0, 0];
% config{2}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiactive{1}        = [-0.05,  0.150];                                         % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{2}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiactive{2}        = [-0.05,  0.05];                                          % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{2}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% % 
% config{2}.circus.channel            = {'mHaT1_7'}; 
% config{2}.circus.reref              = 'no';
% config{2}.circus.refchan            = 'mHaT1_1';
% config{2}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'hspike', 'SpykingCircus');
% config{2}.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
% config{2}.circus.hpfreq             = 0; % even when not using
% config{2}.circus.postfix            = '-1'; % after using circus-gui-matlab's SAVE number

