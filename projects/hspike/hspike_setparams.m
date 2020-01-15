
%% Setting parameters

function [config] = hspike_setparams(config)

disp('setting parameters');

muscale = 50; % multiunit scale

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
    os                  = 'windows';
else
    error('Platform not supported')
end

% Patient 1
config{1}.os                        = os;
config{1}.name                      = {'Hspike','SpikeDetect'};
config{1}.prefix                    = 'P1-'; % edit in code
config{1}.muse.startend             = {'Hspike','Hspike'; 'SpikeDetect','SpikeDetect'};   % start and end Muse marker
config{1}.muse.backupdir             = fullfile(rootpath_analysis, 'markerbackup');

config{1}.patientdir                = fullfile(rootpath_data,       'pat_02711_1193');
config{1}.rawdir                    = fullfile(rootpath_data,       'pat_02711_1193', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis,   'data',   'hspike');         % where to write data
config{1}.imagesavedir              = fullfile(rootpath_analysis,   'images', 'hspike');       % where to print images

% config{1}.directory_searchstring    = '02711_2019-04-*';
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
                                                                   
config{1}.labels.micro              = {'mHaT2_1','mHaT2_3','mHaT2_4','mHaT2_6','mHaT2_8'};
config{1}.labels.macro              = {'_HaT2_1','_HaT2_2'};

config{1}.hyp.imagesavedir          = fullfile(rootpath_analysis, 'images', 'hspike'); 
config{1}.hyp.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{1}.hyp.markerdir             = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.hyp.micromedchannel       = 'F3p6';
config{1}.hyp.contains              = {'NO_SCORE','AWAKE','PHASE_1','PHASE_2','PHASE_3','REM'}; % in order of height in plot
config{1}.hyp.markers               = {'Hspike','SpikeDetect'};
config{1}.hyp.overwrite             = 'append'; % 'append' or 'overwrite'
config{1}.hyp.spikewindow           = 60; % seconds
config{1}.hyp.spikewindowoverlap    = 0.5;
% config{1}.hyp.notcontains           = {"ADStartLoss","ADEndLoss","TTL","StartRecord","StopRecord","NLXEvent","BAD"};

config{1}.align.name                = {'Hspike','SpikeDetect'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.align.channel             = {'_HaT2_1','_HaT2_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.align.abs                 = {'no','no'};
config{1}.align.method              = {'crawlback','min'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{1}.align.filter              = {'bp','bp'};
config{1}.align.freq                = {[1, 10],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
config{1}.align.hilbert             = {'no','no'};
config{1}.align.thresh              = [1, 0];
config{1}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toiactive{1}        = [-0.05,  0.150];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toiactive{2}        = [-0.05,  0.05];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 1000;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow{1}     = [-2, -1];
config{1}.LFP.baselinewindow{2}     = [-2, -1];
config{1}.LFP.slidestep             = [0.01];

config{1}.epoch.toi{1}              = [-0.5  1];                                                                
config{1}.epoch.toi{2}              = [-0.5  1];                                                                
config{1}.epoch.pad                 = [0.5, 0.5, 0.5];

config{1}.circus.channel            = {'mHaT2_1','mHaT2_3','mHaT2_4','mHaT2_6','mHaT2_8'};
config{1}.circus.reref              = 'yes';
config{1}.circus.refchan            = 'mHaT2_2';
config{1}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'hspike', 'SpykingCircus');
config{1}.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
config{1}.circus.hpfreq             = 0; % even when not using
config{1}.circus.postfix            = '-1'; % after using circus-gui-matlab's SAVE number 

config{1}.TFR.channel               = {'_HaT2_1','_HaT2_6'};

config{1}.spike.slidestep           = [0.001];
config{1}.spike.toispikerate{1}     = [-2, 1];         % for plotting spikerate
config{1}.spike.resamplefs          = 1000;
config{1}.spike.width               = 15;
config{1}.spike.ISIbins             = [0 : 0.05 : 0.150];

config{1}.stats.bltoi{1}            = [-0.5, -0.1];  % has to fall within cfg.epoch.toi
config{1}.stats.actoi{1}            = [0,    1];  % has to fall within cfg.epoch.toi 
config{1}.stats.bltoi{2}            = [-0.5, -0.1];  % has to fall within cfg.epoch.toi
config{1}.stats.actoi{2}            = [0,    1];  % has to fall within cfg.epoch.toi 
config{1}.stats.alpha               = 0.025;


%% Patient 2
config{2}.os                        = os;
config{2}.name                      = {'Hspike','SpikeDetect'};
config{2}.prefix                    = '2718-'; % edit in code
config{2}.muse.startend             = {'SpikeHaT1_1','SpikeHaT1_1'; 'SpikeDetect','SpikeDetect'};   % start and end Muse marker
config{2}.muse.backupdir             = fullfile(rootpath_analysis, 'markerbackup');

config{2}.patientdir                = fullfile(rootpath_data,       'pat_02718_1201');
config{2}.rawdir                    = fullfile(rootpath_data,       'pat_02718_1201', 'eeg');
config{2}.datasavedir               = fullfile(rootpath_analysis,   'data',   'hspike');         % where to write data
config{2}.imagesavedir              = fullfile(rootpath_analysis,   'images', 'hspike');       % where to print images

config{2}.directorylist{1}          = {'02718_2019-05-14_12-31',...
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
config{2}.directorylist{2}          = {'02718_2019-05-15_13-52',...
                                       '02718_2019-05-15_15-52',...
                                       '02718_2019-05-15_17-52',...
                                       '02718_2019-05-15_19-52',...
                                       '02718_2019-05-15_21-52',...
                                       '02718_2019-05-15_23-52',...
                                       '02718_2019-05-16_01-52',...
                                       '02718_2019-05-16_03-52',...
                                       '02718_2019-05-16_05-52',...
                                       '02718_2019-05-16_07-52',...
                                       '02718_2019-05-16_09-52',...
                                       '02718_2019-05-16_10-16',...
                                       '02718_2019-05-16_11-15'};
config{2}.directorylist{3}          = {'02718_2019-05-16_13-15',...
                                       '02718_2019-05-16_15-15',...
                                       '02718_2019-05-16_17-15',...
                                       '02718_2019-05-16_19-15',...
                                       '02718_2019-05-16_21-15',...
                                       '02718_2019-05-16_23-15',...
                                       '02718_2019-05-17_01-15',...
                                       '02718_2019-05-17_03-15',...
                                       '02718_2019-05-17_05-15',...
                                       '02718_2019-05-17_07-15',...
                                       '02718_2019-05-17_09-15',...
                                       '02718_2019-05-17_10-12'};          

                                    
config{2}.labels.micro              = {'mHaT1_8'};
config{2}.labels.macro              = {'_HaT1_1','_HaT1_2','_HaT1_3'};

config{2}.align.name                = {'Hspike','SpikeDetect'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.align.channel             = {'_HaT1_1','_HaT1_1'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.align.abs                 = {'no','no'};
config{2}.align.method              = {'crawlback','min'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
config{2}.align.filter              = {'bp','bp'};
config{2}.align.freq                = {[1, 10],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
config{2}.align.hilbert             = {'no','no'};
config{2}.align.thresh              = [1, 0];
config{2}.align.toiplot{1}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toiactive{1}        = [-0.05,  0.150];                                         % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toibaseline{1}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toiplot{2}          = [-0.3,  0.7];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toiactive{2}        = [-0.05,  0.05];                                          % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toibaseline{2}      = [-0.3, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

config{2}.LFP.hpfilter              = 'no';
config{2}.LFP.hpfreq                = 1;
config{2}.LFP.resamplefs            = 1000;
config{2}.LFP.baseline              = 'yes';
config{2}.LFP.baselinewindow{1}     = [-2, -1];
config{2}.LFP.baselinewindow{2}     = [-2, -1];
config{2}.LFP.slidestep             = [0.01];

config{2}.epoch.toi{1}              = [-0.5  1];                                                                
config{2}.epoch.toi{2}              = [-0.5  1];                                                                
config{2}.epoch.pad                 = [0.5, 0.5, 0.5];

config{2}.spikedetect.LS            = 5;    % Left half-wave slope; default: 7
config{2}.spikedetect.RS            = 5;    % Right half-wave slope; default: 7
config{2}.spikedetect.TAMP          = 400;  % Total amplitude; default: 600
config{2}.spikedetect.LD            = 1;    % Left half-wave duration; default = 10
config{2}.spikedetect.RD            = 1;    % Right half-wave duration; default = 10
config{2}.spikedetect.STDCoeff      = 4;    % Chebyshev inequality coefficient (distance from centre point or mean); default 4
config{2}.spikedetect.SCALE         = 70;   % Scaling parameter
config{2}.spikedetect.BlockSize     = 1;    % Data processing block size in minutes
config{2}.spikedetect.TroughSearch  = 40;   % distance in ms to search for a trough on each side of a detected peak
config{2}.spikedetect.DetThresholds = [config{2}.spikedetect.LS; config{2}.spikedetect.RS; config{2}.spikedetect.TAMP; config{2}.spikedetect.LD; config{2}.spikedetect.RD;];
config{2}.spikedetect.FilterSpec    =  [20; 50; 1; 35;];


config{2}.hyp.imagesavedir          = fullfile(rootpath_analysis, 'images', 'hspike'); 
config{2}.hyp.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{2}.hyp.markerdir             = fullfile(rootpath_analysis, 'data',   'hspike');
config{2}.hyp.micromedchannel       = 'HaT1';
config{2}.hyp.contains              = {'NO_SCORE','AWAKE','PHASE_1','PHASE_2','PHASE_3','REM'}; % in order of height in plot
config{2}.hyp.markers               = {'SpikeDetect'};
config{2}.hyp.overwrite             = 'append'; % 'append' or 'overwrite'
config{2}.hyp.spikewindow           = 60; % seconds
config{2}.hyp.spikewindowoverlap    = 0.5;


end

