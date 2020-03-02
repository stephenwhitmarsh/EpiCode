
%% Setting parameters DTX project Paul Baudin

function [config] = dtx_setparams_probe(config)

disp('setting parameters');
muscale = 50; % multiunit scale

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-PROBE\Analyses_Paul';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-raw\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data'); %removed of config{i} so we can more easily modify it 
imagesavedir = rootpath_analysis; 

%% Config common
configcommon.os                        = os;
configcommon.name                      = {'SlowWave','Seizure','InterIctal'};
configcommon.datasavedir               = datasavedir;
configcommon.muse.startend             = {'SlowWave','SlowWave'; 'Crise_Start','Crise_End'; 'Crise_End','SlowWave'};   % 'SlowWave','SlowWave'; for readLFP function : cut data ...s before SlowWave, and ...s after SlowWave
configcommon.align.name                = {'SlowWave'};
configcommon.align.flip                = {'no'};
configcommon.align.abs                 = {'no'};
configcommon.align.method              = {'nearest'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
configcommon.align.filter              = {'lp'};
configcommon.align.freq                = {5};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
configcommon.align.hilbert             = {'no'};
configcommon.align.thresh              = [0];
configcommon.align.toiplot{1}          = [-1,  1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
configcommon.align.toiactive{1}        = [-0.5, 0.5];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
configcommon.align.toibaseline{1}      = [-1, -0.5];         
% baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

configcommon.LFP.name                  = {'SlowWave'};
configcommon.LFP.hpfilter              = 'no';
configcommon.LFP.hpfreq                = 1;
configcommon.LFP.resamplefs            = 320; %because sampling rate is 3200Hz
configcommon.LFP.baseline              = 'yes';
configcommon.LFP.baselinewindow{1}     = [-2, -1];
configcommon.LFP.baselinewindow{2}     = [-2, -1];
configcommon.LFP.baselinewindow{3}     = [0, 1];
configcommon.LFP.slidestep             = 0.01;
configcommon.LFP.electrodeToPlot       = [6, 11]; %index of the channel in configcommon.labels.macro  

% list of onset timing with respect to start-marker (s)
configcommon.epoch.toi{1}              = [-5, 25];  
configcommon.epoch.toi{2}              = [-2, 1];  
configcommon.epoch.toi{3}              = [1, -2];  
configcommon.epoch.pad(1)              = 10;
configcommon.epoch.pad(2)              = 0.5;
configcommon.epoch.pad(3)              = 0.5;

configcommon.circus.reref              = 'no';
configcommon.circus.refchan            = '';
configcommon.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'dtx', 'SpykingCircus');
configcommon.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.suffix             = '-1';

configcommon.stats.bltoi{1}            = [-2, -1];
configcommon.stats.bltoi{1}            = [-2, -1];
configcommon.stats.bltoi{1}            = [0, 1];
configcommon.stats.actoi{1}            = [-1, 0];
configcommon.stats.actoi{1}            = [-1, 0];
configcommon.stats.actoi{1}            = [1, 0];
configcommon.stats.alpha               = 0.025;

configcommon.spike.slidestep           = [0.01];
configcommon.spike.toispikerate{1}     = [-0.1 0.1];           % for plotting spikerate
configcommon.spike.toispikerate{2}     = [-0.1 0.1];           % for plotting spikerate
configcommon.spike.toispikerate{3}     = [-0.1 0.1];           % for plotting spikerate
configcommon.spike.resamplefs          = 1000;
configcommon.spike.bltoi{1}            = [-2, -1];
configcommon.spike.bltoi{1}            = [-2, -1];
configcommon.spike.bltoi{1}            = [0, 1];

configcommon.spike.ISIbins             = [0:0.005:0.150];


%% Rodent 1
config{1}                           = configcommon;
config{1}.prefix                    = 'DTX5-';
config{1}.rawdir                    = fullfile(rootpath_data, 'DTX5-M1-10uM', '2019_03_19_DTX-5');
config{1}.imagesavedir              = fullfile(imagesavedir, 'DTX5');       % where to print images
config{1}.directorylist{1}          =  {'2019-03-19_13-18',...
                                        '2019-03-19_14-30',...
                                        '2019-03-19_16-30',...
                                        '2019-03-19_18-30',...
                                        '2019-03-19_20-30',...
                                        '2019-03-19_22-30',...
                                        '2019-03-20_00-30'};    
                                    
config{1}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{1}.labels.macro              = {'E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{1}.align.channel             = {'E12LFP.ncs'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.LFP.channel               = config{1}.labels.macro;
config{1}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 2
config{2}                           = configcommon;
config{2}.prefix                    = 'DTX2-';
config{2}.rawdir                    = fullfile(rootpath_data, 'DTX2-M1-10uM', '2019_03_01_DTX-2');
config{2}.imagesavedir              = fullfile(imagesavedir, 'DTX2');       % where to print images
config{2}.directorylist{1}          =  {'2019-03-01_12-14',...
                                        '2019-03-01_12-33',...
                                        '2019-03-01_14-14',...
                                        '2019-03-01_16-14',...
                                        '2019-03-01_18-14',...
                                        '2019-03-01_20-14'};    
                                    
config{2}.labels.micro              = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{2}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
        'ECoGS1','ECoGM1','ECoGPtA'}; % /!\ Mistake during acquisition : ECoGS1 = ECoGM1G};
    
config{2}.align.channel                = {'E13LFP.ncs'};
config{2}.LFP.channel                  = config{2}.labels.macro;
config{2}.circus.channel            = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 3
config{3}                           = configcommon;
config{3}.prefix                    = 'DTX4-';
config{3}.rawdir                    = fullfile(rootpath_data, 'DTX4-M1-10uM', '2019_03_08_DTX-4');
config{3}.imagesavedir              = fullfile(imagesavedir,'DTX4');       % where to print images
config{3}.directorylist{1}          =  {'2019-03-08_12-28',...
                                        '2019-03-08_14-05',...
                                        '2019-03-08_16-05',...
                                        '2019-03-08_18-05',...
                                        '2019-03-08_20-05'};

config{3}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{3}.labels.macro              = {'E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP'...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{3}.align.channel                = {'E13LFP.ncs'};
config{3}.LFP.channel                  = config{3}.labels.macro;
config{3}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};

%% Rodent 4
config{4}                           = configcommon;
config{4}.prefix                    = 'DTX10-';
config{4}.rawdir                    = fullfile(rootpath_data, 'DTX10-M1-10uM', '2019_03_28_DTX-10');
config{4}.imagesavedir              = fullfile(imagesavedir,'DTX10');       % where to print images
config{4}.directorylist{1}          =  {'2019-03-28_13-41',...
                                        '2019-03-28_14-50',...
                                        '2019-03-28_15-06',...
                                        '2019-03-28_17-06',...
                                        '2019-03-28_19-06'};
                                    
config{4}.labels.micro              = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{4}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{4}.align.channel                = {'E13LFP.ncs'};
config{4}.LFP.channel                  = config{4}.labels.macro;
config{4}.circus.channel            = {'E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 5
config{5}                           = configcommon;
config{5}.prefix                    = 'DTX7-';
config{5}.rawdir                    = fullfile(rootpath_data, 'DTX7-M1-10uM', '2019_03_22_DTX-7');
config{5}.imagesavedir              = fullfile(imagesavedir,'DTX7');       % where to print images
config{5}.directorylist{1}          =  {'2019-03-22_12-31',...
                                        '2019-03-22_14-31',...
                                        '2019-03-22_16-31',...
                                        '2019-03-22_18-31'};
                                    
config{5}.labels.micro              = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{5}.labels.macro              = {'E06LFP','E07LFP','E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{5}.align.channel                = {'E13LFP.ncs'};
config{5}.LFP.channel                  = config{5}.labels.macro;
config{5}.circus.channel            = {'E06','E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};


%% Rodent 6
config{6}                           = configcommon;
config{6}.prefix                    = 'DTX6-';
config{6}.rawdir                    = fullfile(rootpath_data, 'DTX6-M1-10uM', '2019_03_21_DTX-6');
config{6}.imagesavedir              = fullfile(imagesavedir,'DTX6');       % where to print images
config{6}.directorylist{1}          =  {'2019-03-21_14-12',...
                                        '2019-03-21_16-12',...
                                        '2019-03-21_18-12',...
                                        '2019-03-21_20-12'};
                                    
config{6}.labels.micro              = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
config{6}.labels.macro              = {'E08LFP','E09LFP','E10LFP','E11LFP','E12LFP','E13LFP','E14LFP','E15LFP','E16LFP',...
    'ECoGM1G','ECoGM1D','ECoGPtA'};

config{6}.align.channel                = {'E13LFP.ncs'};
config{6}.LFP.channel                  = config{6}.labels.macro;
config{6}.circus.channel            = {'E07','E08','E09','E10','E11','E12','E13','E14','E15','E16'};
end




