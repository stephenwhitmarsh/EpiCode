function [config] = pnh_setparams(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [config] = pnh_setparams(varargin)
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

pc = ispc;

% overwrite if os is forced
if nargin == 1
    if strcmp(varargin{1}, 'pc')
        pc = true;
    elseif strcmp(varargin{1}, 'unix')
        pc = false;
    else
        error('os not recognized');
    end
end

if ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\vn_pnh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
else
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/vn_pnh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
end

% Patient 1
config{1}.prefix                    = '2230-';
config{1}.rawdir                    = fullfile(rootpath_data,     'pat_02230_0674', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data',   'pnh');
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');

config{1}.muse.name                 = {'PSW', 'FA', 'ES', 'SEIZURE'};
config{1}.muse.startmarker.PSW      = "VF1__START__";
config{1}.muse.endmarker.PSW        = "VF1__START__";
config{1}.muse.startmarker.FA       = "RR__START__";
config{1}.muse.endmarker.FA         = "RR__START__";
config{1}.muse.startmarker.ES       = "P";
config{1}.muse.endmarker.ES         = "P";
config{1}.muse.startmarker.SEIZURE  = "CriseStart";
config{1}.muse.endmarker.SEIZURE    = "CriseEnd";

config{1}.align.name                = {'PSW', 'FA', 'ES'};
config{1}.align.demean              = 'yes';
config{1}.align.channel             = {'m1pNs_1', 'm1pNs_4', 'm1pNs_6', 'm1pNs_7'};
config{1}.align.reref               = 'no';
config{1}.align.refmethod           = 'bipolar';
config{1}.align.latency.PSW         = [-0.2 2];
config{1}.align.latency.FA          = [-0.2 1];
config{1}.align.latency.ES          = [-0.1 0.1];

config{1}.LFP.name                  = {'PSW','FA','ES','SEIZURE'};
config{1}.LFP.channel               = {'m1pNs_1', 'm1pNs_4', 'm1pNs_6', 'm1pNs_7'};
config{1}.LFP.resamplefs            = 500; % 32000/500=64 which is nice
config{1}.LFP.dftfilter             = 'yes'; % notch filter

config{1}.epoch.toi.PSW             = [-2, 2];
config{1}.epoch.toi.FA              = [-2, 2];
config{1}.epoch.toi.ES              = [-1, 1];
config{1}.epoch.toi.SEIZURE         = [-5, 13];
config{1}.epoch.pad.PSW             = 0.5;
config{1}.epoch.pad.ES              = 0.5;
config{1}.epoch.pad.FA              = 0.5;
config{1}.epoch.pad.SEIZURE         = 0.5;

config{1}.TFR.keeptrials            = 'no';
config{1}.TFR.foi.PSW               = 1:200;
config{1}.TFR.foi.FA                = 1:200;
config{1}.TFR.foi.ES                = 1:200;
config{1}.TFR.foi.SEIZURE           = 1:200;
config{1}.TFR.t_ftimwin.PSW         = 10 ./ config{1}.TFR.foi.PSW;
config{1}.TFR.t_ftimwin.FA          = 10 ./ config{1}.TFR.foi.FA;
config{1}.TFR.t_ftimwin.ES          = 10 ./ config{1}.TFR.foi.ES;
config{1}.TFR.t_ftimwin.SEIZURE     = 10 ./ config{1}.TFR.foi.ES;
config{1}.TFR.toi.PSW               = -2 : 0.005 : 2;
config{1}.TFR.toi.FA                = -2 : 0.005 : 2;
config{1}.TFR.toi.ES                = -1 : 0.005 : 1;
config{1}.TFR.toi.SEIZURE           = -5 : 0.005 : 13;
config{1}.TFR.bl.PSW                = [-1.5, -1];
config{1}.TFR.bl.FA                 = [-1.5, -1];
config{1}.TFR.bl.ES                 = [-1, -0.5];
config{1}.TFR.bl.SEIZURE            = [-5, -0.5];

config{1}.circus.channel                        = {'m1pNs_1','m1pNs_2','m1pNs_6','m1pNs_7','m1pNs_8'};
config{1}.circus.reref                          = 'no';
config{1}.circus.refchan                        = '';
config{1}.circus.outputdir                      = 'SpykingCircus';
config{1}.circus.paramfile                      = fullfile(rootpath_analysis, 'data', 'pnh', 'SpykingCircus.params');
config{1}.circus.params.detection.spike_thresh  = '6';
config{1}.circus.params.filtering.cut_off       = '300, auto';
config{1}.circus.params.filtering.remove_median = 'False';
config{1}.circus.params.clustering.max_elts     = '20000';
config{1}.circus.params.detection.peaks         = 'negative';
config{1}.circus.params.data.stream_mode        = 'mapping-file';
config{1}.circus.params.data.mapping_file       = 'filelist.txt';

config{1}.spike.name                = {'PSW', 'FA', 'ES', 'SEIZURE'};
config{1}.spike.slidestep           = [0.01, 0.01, 0.001];
config{1}.spike.toi.PSW             = [-2, 2];           % for plotting spikerate
config{1}.spike.toi.FA              = [-2, 2];           % for p200mslotting spikerate
config{1}.spike.toi.ES              = [-1, 1];          % for plotting spikerate
config{1}.spike.toi.SEIZURE         = [-5, 13];          % for plotting spikerate
config{1}.spike.bl.PSW              = [-1, -0.5];
config{1}.spike.bl.FA               = [-1, -0.5];
config{1}.spike.bl.ES               = [-1, -0.5];
config{1}.spike.bl.SEIZURE          = [-5, -1];
config{1}.spike.pad.PSW             = 0.1;
config{1}.spike.pad.FA              = 0.1;
config{1}.spike.pad.ES              = 0.1;
config{1}.spike.pad.SEIZURE         = 0.1;
config{1}.spike.resamplefs.PSW      = 1000;
config{1}.spike.resamplefs.FA       = 1000;
config{1}.spike.resamplefs.ES       = 1000;
config{1}.spike.resamplefs.SEIZURE  = 1000;
config{1}.spike.pre                 = 0.001;
config{1}.spike.post                = 0.002;
% config{1}.spike.baseline            = [-0.001 -0.0005];
config{1}.spike.ISIbins             = 0 : 0.0005 : 0.150;
config{1}.spike.nrsdfbins           = 100;
% config{1}.spike.psthbin.PSW
% config{1}.spike.psthbin.FA
% config{1}.spike.psthbin.ES
config{1}.spike.sdftimwin.PSW       = [-0.05 0.05];
config{1}.spike.sdftimwin.FA        = [-0.05 0.05];
config{1}.spike.sdftimwin.ES        = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.sdftimwin.SEIZURE   = [-0.01 0.01]; % depends a lot on pattern, default is 0.05
config{1}.spike.psthbin.PSW         = 0.1; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.FA          = 0.1; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.ES          = 0.02; % depends a lot on pattern, default is too large
config{1}.spike.psthbin.SEIZURE     = 0.05; % depends a lot on pattern, default is too large

config{1}.stats.toi.PSW             = [-0.5, 2];
config{1}.stats.toi.FA              = [-0.5, 2];
config{1}.stats.toi.ES              = [-0.5, 1];
config{1}.stats.toi.SEIZURE         = [-0.5, 13];
config{1}.stats.bl.PSW              = [-2, -1];
config{1}.stats.bl.FA               = [-2, -1];
config{1}.stats.bl.ES               = [-1, -0.5];
config{1}.stats.bl.SEIZURE          = [-5, -1];
config{1}.stats.alpha               = 0.025;

config{1}.interval.histbin.ES       = 0.25;
config{1}.interval.histbin.FA       = 1;
config{1}.interval.histbin.PSW      = 1;
config{1}.interval.histbin.SEIZURE  = 1;
config{1}.interval.histlim.ES       = 15;
config{1}.interval.histlim.FA       = 60;
config{1}.interval.histlim.PSW      = 60;
config{1}.interval.histlim.SEIZURE  = 60;

config{1}.plot.dir.PSW              = [1,  2,  3,  4,  5,  6];
config{1}.plot.trial.PSW            = [10, 10, 10, 10, 10, 10];
config{1}.plot.unit.PSW             = [0,  0,  0,  0,  0,  0];
config{1}.plot.dir.FA               = [1,  2,  3,  4,  5,  6];
config{1}.plot.trial.FA             = [10, 10, 10, 10, 10, 10];
config{1}.plot.unit.FA              = [0,  0,  0,  0,  0,  0];
config{1}.plot.dir.ES               = [1,  2,  3,  4,  5,  6];
config{1}.plot.trial.ES             = [10, 10, 10, 10, 10, 10];
config{1}.plot.unit.ES              = [0,  0,  0,  0,  0,  0];
config{1}.plot.dir.SEIZURE          = [3,  4,  4,  4,  9, 11];
config{1}.plot.trial.SEIZURE        = [1,  1,  2,  3,  1, 1];
config{1}.plot.unit.SEIZURE         = [0,  0,  0,  0,  0, 0]; 

% Patient 2, first nodule
config{2}                           = config{1};
config{2}.prefix                    = '2614L-';
config{2}.rawdir                    = fullfile(rootpath_data, 'pat_02614_1073', 'eeg');
config{2}.align.channel             = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_5','mTNmi_6','mTNmi_8'};
config{2}.align.lpfilter            = 'yes';
config{2}.align.lpfreq              = 40;

config{2}.epoch.toi.FA              = [-2, 3];
config{2}.epoch.toi.ES              = [-1, 1];

config{2}.LFP.name                  = {'FA','ES'};
config{2}.LFP.channel               = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_5','mTNmi_6','mTNmi_7','mTNmi_8'};

config{2}.TFR.toi.FA                = -2 : 0.005 : 3;
config{2}.TFR.toi.ES                = -1 : 0.005 : 1;

config{2}.muse.name                 = {'FA','ES'};
config{2}.muse.startmarker.FA       = "p02614_07__START__";
config{2}.muse.endmarker.FA         = "p02614_07__START__";
config{2}.muse.startmarker.ES       = "p02614_09";
config{2}.muse.endmarker.ES         = "p02614_09";
config{2}.circus.channel            = {'mTNmi_2','mTNmi_3','mTNmi_5','mTNmi_6','mTNmi_7','mTNmi_8'}; % redo within 1st and 4th electrode - now removed
config{2}.circus.params.detection.spike_thresh  = '5.5';
config{3}.circus.params.filtering.remove_median = 'False';
config{2}.spike.name                            = {'FA','ES'};
config{2}.spike.toi.FA              = [-2, 3];           % for p200mslotting spikerate
config{2}.stats.toi.FA              = [-0.5, 3];

% Patient 2, second nodule
config{3}                           = config{1};
config{3}.prefix                    = '2614R-';
config{3}.rawdir                    = fullfile(rootpath_data, 'pat_02614_1073', 'eeg');
config{3}.LFP.channel               = {'mCasd_2','mCasd_1','mCasd_4','mCasd_6','mCasd_7'};
config{3}.align.channel             = {'mCasd_2','mCasd_1','mCasd_4','mCasd_6','mCasd_7'};
config{3}.align.lpfilter            = 'yes';
config{3}.align.lpfreq              = 40;

config{3}.epoch.toi.PSW             = [-2, 2];
config{3}.epoch.toi.FA              = [-2, 2];
config{3}.epoch.toi.ES              = [-1, 1];

config{3}.TFR.toi.PSW               = -2 : 0.005 : 2;
config{3}.TFR.toi.FA                = -2 : 0.005 : 2;
config{3}.TFR.toi.ES                = -1 : 0.005 : 1;

config{3}.LFP.name                  = {'PSW','FA','ES'};
config{3}.muse.name                 = {'PSW','FA','ES'};
config{3}.muse.startmarker.PSW      = "p02614_01__START__";
config{3}.muse.endmarker.PSW        = "p02614_01__START__";
config{3}.muse.startmarker.FA       = "p02614_04__START__";
config{3}.muse.endmarker.FA         = "p02614_04__START__";
config{3}.muse.startmarker.ES       = "p02614_02";
config{3}.muse.endmarker.ES         = "p02614_02";
config{3}.circus.params.filtering.remove_median = 'False';
config{3}.circus.channel            = {'mCasd_1','mCasd_2','mCasd_3','mCasd_4','mCasd_5','mCasd_6','mCasd_7'};
config{3}.spike.name                = {'PSW','FA','ES'};


% Patient 3 %
config{4}                           = config{1};
config{4}.prefix                    = '2689-';
config{4}.rawdir                    = fullfile(rootpath_data, 'pat_02689_1168', 'eeg');
config{4}.align.channel             = {'mLMI1_3','mLMI1_2','mLMI1_4','mLMI1_5','mLMI1_7'};
config{4}.align.lpfilter            = 'yes';
config{4}.align.lpfreq              = 40;
config{4}.align.latency.PSW         = [-0.2, 5];

config{4}.epoch.toi.PSW             = [-2, 5];
config{4}.epoch.toi.FA              = [-2, 3];
config{4}.epoch.toi.ES              = [-1, 1];

config{4}.TFR.toi.PSW               = [-2 : 0.005 : 5];
config{4}.TFR.toi.FA                = -2 : 0.005 : 3;
config{4}.TFR.toi.ES                = -1 : 0.005 : 1;

config{4}.LFP.name                  = {'PSW','FA','ES'};
config{4}.LFP.channel               = {'mLMI1_3','mLMI1_2','mLMI1_4','mLMI1_5','mLMI1_7'};
config{4}.muse.name                 = {'PSW','FA','ES'};
config{4}.muse.startmarker.PSW      = "Marker1__START__";
config{4}.muse.endmarker.PSW        = "Marker1__START__";
config{4}.muse.startmarker.FA       = "Marker3__START__";
config{4}.muse.endmarker.FA         = "Marker3__START__";
config{4}.muse.startmarker.ES       = "Marker2";
config{4}.muse.endmarker.ES         = "Marker2";
config{4}.circus.channel            = {'mLMI1_2','mLMI1_3','mLMI1_4','mLMI1_7'};
config{4}.spike.name                = {'PSW','FA','ES'};
config{4}.spike.toi.PSW             = [-2, 5];           % for p200mslotting spikerate
config{4}.stats.toi.PSW             = [-0.5, 5];
config{4}.spike.toi.FA              = [-2, 3];           % for p200mslotting spikerate
config{4}.stats.toi.FA              = [-0.5, 3];

% DATA

config{1}.directorylist             = [];
config{1}.directorylist{1}          =  {'02230_2015-02-25_12-36'...
    '02230_2015-02-25_14-36'...
    '02230_2015-02-25_15-16'...
    '02230_2015-02-25_17-16'...
    '02230_2015-02-25_19-16'...
    '02230_2015-02-25_21-16'...
    '02230_2015-02-25_23-16'...
    '02230_2015-02-26_01-16'...
    '02230_2015-02-26_03-16'...
    '02230_2015-02-26_05-16'...
    '02230_2015-02-26_07-16'};

%                                         '02230_2015-02-26_09-16'...
%                                         '02230_2015-02-26_10-03'...
%                                         '02230_2015-02-26_10-31'};

config{2}.directorylist             = [];
config{2}.directorylist{1}          =  {'02614_2018-06-12_15-23'...
    '02614_2018-06-12_17-23'...
    '02614_2018-06-12_21-23'};

config{3}.directorylist             = [];
config{3}.directorylist{1}          =  {'02614_2018-06-12_15-23'...
    '02614_2018-06-12_17-23'...
    '02614_2018-06-12_21-23'};

config{4}.directorylist             = [];
config{4}.directorylist{1}          =  {'02689_2019-02-14_16-25'...
    '02689_2019-02-14_18-25'...
    '02689_2019-02-14_20-25'...
    '02689_2019-02-14_22-25'...
    '02689_2019-02-15_00-25'...
    '02689_2019-02-15_02-25'...
    '02689_2019-02-15_04-25'...
    '02689_2019-02-15_06-25'...
    '02689_2019-02-15_08-25'};


muscale = 50;

config{1}.plot.toi.patternA               = [1500.8, 1503.3, 0];
config{1}.plot.hpfilter.patternA{1}       = 'no';
config{1}.plot.hpfilter.patternA{2}       = 'no';
config{1}.plot.hpfilter.patternA{3}       = 'no';
config{1}.plot.hpfilter.patternA{4}       = 'no';
config{1}.plot.hpfilter.patternA{5}       = 'yes';
config{1}.plot.hpfilter.patternA{6}       = 'yes';
config{1}.plot.hpfreq.patternA{1}         = [];
config{1}.plot.hpfreq.patternA{2}         = [];
config{1}.plot.hpfreq.patternA{3}         = [];
config{1}.plot.hpfreq.patternA{4}         = [];
config{1}.plot.hpfreq.patternA{5}         = 300;
config{1}.plot.hpfreq.patternA{6}         = 300;
config{1}.plot.lpfilter.patternA{1}       = 'yes';
config{1}.plot.lpfilter.patternA{2}       = 'yes';
config{1}.plot.lpfilter.patternA{3}       = 'yes';
config{1}.plot.lpfilter.patternA{4}       = 'yes';
config{1}.plot.lpfilter.patternA{5}       = 'yes';
config{1}.plot.lpfilter.patternA{6}       = 'yes';
config{1}.plot.lpfreq.patternA{1}         = 300;
config{1}.plot.lpfreq.patternA{2}         = 300;
config{1}.plot.lpfreq.patternA{3}         = 300;
config{1}.plot.lpfreq.patternA{4}         = 300;
config{1}.plot.lpfreq.patternA{5}         = 3000;
config{1}.plot.lpfreq.patternA{6}         = 3000;
config{1}.plot.fname.patternA{1}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_2.ncs');
config{1}.plot.fname.patternA{2}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_1.ncs');
config{1}.plot.fname.patternA{3}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.fname.patternA{4}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_6.ncs');
config{1}.plot.fname.patternA{5}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_1.ncs');
config{1}.plot.fname.patternA{6}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_2.ncs');
config{1}.plot.refname.patternA{1}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_3.ncs');
config{1}.plot.refname.patternA{2}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_2.ncs');
config{1}.plot.refname.patternA{3}        = [];
config{1}.plot.refname.patternA{4}        = [];
config{1}.plot.refname.patternA{5}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.refname.patternA{6}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.scale.patternA{1}          = [-200 200];
config{1}.plot.scale.patternA{2}          = [-200 200];
config{1}.plot.scale.patternA{3}          = [-200 200];
config{1}.plot.scale.patternA{4}          = [-200 200];
config{1}.plot.scale.patternA{5}          = [-muscale muscale];
config{1}.plot.scale.patternA{6}          = [-muscale muscale];

config{1}.plot.toi.patternB               = [1976.2-0.5, 1976.2+2, -0.5];
config{1}.plot.hpfilter.patternB{1}       = 'no';
config{1}.plot.hpfilter.patternB{2}       = 'no';
config{1}.plot.hpfilter.patternB{3}       = 'no';
config{1}.plot.hpfilter.patternB{4}       = 'no';
config{1}.plot.hpfilter.patternB{5}       = 'yes';
config{1}.plot.hpfilter.patternB{6}       = 'yes';
config{1}.plot.hpfreq.patternB{1}         = [];
config{1}.plot.hpfreq.patternB{2}         = [];
config{1}.plot.hpfreq.patternB{3}         = [];
config{1}.plot.hpfreq.patternB{4}         = [];
config{1}.plot.hpfreq.patternB{5}         = 300;
config{1}.plot.hpfreq.patternB{6}         = 300;
config{1}.plot.lpfilter.patternB{1}       = 'yes';
config{1}.plot.lpfilter.patternB{2}       = 'yes';
config{1}.plot.lpfilter.patternB{3}       = 'yes';
config{1}.plot.lpfilter.patternB{4}       = 'yes';
config{1}.plot.lpfilter.patternB{5}       = 'yes';
config{1}.plot.lpfilter.patternB{6}       = 'yes';
config{1}.plot.lpfreq.patternB{1}         = 300;
config{1}.plot.lpfreq.patternB{2}         = 300;
config{1}.plot.lpfreq.patternB{3}         = 300;
config{1}.plot.lpfreq.patternB{4}         = 300;
config{1}.plot.lpfreq.patternB{5}         = 3000;
config{1}.plot.lpfreq.patternB{6}         = 3000;
config{1}.plot.fname.patternB{1}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_2.ncs');
config{1}.plot.fname.patternB{2}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_1.ncs');
config{1}.plot.fname.patternB{3}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.fname.patternB{4}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_6.ncs');
config{1}.plot.fname.patternB{5}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_1.ncs');
config{1}.plot.fname.patternB{6}          = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_8.ncs');
config{1}.plot.refname.patternB{1}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_3.ncs');
config{1}.plot.refname.patternB{2}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_1pNs_2.ncs');
config{1}.plot.refname.patternB{3}        = [];
config{1}.plot.refname.patternB{4}        = [];
config{1}.plot.refname.patternB{5}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.refname.patternB{6}        = fullfile(config{1}.rawdir, '02230_2015-02-25_17-16', '02230_2015-02-25_17-16_m1pNs_4.ncs');
config{1}.plot.scale.patternB{1}          = [-200 200];
config{1}.plot.scale.patternB{2}          = [-200 200];
config{1}.plot.scale.patternB{3}          = [-200 200];
config{1}.plot.scale.patternB{4}          = [-200 200];
config{1}.plot.scale.patternB{5}          = [-muscale muscale];
config{1}.plot.scale.patternB{6}          = [-muscale muscale];

config{1}.plot.toi.patternC               = [894.659-0.4, 894.659+0.4, -0.4];
config{1}.plot.hpfilter.patternC{1}       = 'no';
config{1}.plot.hpfilter.patternC{2}       = 'no';
config{1}.plot.hpfilter.patternC{3}       = 'no';
config{1}.plot.hpfilter.patternC{4}       = 'no';
config{1}.plot.hpfilter.patternC{5}       = 'yes';
config{1}.plot.hpfilter.patternC{6}       = 'yes';
config{1}.plot.hpfreq.patternC{1}         = [];
config{1}.plot.hpfreq.patternC{2}         = [];
config{1}.plot.hpfreq.patternC{3}         = [];
config{1}.plot.hpfreq.patternC{4}         = [];
config{1}.plot.hpfreq.patternC{5}         = 300;
config{1}.plot.hpfreq.patternC{6}         = 300;
config{1}.plot.lpfilter.patternC{1}       = 'yes';
config{1}.plot.lpfilter.patternC{2}       = 'yes';
config{1}.plot.lpfilter.patternC{3}       = 'yes';
config{1}.plot.lpfilter.patternC{4}       = 'yes';
config{1}.plot.lpfilter.patternC{5}       = 'yes';
config{1}.plot.lpfilter.patternC{6}       = 'yes';
config{1}.plot.lpfreq.patternC{1}         = 300;
config{1}.plot.lpfreq.patternC{2}         = 300;
config{1}.plot.lpfreq.patternC{3}         = 300;
config{1}.plot.lpfreq.patternC{4}         = 300;
config{1}.plot.lpfreq.patternC{5}         = 3000;
config{1}.plot.lpfreq.patternC{6}         = 3000;
config{1}.plot.fname.patternC{1}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_1pNs_2.ncs');
config{1}.plot.fname.patternC{2}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_1pNs_1.ncs');
config{1}.plot.fname.patternC{3}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_4.ncs');
config{1}.plot.fname.patternC{4}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_6.ncs');
config{1}.plot.fname.patternC{5}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_1.ncs');
config{1}.plot.fname.patternC{6}          = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_8.ncs');
config{1}.plot.refname.patternC{1}        = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_1pNs_3.ncs');
config{1}.plot.refname.patternC{2}        = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_1pNs_2.ncs');
config{1}.plot.refname.patternC{3}        = [];
config{1}.plot.refname.patternC{4}        = [];
config{1}.plot.refname.patternC{5}        = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_4.ncs');
config{1}.plot.refname.patternC{6}        = fullfile(config{1}.rawdir, '02230_2015-02-25_14-36', '02230_2015-02-25_14-36_m1pNs_4.ncs');
config{1}.plot.scale.patternC{1}          = [-200 200];
config{1}.plot.scale.patternC{2}          = [-200 200];
config{1}.plot.scale.patternC{3}          = [-200 200];
config{1}.plot.scale.patternC{4}          = [-200 200];
config{1}.plot.scale.patternC{5}          = [-muscale muscale];
config{1}.plot.scale.patternC{6}          = [-muscale muscale];

%% NODULE 2
config{2}.plot                            = config{1}.plot;
config{2}.plot.toi.patternA               = [1004.0, 1008.0, 0];
config{2}.plot.fname.patternA{1}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_2.ncs');
config{2}.plot.fname.patternA{2}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_1.ncs');
config{2}.plot.fname.patternA{3}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_2.ncs');
config{2}.plot.fname.patternA{4}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_5.ncs');
config{2}.plot.fname.patternA{5}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_2.ncs');
config{2}.plot.fname.patternA{6}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_6.ncs');
config{2}.plot.refname.patternA{1}        = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_3.ncs');
config{2}.plot.refname.patternA{2}        = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_2.ncs');

config{2}.plot.toi.patternB               = [412.9, 415.4, 0];
config{2}.plot.fname.patternB{1}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_2.ncs');
config{2}.plot.fname.patternB{2}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_1.ncs');
config{2}.plot.fname.patternB{3}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_2.ncs');
config{2}.plot.fname.patternB{4}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_5.ncs');
config{2}.plot.fname.patternB{5}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_2.ncs');
config{2}.plot.fname.patternB{6}          = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_mCasd_4.ncs');
config{2}.plot.refname.patternB{1}        = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_3.ncs');
config{2}.plot.refname.patternB{2}        = fullfile(config{2}.rawdir, '02614_2018-06-12_17-23', '02614_2018-06-12_17-23_Casd_2.ncs');

config{2}.plot.toi.patternC               = [1664.850-0.20, 1664.850+0.4, -0.20];
config{2}.plot.fname.patternC{1}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_Casd_2.ncs');
config{2}.plot.fname.patternC{2}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_Casd_1.ncs');
config{2}.plot.fname.patternC{3}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mCasd_2.ncs');
config{2}.plot.fname.patternC{4}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mCasd_5.ncs');
config{2}.plot.fname.patternC{5}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mCasd_2.ncs');
config{2}.plot.fname.patternC{6}          = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mCasd_5.ncs');
config{2}.plot.refname.patternC{1}        = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_Casd_3.ncs');
config{2}.plot.refname.patternC{2}        = fullfile(config{2}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_Casd_2.ncs');

%% nodule 3
config{3}.plot                            = config{1}.plot;
config{3}.plot.toi                        = rmfield(config{3}.plot.toi, 'patternC');
config{3}.plot.toi.patternA               = [1181, 1185, 0];
config{3}.plot.fname.patternA{1}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_2.ncs');
config{3}.plot.fname.patternA{2}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_1.ncs');
config{3}.plot.fname.patternA{3}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_3.ncs');
config{3}.plot.fname.patternA{4}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_5.ncs');
config{3}.plot.fname.patternA{5}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_3.ncs');
config{3}.plot.fname.patternA{6}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_7.ncs');
config{3}.plot.refname.patternA{1}        = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_3.ncs');
config{3}.plot.refname.patternA{2}        = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_2.ncs');

config{3}.plot.toi.patternB               = [722.5, 723.3, -0.4];
config{3}.plot.fname.patternB{1}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_2.ncs');
config{3}.plot.fname.patternB{2}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_1.ncs');
config{3}.plot.fname.patternB{3}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_3.ncs');
config{3}.plot.fname.patternB{4}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_4.ncs');
config{3}.plot.fname.patternB{5}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_3.ncs');
config{3}.plot.fname.patternB{6}          = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_mTNmi_5.ncs');
config{3}.plot.refname.patternB{1}        = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_3.ncs');
config{3}.plot.refname.patternB{2}        = fullfile(config{3}.rawdir, '02614_2018-06-12_15-23', '02614_2018-06-12_15-23_TNmi_2.ncs');

%% nodule 3
config{4}.plot                            = config{1}.plot;
config{4}.plot.toi.patternA               = [919.979, 927.535, 0];
config{4}.plot.fname.patternA{1}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_LMI1_2.ncs');
config{4}.plot.fname.patternA{2}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_LMI1_1.ncs');
config{4}.plot.fname.patternA{3}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_mLMI1_3.ncs');
config{4}.plot.fname.patternA{4}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_mLMI1_4.ncs');
config{4}.plot.fname.patternA{5}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_mLMI1_3.ncs');
config{4}.plot.fname.patternA{6}          = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_mLMI1_4.ncs');
config{4}.plot.refname.patternA{1}        = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_LMI1_3.ncs');
config{4}.plot.refname.patternA{2}        = fullfile(config{4}.rawdir, '02689_2019-02-15_02-25', '02689_2019-02-15_02-25_LMI1_2.ncs');

config{4}.plot.toi.patternB               = [648.952-0.4, 648.952+1, -0.4];
config{4}.plot.fname.patternB{1}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_2.ncs');
config{4}.plot.fname.patternB{2}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_1.ncs');
config{4}.plot.fname.patternB{3}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_3.ncs');
config{4}.plot.fname.patternB{4}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_4.ncs');
config{4}.plot.fname.patternB{5}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_3.ncs');
config{4}.plot.fname.patternB{6}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_4.ncs');
config{4}.plot.refname.patternB{1}        = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_3.ncs');
config{4}.plot.refname.patternB{2}        = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_2.ncs');

config{4}.plot.toi.patternC               = [26.9326-0.4, 26.9326+0.6, -0.4];
config{4}.plot.fname.patternC{1}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_2.ncs');
config{4}.plot.fname.patternC{2}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_1.ncs');
config{4}.plot.fname.patternC{3}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_3.ncs');
config{4}.plot.fname.patternC{4}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_4.ncs');
config{4}.plot.fname.patternC{5}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_3.ncs');
config{4}.plot.fname.patternC{6}          = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_mLMI1_4.ncs');
config{4}.plot.refname.patternC{1}        = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_3.ncs');
config{4}.plot.refname.patternC{2}        = fullfile(config{4}.rawdir, '02689_2019-02-14_16-25', '02689_2019-02-14_16-25_LMI1_2.ncs');
    
% 




