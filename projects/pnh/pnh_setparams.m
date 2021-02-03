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
config{1}.muse.name                 = {'PSW','FA','ES'};
config{1}.muse.startmarker.PSW      = "VF1__START__";
config{1}.muse.endmarker.PSW        = "VF1__START__";
config{1}.muse.startmarker.FA       = "RR__START__";
config{1}.muse.endmarker.FA         = "RR__START__";
config{1}.muse.startmarker.ES       = "P";
config{1}.muse.endmarker.ES         = "P";
config{1}.align.name                = {'PSW','FA','ES'};
config{1}.align.demean              = 'yes';
config{1}.align.channel             = {'m1pNs_1', 'm1pNs_4', 'm1pNs_6', 'm1pNs_7'};
config{1}.align.reref               = 'no';
config{1}.align.refmethod           = 'bipolar';
config{1}.align.latency.PSW         = [-0.2 2];
config{1}.align.latency.FA          = [-0.2 1];
config{1}.align.latency.ES          = [-0.1 0.1];
config{1}.LFP.name                  = {'PSW','FA','ES'};
config{1}.LFP.channel               = {'m1pNs_1', 'm1pNs_4', 'm1pNs_6', 'm1pNs_7'};
config{1}.LFP.resamplefs            = 500; % 32000/500=64 which is nice
config{1}.LFP.dftfilter             = 'yes'; % notch filter
config{1}.epoch.toi.PSW             = [-2  2];
config{1}.epoch.pad.PSW             = 0.5;
config{1}.epoch.toi.FA              = [-2  2];
config{1}.epoch.pad.FA              = 0.5;
config{1}.epoch.toi.ES              = [-1  1];
config{1}.epoch.pad.ES              = 0.5;
config{1}.TFR.keeptrials            = 'no';
config{1}.TFR.foi.PSW               = [10:200];
config{1}.TFR.foi.FA                = [10:200];
config{1}.TFR.foi.ES                = [10:200];
config{1}.TFR.t_ftimwin.PSW         = 7 ./ config{1}.TFR.foi.PSW;
config{1}.TFR.t_ftimwin.FA          = 7 ./ config{1}.TFR.foi.FA;
config{1}.TFR.t_ftimwin.ES          = 7 ./ config{1}.TFR.foi.ES;
config{1}.TFR.toi.PSW               = [-2 : 0.005 : 2];
config{1}.TFR.toi.FA                = [-2 : 0.005 : 2];
config{1}.TFR.toi.ES                = [-1 : 0.005 : 1];
config{1}.TFR.bl.PSW               =  [-1.5, -1];
config{1}.TFR.bl.FA                =  [-1.5, -1];
config{1}.TFR.bl.ES                =  [-1, -0.5];

% Patient 2, first nodule
config{2}                           = config{1};
config{2}.prefix                    = '2614L-';
config{2}.rawdir                    = fullfile(rootpath_data, 'pat_02614_1073', 'eeg');
config{2}.align.channel             = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_5','mTNmi_6','mTNmi_8'};
config{2}.align.lpfilter            = 'yes';
config{2}.align.lpfreq              = 40;
config{2}.LFP.name                  = {'FA','ES'};
config{2}.LFP.channel               = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_5','mTNmi_6','mTNmi_7','mTNmi_8'};
config{2}.muse.name                 = {'FA','ES'};
config{2}.muse.startmarker.FA       = "p02614_07__START__";
config{2}.muse.endmarker.FA         = "p02614_07__START__";
config{2}.muse.startmarker.ES       = "p02614_09";
config{2}.muse.endmarker.ES         = "p02614_09";


% Patient 2, second nodule
config{3}                           = config{1};
config{3}.prefix                    = '2614R-';
config{3}.rawdir                    = fullfile(rootpath_data, 'pat_02614_1073', 'eeg');
config{3}.LFP.channel               = {'mCasd_2','mCasd_1','mCasd_4','mCasd_6','mCasd_7'};
config{3}.LFP.channel               = {'mCasd_2','mCasd_1','mCasd_4','mCasd_6','mCasd_7'};
config{3}.align.channel             = {'mCasd_2','mCasd_1','mCasd_4','mCasd_6','mCasd_7'};
config{3}.align.lpfilter            = 'yes';
config{3}.align.lpfreq              = 40;
config{3}.LFP.name                  = {'PSW','FA','ES'};
config{3}.muse.name                 = {'PSW','FA','ES'};
config{3}.muse.startmarker.PSW      = "p02614_01__START__";
config{3}.muse.endmarker.PSW        = "p02614_01__START__";
config{3}.muse.startmarker.FA       = "p02614_04__START__";
config{3}.muse.endmarker.FA         = "p02614_04__START__";
config{3}.muse.startmarker.ES       = "p02614_02";
config{3}.muse.endmarker.ES         = "p02614_02";
                              
                                      
% Patient 3 %
config{4}                           = config{1};
config{4}.prefix                    = '2689-';
config{4}.rawdir                    = fullfile(rootpath_data, 'pat_02689_1168', 'eeg');
config{4}.align.channel             = {'mLMI1_3','mLMI1_2','mLMI1_4','mLMI1_5','mLMI1_7'};
config{4}.align.lpfilter            = 'yes';
config{4}.align.lpfreq              = 40;
config{4}.LFP.channel               = {'mLMI1_3','mLMI1_2','mLMI1_4','mLMI1_5','mLMI1_7'};
config{4}.muse.name                 = {'PSW','FA','ES'};
config{4}.muse.startmarker.PSW      = "Marker1__START__";
config{4}.muse.endmarker.PSW        = "Marker1__START__";
config{4}.muse.startmarker.FA       = "Marker3__START__";
config{4}.muse.endmarker.FA         = "Marker3__START__";
config{4}.muse.startmarker.ES       = "Marker2";
config{4}.muse.endmarker.ES         = "Marker2";
                              
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
                                        '02230_2015-02-26_07-16'...
                                        '02230_2015-02-26_09-16'...
                                        '02230_2015-02-26_10-03'...
                                        '02230_2015-02-26_10-31'};
     
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
                                        '02689_2019-02-15_02-25'};
                                                                        