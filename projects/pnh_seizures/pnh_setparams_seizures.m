function [config] = pnh_setparams_seizures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [config] = pnh_setparams
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

if ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
else       
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
end
    
if ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\vn_pnh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
else       
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/vn_pnh';
    rootpath_data       = '/network/lustre/iss01/epimicro/patients/raw';
end
    
%% 

% Patient 1, perivetricular heterotopia #1
disp('setting parameters');

% patient 1
muscale = 50; % multiunit scale

% Patient 1, perivetricular heterotopia #1
config{1}.prefix                        = '2230-seizures-';
config{1}.rawdir                        = fullfile(rootpath_data,     'pat_02230_0674', 'eeg');
config{1}.datasavedir                   = fullfile(rootpath_analysis, 'data');         % where to write data
config{1}.imagesavedir                  = fullfile(rootpath_analysis, 'images');       % where to print images

config{1}.name                          = {'seizure'};
config{1}.muse.startmarker.seizure      = "CriseStart";
config{1}.muse.endmarker.seizure        = "CriseStart";
config{1}.muse.backupdir                = fullfile(rootpath_analysis, 'markerbackup');

config{1}.LFP.name                      = 'seizure';
config{1}.LFP.hpfilter                  = 'no';
config{1}.LFP.hpfreq                    = 1;
config{1}.LFP.resamplefs                = 1000;
config{1}.LFP.rerefmethod               = 'bipolar';
config{1}.LFP.baseline                  = 'yes';
config{1}.LFP.baselinewindow.seizure    = [-1, -0.5];
config{1}.LFP.slidestep                 = [0.01, 0.01, 0.01];
config{1}.LFP.channel                   = {'_2pNi_1', '_1pNs_1', '_1pHe_1'};

config{1}.LFP.reref                     = 'yes';
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
% config{1}.refmethod                     = 'bipolar';


config{1}.epoch.toi.seizure             = [-1  10];
config{1}.epoch.pad.seizure             = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patient 2, first nodule %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

config{2}                           = config{1};
config{2}.prefix                    = '2614-seizures-';
config{2}.rawdir                    = fullfile(rootpath_data, 'pat_02614_1073', 'eeg');
config{2}.LFP.channel               = {'_TNmg_1', '_TNmi_1'};

%%%%%%%%%%%%%
% Patient 3 %
%%%%%%%%%%%%%

config{3}                           = config{1};
config{3}.prefix                    = '2689-seizures-';
config{3}.rawdir                    = fullfile(rootpath_data, 'pat_02689_1168', 'eeg');
config{3}.LFP.channel               = {'_TNmg_1', '_TNmi_1'};


