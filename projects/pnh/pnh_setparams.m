function [config] = pnh_setparams

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

%% Patient 1

config{1}.prefix                    = '2711-';
config{1}.rawdir                    = fullfile(rootpath_data,     'pat_02711_1193', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'hspike');
config{1}.visible                   = 'on';

config{1}.name                      = {'Hspike'};
config{1}.muse.startmarker.Hspike   = "Hspike";
config{1}.muse.endmarker.Hspike     = "Hspike";
config{1}.muse.backupdir            = fullfile(rootpath_analysis, 'markerbackup');

config{1}.hyp.imagesavedir          = fullfile(rootpath_analysis, 'images', 'hspike');
config{1}.hyp.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{1}.hyp.markerdir             = fullfile(rootpath_analysis, 'data',   'hspike');
config{1}.hyp.micromedchannel       = 'F3p6';
config{1}.hyp.markers               = [];
config{1}.hyp.overwrite             = 'append';
config{1}.hyp.spikewindow           = 60;
config{1}.hyp.spikewindowoverlap    = 0.5;

config{1}.epoch.toi.Hspike          = [-0.5  1];
config{1}.epoch.pad.Hspike          = 0.5;

config{1}.LFP.name                  = {'Hspike'};
config{1}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 250;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow.Hspike = [-1, -0.5];