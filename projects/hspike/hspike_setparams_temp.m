function [config] = hspike_setparams_temp

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
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.epimicro\patients\raw';
    os                  = 'windows';
else
    error('Platform not supported')
end

%% Patient 1
config{1}.name                      = {'Hspike','SpikeDetect'};
config{1}.prefix                    = '2711-'; % edit in code
config{1}.muse.startend             = {'Hspike','Hspike'; 'SpikeDetect','SpikeDetect'};   % start and end Muse marker
config{1}.muse.backupdir             = fullfile(rootpath_analysis, 'markerbackup');
config{1}.muse.write                = false;
config{1}.rawdir                    = fullfile(rootpath_data,       'pat_02711_1193', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis,   'data',   'hspike');         % where to write data
config{1}.imagesavedir              = fullfile(rootpath_analysis,   'images', 'hspike');       % where to print images
config{1}.directorylist             = [];
config{1}.directorylist{1}          =  {'02711_2019-04-17_12-29',...
    '02711_2019-04-17_14-29',...
    '02711_2019-04-17_16-29',...
    '02711_2019-04-17_18-29'};

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

config{1}.LFP.name                  = {'Hspike'};
config{1}.LFP.channel               = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 250;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow{1}     = [-2, -1];
config{1}.LFP.baselinewindow{2}     = [-2, -1];
config{1}.LFP.slidestep             = [0.01];
config{1}.LFP.write                 = false;

config{1}.epoch.toi{1}              = [-0.5  1];
config{1}.epoch.toi{2}              = [-0.5  1];
config{1}.epoch.pad                 = {0.5, 0.5, 0.5};

config{1}.align.channel             = {'_HaT2_1','_HaT2_2','_HaT2_3','_HaT2_4','_HaT2_5'};
config{1}.align.reref               = 'yes';
config{1}.align.refmethod           = 'bipolar';

