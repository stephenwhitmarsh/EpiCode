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

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
    rootpath_data       = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data';
else
    error('Platform not supported')
end

%% 

% Patient 1, perivetricular heterotopia #1
disp('setting parameters');

% patient 1
muscale = 50; % multiunit scale

% Patient 1, perivetricular heterotopia #1
config{1}.os                        = os;
config{1}.name                      = {'PSW','FA','ES'};
config{1}.prefix                    = 'N1-';
config{1}.muse.startend             = { 'VF1__START__','VF1__END__';'RR__START__','RR__END__';'P','P'};   % start and end Muse marker

% config{1}.datasavedir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/data';         % where to write data
% config{1}.imagesavedir              = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/images';       % where to print images
% config{1}.rawdir                    = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg';

config{1}.patientdir                = fullfile(rootpath_data,     'pnh', 'pat_02230_0674', 'eeg');
config{1}.rawdir                    = fullfile(rootpath_data,     'pnh', 'pat_02230_0674', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images
config{1}.directory_searchstring    = '02230_2015-02-*';

config{1}.labels.micro              = {'m1pNs_1','m1pNs_2','m1pNs_4','m1pNs_6','m1pNs_7','m1pNs_8'};
config{1}.labels.macro              = {'_1pNs_1','_1pNs_2','_1pNs_3','_1pNs_4','_1pNs_5','_1pNs_6','_1pNs_7','_1pNs_8'};

config{1}.align.channel             = {'m1pNs_4','m1pNs_4','m1pNs_4'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{1}.align.method              = {'max','first','max'};                                                              % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
config{1}.align.filter              = {'bp','bp','bp'};
config{1}.align.freq                = {[1, 4],[1, 4],[1, 40]};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
config{1}.align.hilbert             = {'no','no','no'};
config{1}.align.thresh              = [0,0.25,0.25];
config{1}.align.toiactive{1}        = [-0.1,  0.4];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toiactive{2}        = [-0.1,  0.3];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toiactive{3}        = [-0.1,  0.1];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{1}.align.toibaseline{1}      = [-1.0, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toibaseline{2}      = [-1.0, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{1}.align.toibaseline{3}      = [-1.0  -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 1000;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow{1}     = [-2, -1];
config{1}.LFP.baselinewindow{2}     = [-2, -1];
config{1}.LFP.baselinewindow{3}     = [-1, -0.5];
config{1}.LFP.slidestep             = [0.01, 0.01, 0.01];

config{1}.epoch.toi{1}              = [-2.00,  2.00];                                                                  % list of onset timing with respect to start-marker (s)
config{1}.epoch.toi{2}              = [-2.00,  2.00];                                                                  % list of onset timing with respect to start-marker (s)
config{1}.epoch.toi{3}              = [-1.00,  1.00];                                                                  % list of onset timing with respect to start-marker (s)
config{1}.epoch.pad                 = [0.5, 0.5, 0.5];

config{1}.circus.channel            = {'m1pNs_1','m1pNs_2','m1pNs_6','m1pNs_7','m1pNs_8'};
config{1}.circus.reref              = 'yes';
config{1}.circus.refchan            = 'm1pNs_4';
config{1}.circus.outputdir          = 'SpykingCircus';
config{1}.circus.suffix             = '-2';

config{1}.spike.slidestep           = [0.01, 0.01, 0.001];
config{1}.spike.toispikerate{1}     = [-0.1 0.1];           % for plotting spikerate
config{1}.spike.toispikerate{2}     = [-0.1 0.1];           % for p200mslotting spikerate
config{1}.spike.toispikerate{3}     = [-0.025, 0.025];      % for plotting spikerate
config{1}.spike.resamplefs          = 1000;
config{1}.spike.pre                 = 0.001;
config{1}.spike.post                = 0.002;
config{1}.spike.baseline            = [-0.001 -0.0005];
config{1}.spike.bltoi{1}            = [-2,    -1];
config{1}.spike.bltoi{2}            = [-2,    -1];
config{1}.spike.bltoi{3}            = [-1,  -0.5];
config{1}.spike.ISIbins             = [0:0.003:0.150]; %in s

config{1}.stats.actoi{1}            = [-0.15, 0.15];
config{1}.stats.actoi{2}            = [-0.15, 0.15];
config{1}.stats.actoi{3}            = [-0.15, 0.15]; 
config{1}.stats.bltoi{1}            = [-2,    -1];
config{1}.stats.bltoi{2}            = [-2,    -1];
config{1}.stats.bltoi{3}            = [-1,  -0.5];
config{1}.stats.alpha               = 0.025;

config{1}.plot.toi{1}               = [1500.8, 1503.3, 0];
config{1}.plot.hpfilter{1}{1}       = 'no';
config{1}.plot.hpfilter{1}{2}       = 'no';
config{1}.plot.hpfilter{1}{3}       = 'no';
config{1}.plot.hpfilter{1}{4}       = 'no';
config{1}.plot.hpfilter{1}{5}       = 'yes';
config{1}.plot.hpfilter{1}{6}       = 'yes';
config{1}.plot.hpfreq{1}{1}         = [];
config{1}.plot.hpfreq{1}{2}         = [];
config{1}.plot.hpfreq{1}{3}         = [];
config{1}.plot.hpfreq{1}{4}         = [];
config{1}.plot.hpfreq{1}{5}         = 500;
config{1}.plot.hpfreq{1}{6}         = 500;
config{1}.plot.lpfilter{1}{1}       = 'yes';
config{1}.plot.lpfilter{1}{2}       = 'yes';
config{1}.plot.lpfilter{1}{3}       = 'yes';
config{1}.plot.lpfilter{1}{4}       = 'yes';
config{1}.plot.lpfilter{1}{5}       = 'yes';
config{1}.plot.lpfilter{1}{6}       = 'yes';
config{1}.plot.lpfreq{1}{1}         = 500;
config{1}.plot.lpfreq{1}{2}         = 500;
config{1}.plot.lpfreq{1}{3}         = 500;
config{1}.plot.lpfreq{1}{4}         = 500;
config{1}.plot.lpfreq{1}{5}         = 3000;
config{1}.plot.lpfreq{1}{6}         = 3000;
config{1}.plot.fname{1}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_2.ncs';
config{1}.plot.fname{1}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_1.ncs';
config{1}.plot.fname{1}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.fname{1}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_6.ncs';
config{1}.plot.fname{1}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_1.ncs';
config{1}.plot.fname{1}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_2.ncs';
config{1}.plot.refname{1}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_3.ncs';
config{1}.plot.refname{1}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_2.ncs';
config{1}.plot.refname{1}{3}        = [];
config{1}.plot.refname{1}{4}        = [];
config{1}.plot.refname{1}{5}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.refname{1}{6}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.scale{1}{1}          = [-200 200];
config{1}.plot.scale{1}{2}          = [-200 200];
config{1}.plot.scale{1}{3}          = [-200 200];
config{1}.plot.scale{1}{4}          = [-200 200];
config{1}.plot.scale{1}{5}          = [-muscale muscale];
config{1}.plot.scale{1}{6}          = [-muscale muscale];

config{1}.plot.toi{2}               = [1976.2-0.5, 1976.2+2, -0.5];
config{1}.plot.hpfilter{2}{1}       = 'no';
config{1}.plot.hpfilter{2}{2}       = 'no';
config{1}.plot.hpfilter{2}{3}       = 'no';
config{1}.plot.hpfilter{2}{4}       = 'no';
config{1}.plot.hpfilter{2}{5}       = 'yes';
config{1}.plot.hpfilter{2}{6}       = 'yes';
config{1}.plot.hpfreq{2}{1}         = [];
config{1}.plot.hpfreq{2}{2}         = [];
config{1}.plot.hpfreq{2}{3}         = [];
config{1}.plot.hpfreq{2}{4}         = [];
config{1}.plot.hpfreq{2}{5}         = 500;
config{1}.plot.hpfreq{2}{6}         = 500;
config{1}.plot.lpfilter{2}{1}       = 'yes';
config{1}.plot.lpfilter{2}{2}       = 'yes';
config{1}.plot.lpfilter{2}{3}       = 'yes';
config{1}.plot.lpfilter{2}{4}       = 'yes';
config{1}.plot.lpfilter{2}{5}       = 'yes';
config{1}.plot.lpfilter{2}{6}       = 'yes';
config{1}.plot.lpfreq{2}{1}         = 500;
config{1}.plot.lpfreq{2}{2}         = 500;
config{1}.plot.lpfreq{2}{3}         = 500;
config{1}.plot.lpfreq{2}{4}         = 500;
config{1}.plot.lpfreq{2}{5}         = 3000;
config{1}.plot.lpfreq{2}{6}         = 3000;
config{1}.plot.fname{2}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_2.ncs';
config{1}.plot.fname{2}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_1.ncs';
config{1}.plot.fname{2}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.fname{2}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_6.ncs';
config{1}.plot.fname{2}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_1.ncs';
config{1}.plot.fname{2}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_8.ncs';
config{1}.plot.refname{2}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_3.ncs';
config{1}.plot.refname{2}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_1pNs_2.ncs';
config{1}.plot.refname{2}{3}        = [];
config{1}.plot.refname{2}{4}        = [];
config{1}.plot.refname{2}{5}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.refname{2}{6}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_17-16/02230_2015-02-25_17-16_m1pNs_4.ncs';
config{1}.plot.scale{2}{1}          = [-200 200];
config{1}.plot.scale{2}{2}          = [-200 200];
config{1}.plot.scale{2}{3}          = [-200 200];
config{1}.plot.scale{2}{4}          = [-200 200];
config{1}.plot.scale{2}{5}          = [-muscale muscale];
config{1}.plot.scale{2}{6}          = [-muscale muscale];

config{1}.plot.toi{3}               = [894.659-0.4, 894.659+0.4, -0.4];
config{1}.plot.hpfilter{3}{1}       = 'no';
config{1}.plot.hpfilter{3}{2}       = 'no';
config{1}.plot.hpfilter{3}{3}       = 'no';
config{1}.plot.hpfilter{3}{4}       = 'no';
config{1}.plot.hpfilter{3}{5}       = 'yes';
config{1}.plot.hpfilter{3}{6}       = 'yes';
config{1}.plot.hpfreq{3}{1}         = [];
config{1}.plot.hpfreq{3}{2}         = [];
config{1}.plot.hpfreq{3}{3}         = [];
config{1}.plot.hpfreq{3}{4}         = [];
config{1}.plot.hpfreq{3}{5}         = 500;
config{1}.plot.hpfreq{3}{6}         = 500;
config{1}.plot.lpfilter{3}{1}       = 'yes';
config{1}.plot.lpfilter{3}{2}       = 'yes';
config{1}.plot.lpfilter{3}{3}       = 'yes';
config{1}.plot.lpfilter{3}{4}       = 'yes';
config{1}.plot.lpfilter{3}{5}       = 'yes';
config{1}.plot.lpfilter{3}{6}       = 'yes';
config{1}.plot.lpfreq{3}{1}         = 500;
config{1}.plot.lpfreq{3}{2}         = 500;
config{1}.plot.lpfreq{3}{3}         = 500;
config{1}.plot.lpfreq{3}{4}         = 500;
config{1}.plot.lpfreq{3}{5}         = 3000;
config{1}.plot.lpfreq{3}{6}         = 3000;
config{1}.plot.fname{3}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_1pNs_2.ncs';
config{1}.plot.fname{3}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_1pNs_1.ncs';
config{1}.plot.fname{3}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_4.ncs';
config{1}.plot.fname{3}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_6.ncs';
config{1}.plot.fname{3}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_1.ncs';
config{1}.plot.fname{3}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_8.ncs';

config{1}.plot.refname{3}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_1pNs_3.ncs';
config{1}.plot.refname{3}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_1pNs_2.ncs';
config{1}.plot.refname{3}{3}        = [];
config{1}.plot.refname{3}{4}        = [];
config{1}.plot.refname{3}{5}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_4.ncs';
config{1}.plot.refname{3}{6}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg/02230_2015-02-25_14-36/02230_2015-02-25_14-36_m1pNs_4.ncs';
config{1}.plot.scale{3}{1}          = [-200 200];
config{1}.plot.scale{3}{2}          = [-200 200];
config{1}.plot.scale{3}{3}          = [-200 200];
config{1}.plot.scale{3}{4}          = [-200 200];
config{1}.plot.scale{3}{5}          = [-muscale muscale];
config{1}.plot.scale{3}{6}          = [-muscale muscale];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patient 2, first nodule %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
config{2}.os                        = os;
config{2}.name                      = { 'PSW','FA','ES'};
config{2}.prefix                    = 'N2-';

% config{2}.datasavedir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';        % where to write data
% config{2}.imagesavedir              = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';      % where to print images
% config{2}.rawdir                    = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg';

config{2}.patientdir                = fullfile(rootpath_data,     'pnh', 'pat_02614_1073', 'eeg');
config{2}.rawdir                    = fullfile(rootpath_data,     'pnh', 'pat_02614_1073', 'eeg');
config{2}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{2}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images

config{2}.directory_searchstring    = '02614_2018-06-*';
config{2}.labels.micro              = 
config{2}.labels.macro              = {'_Casd_1','_Casd_2','_Casd_3','_Casd_4','_Casd_5','_Casd_6','_Casd_7','_Casd_8'};

config{2}.muse.startend             = {'H02614_01__START__','H02614_01__END__';'H02614_04__START__','H02614_04__END__';'H02614_02','H02614_02'};   % start and end Muse marker

config{2}.align.channel             = {'mCasd_2','mCasd_2','mCasd_2'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{2}.align.method              = {'first','first','max'};                                 % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
config{2}.align.filter              = {'lp','bp','bp'};                                        % hipass part of bandpass seems to create a huge phase shift!!!

config{2}.align.freq                = {4, [5, 20], [1 ,40]};                                    % lowpass filter freq to smooth peak detection (Hz)
config{2}.align.hilbert             = {'no','no','no'};

config{2}.align.thresh              = [0.25,1,0.25];
config{2}.align.toiactive{1}        = [-0.1,  0.4];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toiactive{2}        = [-0.1,  0.5];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toiactive{3}        = [-0.1,  0.2];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{2}.align.toibaseline{1}      = [-1.0, -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toibaseline{2}      = [-1.0  -0.1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{2}.align.toibaseline{3}      = [-1.0  -0.2];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

config{2}.epoch.toi{1}              = [-2.00, 2.00];                                           % list of onset timing with respect to start-marker (s)
config{2}.epoch.toi{2}              = [-2.00, 2.50];                                           % list of onset timing with respect to start-marker (s)
config{2}.epoch.toi{3}              = [-2.00, 2.00];
config{2}.epoch.pad                 = [ 0.5, 0.5, 0.5];

config{2}.LFP.hpfilter              = 'no';
config{2}.LFP.hpfreq                = 1;
config{2}.LFP.resamplefs            = 1000;
config{2}.LFP.baseline              = 'yes';
config{2}.LFP.baselinewindow{1}     = [-2, -1];
config{2}.LFP.baselinewindow{2}     = [-2, -1];
config{2}.LFP.baselinewindow{3}     = [-2, -1];
config{2}.LFP.slidestep             = [0.01, 0.01, 0.01];

config{2}.circus.channel            = {'mCasd_1','mCasd_2','mCasd_3','mCasd_4','mCasd_5','mCasd_6','mCasd_7'};
config{2}.circus.reref              = 'no';
config{2}.circus.refchan            = '';
config{2}.circus.outputdir          = 'SpykingCircus';
config{2}.circus.suffix             = '-2';

config{2}.spike.slidestep           = [0.01,0.01,0.001];
config{2}.spike.toispikerate{1}     = [-0.1,  0.1];         % for plotting spikerate
config{2}.spike.toispikerate{2}     = [-0.1,  0.1];         % for plotting spikerate
config{2}.spike.toispikerate{3}     = [-0.025,  0.025];       % for plotting spikerate
config{2}.spike.resamplefs          = 1000;
config{2}.spike.pre                 = 0.001;
config{2}.spike.post                = 0.002;
config{2}.spike.baseline            = [-0.001 -0.0005];
config{2}.spike.bltoi{1}            = [-2, -1];
config{2}.spike.bltoi{2}            = [-2, -1];
config{2}.spike.bltoi{3}            = [-2, -1];

config{2}.stats.actoi{1}            = [-0.15, 0.15];
config{2}.stats.actoi{2}            = [-0.15, 0.15];
config{2}.stats.actoi{3}            = [-0.15, 0.15]; 
config{2}.stats.bltoi{1}            = [-2, -1];  % first is baseline
config{2}.stats.bltoi{2}            = [-2, -1];  % first is baseline
config{2}.stats.bltoi{3}            = [-2, -1];  % first is baseline
config{2}.stats.alpha               = 0.025;  % first is baseline

config{2}.plot.toi{1}               = [1004.0, 1008.0, 0];
config{2}.plot.hpfilter{1}{1}       = 'no';
config{2}.plot.hpfilter{1}{2}       = 'no';
config{2}.plot.hpfilter{1}{3}       = 'no';
config{2}.plot.hpfilter{1}{4}       = 'no';
config{2}.plot.hpfilter{1}{5}       = 'yes';
config{2}.plot.hpfilter{1}{6}       = 'yes';
config{2}.plot.hpfreq{1}{1}         = [];
config{2}.plot.hpfreq{1}{2}         = [];
config{2}.plot.hpfreq{1}{3}         = [];
config{2}.plot.hpfreq{1}{4}         = [];
config{2}.plot.hpfreq{1}{5}         = 500;
config{2}.plot.hpfreq{1}{6}         = 500;
config{2}.plot.lpfilter{1}{1}       = 'yes';
config{2}.plot.lpfilter{1}{2}       = 'yes';
config{2}.plot.lpfilter{1}{3}       = 'yes';
config{2}.plot.lpfilter{1}{4}       = 'yes';
config{2}.plot.lpfilter{1}{5}       = 'yes';
config{2}.plot.lpfilter{1}{6}       = 'yes';
config{2}.plot.lpfreq{1}{1}         = 500;
config{2}.plot.lpfreq{1}{2}         = 500;
config{2}.plot.lpfreq{1}{3}         = 500;
config{2}.plot.lpfreq{1}{4}         = 500;
config{2}.plot.lpfreq{1}{5}         = 3000;
config{2}.plot.lpfreq{1}{6}         = 3000;

config{2}.plot.fname{1}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_2.ncs';
config{2}.plot.fname{1}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_1.ncs';
config{2}.plot.fname{1}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_2.ncs';
config{2}.plot.fname{1}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_5.ncs';
config{2}.plot.fname{1}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_2.ncs';
config{2}.plot.fname{1}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_6.ncs';
config{2}.plot.refname{1}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_3.ncs';
config{2}.plot.refname{1}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_2.ncs';
config{2}.plot.refname{1}{3}        = [];
config{2}.plot.refname{1}{4}        = [];
config{2}.plot.refname{1}{5}        = [];
config{2}.plot.refname{1}{6}        = [];
config{2}.plot.scale{1}{1}          = [-200 200];
config{2}.plot.scale{1}{2}          = [-200 200];
config{2}.plot.scale{1}{3}          = [-200 200];
config{2}.plot.scale{1}{4}          = [-200 200];
config{2}.plot.scale{1}{5}          = [-muscale muscale];
config{2}.plot.scale{1}{6}          = [-muscale muscale];

config{2}.plot.toi{2}               = [412.9, 415.4, 0];
config{2}.plot.hpfilter{2}{1}       = 'no';
config{2}.plot.hpfilter{2}{2}       = 'no';
config{2}.plot.hpfilter{2}{3}       = 'no';
config{2}.plot.hpfilter{2}{4}       = 'no';
config{2}.plot.hpfilter{2}{5}       = 'yes';
config{2}.plot.hpfilter{2}{6}       = 'yes';
config{2}.plot.hpfreq{2}{1}         = [];
config{2}.plot.hpfreq{2}{2}         = [];
config{2}.plot.hpfreq{2}{3}         = [];
config{2}.plot.hpfreq{2}{4}         = [];
config{2}.plot.hpfreq{2}{5}         = 500;
config{2}.plot.hpfreq{2}{6}         = 500;
config{2}.plot.lpfilter{2}{1}       = 'yes';
config{2}.plot.lpfilter{2}{2}       = 'yes';
config{2}.plot.lpfilter{2}{3}       = 'yes';
config{2}.plot.lpfilter{2}{4}       = 'yes';
config{2}.plot.lpfilter{2}{5}       = 'yes';
config{2}.plot.lpfilter{2}{6}       = 'yes';
config{2}.plot.lpfreq{2}{1}         = 500;
config{2}.plot.lpfreq{2}{2}         = 500;
config{2}.plot.lpfreq{2}{3}         = 500;
config{2}.plot.lpfreq{2}{4}         = 500;
config{2}.plot.lpfreq{2}{5}         = 3000;
config{2}.plot.lpfreq{2}{6}         = 3000;
config{2}.plot.fname{2}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_2.ncs';
config{2}.plot.fname{2}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_1.ncs';
config{2}.plot.fname{2}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_2.ncs';
config{2}.plot.fname{2}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_5.ncs';
config{2}.plot.fname{2}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_2.ncs';
config{2}.plot.fname{2}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_mCasd_4.ncs';
config{2}.plot.refname{2}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_3.ncs';
config{2}.plot.refname{2}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_17-23/02614_2018-06-12_17-23_Casd_2.ncs';
config{2}.plot.refname{2}{3}        = [];
config{2}.plot.refname{2}{4}        = [];
config{2}.plot.refname{2}{5}        = [];
config{2}.plot.refname{2}{6}        = [];
config{2}.plot.scale{2}{1}          = [-200 200];
config{2}.plot.scale{2}{2}          = [-200 200];
config{2}.plot.scale{2}{3}          = [-200 200];
config{2}.plot.scale{2}{4}          = [-200 200];
config{2}.plot.scale{2}{5}          = [-muscale muscale];
config{2}.plot.scale{2}{6}          = [-muscale muscale];
% 
% config{2}.plot.toi{3}               = [8.06, 8.36, 0];
% config{2}.plot.toi{3}               = [5865.233-1, 894.659+1, -1];
config{2}.plot.toi{3}               = [1664.850-0.20, 1664.850+0.4, -0.20];

config{2}.plot.hpfilter{3}{1}       = 'no';
config{2}.plot.hpfilter{3}{2}       = 'no';
config{2}.plot.hpfilter{3}{3}       = 'no';
config{2}.plot.hpfilter{3}{4}       = 'no';
config{2}.plot.hpfilter{3}{5}       = 'yes';
config{2}.plot.hpfilter{3}{6}       = 'yes';
config{2}.plot.hpfreq{3}{1}         = [];
config{2}.plot.hpfreq{3}{2}         = [];
config{2}.plot.hpfreq{3}{3}         = [];
config{2}.plot.hpfreq{3}{4}         = [];
config{2}.plot.hpfreq{3}{5}         = 500;
config{2}.plot.hpfreq{3}{6}         = 500;
config{2}.plot.lpfilter{3}{1}       = 'yes';
config{2}.plot.lpfilter{3}{2}       = 'yes';
config{2}.plot.lpfilter{3}{3}       = 'yes';
config{2}.plot.lpfilter{3}{4}       = 'yes';
config{2}.plot.lpfilter{3}{5}       = 'yes';
config{2}.plot.lpfilter{3}{6}       = 'yes';
config{2}.plot.lpfreq{3}{1}         = 500;
config{2}.plot.lpfreq{3}{2}         = 500;
config{2}.plot.lpfreq{3}{3}         = 500;
config{2}.plot.lpfreq{3}{4}         = 500;
config{2}.plot.lpfreq{3}{5}         = 3000;
config{2}.plot.lpfreq{3}{6}         = 3000;
config{2}.plot.fname{3}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_Casd_2.ncs';
config{2}.plot.fname{3}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_Casd_1.ncs';
config{2}.plot.fname{3}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mCasd_2.ncs';
config{2}.plot.fname{3}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mCasd_5.ncs';
config{2}.plot.fname{3}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mCasd_2.ncs';
config{2}.plot.fname{3}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mCasd_5.ncs';
config{2}.plot.refname{3}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_Casd_3.ncs';
config{2}.plot.refname{3}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_Casd_2.ncs';
config{2}.plot.refname{3}{3}        = [];
config{2}.plot.refname{3}{4}        = [];
config{2}.plot.refname{3}{5}        = [];
config{2}.plot.refname{3}{6}        = [];
config{2}.plot.scale{3}{1}          = [-200 200];
config{2}.plot.scale{3}{2}          = [-200 200];
config{2}.plot.scale{3}{3}          = [-200 200];
config{2}.plot.scale{3}{4}          = [-200 200];
config{2}.plot.scale{3}{5}          = [-muscale muscale];
config{2}.plot.scale{3}{6}          = [-muscale muscale];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patient 2, second nodule  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config{3}.os                        = os;
config{3}.name                      = {'FA','ES'};
config{3}.prefix                    = 'N3-';

% config{3}.datasavedir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';
% config{3}.imagesavedir              = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';
% config{3}.rawdir                    = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg';

config{3}.patientdir                = fullfile(rootpath_data,     'pnh', 'pat_02614_1073', 'eeg');
config{3}.rawdir                    = fullfile(rootpath_data,     'pnh', 'pat_02614_1073', 'eeg');
config{3}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{3}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images

config{3}.directory_searchstring    = '02614_2018-06-*';
config{3}.labels.micro              = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_5','mTNmi_6','mTNmi_7','mTNmi_8'};
config{3}.labels.macro              = {'_TNmi_1','_TNmi_2','_TNmi_3','_TNmi_4','_TNmi_5','_TNmi_6','_TNmi_7','_TNmi_8'};

config{3}.muse.startend             = {'H02614_07__START__','H02614_07__END__';'H02614_09','H02614_09'};  % start and end Muse marker

config{3}.align.channel             = {'mTNmi_3','mTNmi_3'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
config{3}.align.method              = {'max','min'};                                                            % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
config{3}.align.filter              = {'bp','lp'};
config{3}.align.freq                = {[5, 50], 50};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
config{3}.align.hilbert             = {'yes','no'};
config{3}.align.thresh              = [ 0.25,  -inf];
config{3}.align.toiactive{1}        = [-0.10,  0.40];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{3}.align.toiactive{2}        = [-0.025, 0.05];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
config{3}.align.toibaseline{1}      = [-1.00, -0.10];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
config{3}.align.toibaseline{2}      = [-0.20  -0.10];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];

config{3}.epoch.toi{1}              = [-2.00,  2.50];                                                                  % list of onset timing with respect to start-marker (s)
config{3}.epoch.toi{2}              = [-1.00,  1.00];                                                                  % list of onset timing with respect to start-marker (s)
config{3}.epoch.pad                 = [ 0.5,  0.5];

config{3}.LFP.hpfilter              = 'no';
config{3}.LFP.hpfreq                = 1;
config{3}.LFP.resamplefs            = 1000;
config{3}.LFP.baseline              = 'yes';
config{3}.LFP.baselinewindow{1}     = [-2,-1];
config{3}.LFP.baselinewindow{2}     = [-1, -0.5];
config{3}.LFP.slidestep             = [0.01, 0.01];

config{3}.circus.channel            = {'mTNmi_1','mTNmi_2','mTNmi_3','mTNmi_4','mTNmi_5','mTNmi_6','mTNmi_7','mTNmi_8'};
config{3}.circus.reref              = 'no';
config{3}.circus.refchan            = '';
config{3}.circus.outputdir          = 'SpykingCircus';
config{3}.circus.suffix             = '-2';

config{3}.spike.slidestep           = [0.01,   0.001];
config{3}.spike.toispikerate{1}     = [-0.1,   0.1];        % for plotting spikerate
config{3}.spike.toispikerate{2}     = [-0.025,  0.025];        % for plotting spikerate
config{3}.spike.resamplefs          = 1000;
config{3}.spike.pre                 = 0.001;
config{3}.spike.post                = 0.002;
config{3}.spike.baseline            = [-0.001 -0.0005];
config{3}.spike.bltoi{1}            = [-2,    -1];
config{3}.spike.bltoi{2}            = [-1,  -0.5];

config{3}.stats.actoi{1}            = [-0.15, 0.15];
config{3}.stats.actoi{2}            = [-0.15, 0.15]; 
config{3}.stats.bltoi{1}            = [-2, -1];
config{3}.stats.bltoi{2}            = [-1, -0.5];
config{3}.stats.alpha               = 0.025;

config{3}.plot.toi{1}               = [1181, 1185, 0];
config{3}.plot.hpfilter{1}{1}       = 'no';
config{3}.plot.hpfilter{1}{2}       = 'no';
config{3}.plot.hpfilter{1}{3}       = 'no';
config{3}.plot.hpfilter{1}{4}       = 'no';
config{3}.plot.hpfilter{1}{5}       = 'yes';
config{3}.plot.hpfilter{1}{6}       = 'yes';
config{3}.plot.hpfreq{1}{1}         = [];
config{3}.plot.hpfreq{1}{2}         = [];
config{3}.plot.hpfreq{1}{3}         = [];
config{3}.plot.hpfreq{1}{4}         = [];
config{3}.plot.hpfreq{1}{5}         = 500;
config{3}.plot.hpfreq{1}{6}         = 500;
config{3}.plot.lpfilter{1}{1}       = 'yes';
config{3}.plot.lpfilter{1}{2}       = 'yes';
config{3}.plot.lpfilter{1}{3}       = 'yes';
config{3}.plot.lpfilter{1}{4}       = 'yes';
config{3}.plot.lpfilter{1}{5}       = 'yes';
config{3}.plot.lpfilter{1}{6}       = 'yes';
config{3}.plot.lpfreq{1}{1}         = 500;
config{3}.plot.lpfreq{1}{2}         = 500;
config{3}.plot.lpfreq{1}{3}         = 500;
config{3}.plot.lpfreq{1}{4}         = 500;
config{3}.plot.lpfreq{1}{5}         = 3000;
config{3}.plot.lpfreq{1}{6}         = 3000;
config{3}.plot.fname{1}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_1.ncs';
config{3}.plot.fname{1}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_2.ncs';
config{3}.plot.fname{1}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_3.ncs';
config{3}.plot.fname{1}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_5.ncs';
config{3}.plot.fname{1}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_3.ncs';
config{3}.plot.fname{1}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_7.ncs';
config{3}.plot.refname{1}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_2.ncs';
config{3}.plot.refname{1}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_3.ncs';
config{3}.plot.refname{1}{3}        = [];
config{3}.plot.refname{1}{4}        = [];
config{3}.plot.refname{1}{5}        = [];
config{3}.plot.refname{1}{6}        = [];
config{3}.plot.scale{1}{1}          = [-200 200];
config{3}.plot.scale{1}{2}          = [-200 200];
config{3}.plot.scale{1}{3}          = [-200 200];
config{3}.plot.scale{1}{4}          = [-200 200];
config{3}.plot.scale{1}{5}          = [-muscale muscale];
config{3}.plot.scale{1}{6}          = [-muscale muscale];

config{3}.plot.toi{2}               = [722.5, 723.3, -0.4];
config{3}.plot.hpfilter{2}{1}       = 'no';
config{3}.plot.hpfilter{2}{2}       = 'no';
config{3}.plot.hpfilter{2}{3}       = 'no';
config{3}.plot.hpfilter{2}{4}       = 'no';
config{3}.plot.hpfilter{2}{5}       = 'yes';
config{3}.plot.hpfilter{2}{6}       = 'yes';
config{3}.plot.hpfreq{2}{1}         = [];
config{3}.plot.hpfreq{2}{2}         = [];
config{3}.plot.hpfreq{2}{3}         = [];
config{3}.plot.hpfreq{2}{4}         = [];
config{3}.plot.hpfreq{2}{5}         = 500;
config{3}.plot.hpfreq{2}{6}         = 500;
config{3}.plot.lpfilter{2}{1}       = 'yes';
config{3}.plot.lpfilter{2}{2}       = 'yes';
config{3}.plot.lpfilter{2}{3}       = 'yes';
config{3}.plot.lpfilter{2}{4}       = 'yes';
config{3}.plot.lpfilter{2}{5}       = 'yes';
config{3}.plot.lpfilter{2}{6}       = 'yes';
config{3}.plot.lpfreq{2}{1}         = 500;
config{3}.plot.lpfreq{2}{2}         = 500;
config{3}.plot.lpfreq{2}{3}         = 500;
config{3}.plot.lpfreq{2}{4}         = 500;
config{3}.plot.lpfreq{2}{5}         = 3000;
config{3}.plot.lpfreq{2}{6}         = 3000;
config{3}.plot.fname{2}{1}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_1.ncs';
config{3}.plot.fname{2}{2}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_2.ncs';
config{3}.plot.fname{2}{3}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_3.ncs';
config{3}.plot.fname{2}{4}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_4.ncs';
config{3}.plot.fname{2}{5}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_3.ncs';
config{3}.plot.fname{2}{6}          = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_mTNmi_5.ncs';
config{3}.plot.refname{2}{1}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_2.ncs';
config{3}.plot.refname{2}{2}        = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg/02614_2018-06-12_15-23/02614_2018-06-12_15-23_TNmi_3.ncs';
config{3}.plot.refname{2}{3}        = [];
config{3}.plot.refname{2}{4}        = [];
config{3}.plot.refname{2}{5}        = [];
config{3}.plot.refname{2}{6}        = [];
config{3}.plot.scale{2}{1}          = [-200 200];
config{3}.plot.scale{2}{2}          = [-200 200];
config{3}.plot.scale{2}{3}          = [-200 200];
config{3}.plot.scale{2}{4}          = [-200 200];
config{3}.plot.scale{2}{5}          = [-muscale muscale];
config{3}.plot.scale{2}{6}          = [-muscale muscale];


% plot seizures:

config{4}.os                        = os;
config{4}.name                      = {'Seizure_N1'};
config{4}.prefix                    = 'Seizures-N1-';
% config{4}.imagesavedir              = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/images';
% config{4}.datasavedir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/data';         % where to write data
config{4}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{4}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images

config{4}.plot.toi{1}               = [6303.04, 6317.04, -1];
config{4}.plot.hpfilter{1}{1}       = 'no';
config{4}.plot.hpfilter{1}{2}       = 'no';
config{4}.plot.hpfilter{1}{3}       = 'no';
config{4}.plot.hpfilter{1}{4}       = 'no';
config{4}.plot.hpfilter{1}{5}       = 'no';
config{4}.plot.hpfilter{1}{6}       = 'no';
config{4}.plot.hpfilter{1}{7}       = 'yes';
config{4}.plot.hpfilter{1}{8}      = 'yes';
config{4}.plot.hpfreq{1}{1}         = [];
config{4}.plot.hpfreq{1}{2}         = [];
config{4}.plot.hpfreq{1}{3}         = [];
config{4}.plot.hpfreq{1}{4}         = [];
config{4}.plot.hpfreq{1}{5}         = [];
config{4}.plot.hpfreq{1}{6}         = [];
config{4}.plot.hpfreq{1}{7}         = 500;
config{4}.plot.hpfreq{1}{8}         = 500;

config{4}.plot.lpfilter{1}{1}       = 'yes';
config{4}.plot.lpfilter{1}{2}       = 'yes';
config{4}.plot.lpfilter{1}{3}       = 'yes';
config{4}.plot.lpfilter{1}{4}       = 'yes';
config{4}.plot.lpfilter{1}{5}       = 'yes';
config{4}.plot.lpfilter{1}{6}       = 'yes';
config{4}.plot.lpfilter{1}{7}       = 'yes';
config{4}.plot.lpfilter{1}{8}       = 'yes';
config{4}.plot.lpfreq{1}{1}         = 500;
config{4}.plot.lpfreq{1}{2}         = 500;
config{4}.plot.lpfreq{1}{3}         = 500;
config{4}.plot.lpfreq{1}{4}         = 500;
config{4}.plot.lpfreq{1}{5}         = 500;
config{4}.plot.lpfreq{1}{6}         = 500;
config{4}.plot.lpfreq{1}{7}         = 3000;
config{4}.plot.lpfreq{1}{8}         = 3000;

config{4}.plot.fname{1}{1}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_4.ncs';
config{4}.plot.fname{1}{2}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_3.ncs';
config{4}.plot.fname{1}{3}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_2.ncs';
config{4}.plot.fname{1}{4}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_1.ncs';
config{4}.plot.fname{1}{5}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_m1pNs_1.ncs';
config{4}.plot.fname{1}{6}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_m1pNs_6.ncs';
config{4}.plot.fname{1}{7}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_m1pNs_1.ncs';
config{4}.plot.fname{1}{8}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_m1pNs_6.ncs';
config{4}.plot.refname{1}{1}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_5.ncs';
config{4}.plot.refname{1}{2}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_4.ncs';
config{4}.plot.refname{1}{3}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_3.ncs';
config{4}.plot.refname{1}{4}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02230_0674/eeg/02230_2015-02-25_15-16/02230_2015-02-25_15-16_1pNs_2.ncs';
config{4}.plot.refname{1}{5}        = [];
config{4}.plot.refname{1}{6}        = [];
config{4}.plot.refname{1}{7}        = [];
config{4}.plot.refname{1}{8}        = [];

config{4}.plot.scale{1}{1}          = [-250 250];
config{4}.plot.scale{1}{2}          = [-250 250];
config{4}.plot.scale{1}{3}          = [-250 250];
config{4}.plot.scale{1}{4}          = [-250 250];
config{4}.plot.scale{1}{5}          = [-250 250];
config{4}.plot.scale{1}{6}          = [-250 250];
config{4}.plot.scale{1}{7}          = [-muscale muscale];
config{4}.plot.scale{1}{8}         = [-muscale muscale];

%%%%
config{5}.os                        = os;
config{5}.name                      = {'Seizure_N3'};
config{5}.prefix                    = 'Seizures-N3-';
% config{5}.imagesavedir              = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/images';
% config{5}.datasavedir               = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';        % where to write data
config{5}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{5}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images

config{5}.plot.toi{1}               = [2054.65-2, 2068.65-2, -1];
config{5}.plot.hpfilter{1}{1}       = 'no';
config{5}.plot.hpfilter{1}{2}       = 'no';
config{5}.plot.hpfilter{1}{3}       = 'no';
config{5}.plot.hpfilter{1}{4}       = 'no';
config{5}.plot.hpfilter{1}{5}       = 'no';
config{5}.plot.hpfilter{1}{6}       = 'no';
config{5}.plot.hpfilter{1}{7}       = 'yes';
config{5}.plot.hpfilter{1}{8}       = 'yes';
config{5}.plot.hpfreq{1}{1}         = [];
config{5}.plot.hpfreq{1}{2}         = [];
config{5}.plot.hpfreq{1}{3}         = [];
config{5}.plot.hpfreq{1}{4}         = [];
config{5}.plot.hpfreq{1}{5}         = [];
config{5}.plot.hpfreq{1}{6}         = [];
config{5}.plot.hpfreq{1}{7}        = 500;
config{5}.plot.hpfreq{1}{8}        = 500;


config{5}.plot.lpfilter{1}{1}       = 'yes';
config{5}.plot.lpfilter{1}{2}       = 'yes';
config{5}.plot.lpfilter{1}{3}       = 'yes';
config{5}.plot.lpfilter{1}{4}       = 'yes';
config{5}.plot.lpfilter{1}{5}       = 'yes';
config{5}.plot.lpfilter{1}{6}       = 'yes';
config{5}.plot.lpfilter{1}{7}       = 'yes';
config{5}.plot.lpfilter{1}{8}       = 'yes';
config{5}.plot.lpfreq{1}{1}         = 500;
config{5}.plot.lpfreq{1}{2}         = 500;
config{5}.plot.lpfreq{1}{3}         = 500;
config{5}.plot.lpfreq{1}{4}         = 500;
config{5}.plot.lpfreq{1}{5}         = 500;
config{5}.plot.lpfreq{1}{6}         = 500;
config{5}.plot.lpfreq{1}{7}        = 3000;
config{5}.plot.lpfreq{1}{8}        = 3000;

config{5}.plot.fname{1}{1}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_4.ncs';
config{5}.plot.fname{1}{2}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_3.ncs';
config{5}.plot.fname{1}{3}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_2.ncs';
config{5}.plot.fname{1}{4}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_1.ncs';
config{5}.plot.fname{1}{5}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_mTNmi_5.ncs';
config{5}.plot.fname{1}{6}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_mTNmi_7.ncs';
config{5}.plot.fname{1}{7}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_mTNmi_5.ncs';
config{5}.plot.fname{1}{8}          = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_mTNmi_7.ncs';
config{5}.plot.refname{1}{1}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_5.ncs';
config{5}.plot.refname{1}{2}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_4.ncs';
config{5}.plot.refname{1}{3}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_3.ncs';
config{5}.plot.refname{1}{4}        = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/02614_2018-06-25_10-47/02614_2018-06-25_10-47_TNmi_2.ncs';
config{5}.plot.refname{1}{5}        = [];
config{5}.plot.refname{1}{6}        = [];
config{5}.plot.refname{1}{7}        = [];
config{5}.plot.refname{1}{8}        = [];

config{5}.plot.scale{1}{1}          = [-350 350];
config{5}.plot.scale{1}{2}          = [-350 350];
config{5}.plot.scale{1}{3}          = [-350 350];
config{5}.plot.scale{1}{4}          = [-350 350];
config{5}.plot.scale{1}{5}          = [-350 350];
config{5}.plot.scale{1}{6}          = [-350 350];
config{5}.plot.scale{1}{7}          = [-muscale muscale];
config{5}.plot.scale{1}{8}          = [-muscale muscale];
end
