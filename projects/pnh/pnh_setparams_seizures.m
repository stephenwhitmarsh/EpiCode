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
    
% if ismac
%     error('Platform not supported')
% elseif isunix
%     rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
%     rootpath_data       = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data';
% elseif ispc
%     rootpath_analysis	= '\\lexport\iss01.charpier\analyses\stephen.whitmarsh';
%     rootpath_data       = '\\lexport\iss01.charpier\analyses\stephen.whitmarsh\data';
% else
%     error('Platform not supported')
% end

%% 

% Patient 1, perivetricular heterotopia #1
disp('setting parameters');

% patient 1
muscale = 50; % multiunit scale

% Patient 1, perivetricular heterotopia #1
config{1}.prefix                    = '2230-N1-seizures-';
config{1}.rawdir                    = fullfile(rootpath_data,     'pat_02230_0674', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data', 'pnh');         % where to write data
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'pnh');       % where to print images

config{1}.name                      = {'seizure'};
config{1}.muse.startmarker.seizure  = "CriseStart";
config{1}.muse.endmarker.seizure    = "CriseEnd";
config{1}.muse.backupdir            = fullfile(rootpath_analysis, 'markerbackup');

config{1}.LFP.name                  = 'seizure';
config{1}.LFP.hpfilter              = 'no';
config{1}.LFP.hpfreq                = 1;
config{1}.LFP.resamplefs            = 1000;
config{1}.LFP.baseline              = 'yes';
config{1}.LFP.baselinewindow.seizure = [-1, -0.5];
config{1}.LFP.slidestep             = [0.01, 0.01, 0.01];
config{1}.LFP.channel               = {'_2pNi_1','_1pNs_1','_1pHe_1'};


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

config{1}.epoch.toi.seizure         = [-1  10];
config{1}.epoch.pad.seizure         = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patient 2, first nodule %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

config{2}                           = config{1};
config{2}.prefix                    = '2614-N2-seizures-';
config{2}.rawdir                    = fullfile(rootpath_data,     'pat_02614_1073', 'eeg');
config{2}.LFP.channel               = {'_TNmg_1','_TNmi_1'};

%%% DATA %%%
config{1}.directorylist             = [];
config{1}.directorylist{1}          = { '02230_2015-02-25_15-16',...
                                        '02230_2015-02-25_17-16',...
                                        '02230_2015-02-26_03-16',...
                                        '02230_2015-02-26_09-16',...
                                        '02230_2015-02-28_00-34',...
                                        '02230_2015-02-28_04-34',...
                                        '02230_2015-03-02_11-37',...
                                        '02230_2015-03-04_01-37',...
                                        '02230_2015-03-04_05-37',...
                                        '02230_2015-03-04_09-37',...
                                        '02230_2015-03-06_13-37',...
                                        '02230_2015-03-06_15-37',...
                                        '02230_2015-03-07_15-37',...
                                        '02230_2015-03-07_17-37',...
                                        '02230_2015-03-08_17-37',...
                                        '02230_2015-03-08_17-37',...
                                        '02230_2015-03-09_13-38',...
                                        '02230_2015-03-09_15-38',...
                                        '02230_2015-03-09_17-38',...
                                        '02230_2015-03-10_07-38',...
                                        '02230_2015-03-10_09-38',...
                                        '02230_2015-03-11_07-38'};
                                    

config{2}.directorylist             = [];
config{2}.directorylist{1}          = { 
                                        '02614_2018-06-23_01-23',...
                                        '02614_2018-06-23_03-23',...
                                        '02614_2018-06-23_05-23',...
                                        '02614_2018-06-24_19-23',...
                                        '02614_2018-06-24_23-23',...
                                        '02614_2018-06-25_01-23',...
                                        '02614_2018-06-25_07-23',...
                                        '02614_2018-06-25_10-47',...
                                        '02614_2018-06-25_14-47',...
                                        '02614_2018-06-25_16-47',...
                                        '02614_2018-06-26_04-47'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patient 2, second nodule  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
