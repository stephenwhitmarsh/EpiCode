
%% Setting parameters

function [config] = seg_setparams(config)

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

% Patient 1, perivetricular heterotopia #1
disp('setting parameters');

% patient 1
config{1}.os                        = os;
config{1}.prefix                    = 'P1-';
config{1}.type                      = 'parts';
config{1}.name                      = {'Crise'};
config{1}.muse.startend             = {'CriseStart','CriseEnd'};   % start and end Muse marker
config{1}.patientdir                = fullfile(rootpath_data, 'pat_02599_1057');
config{1}.rawdir                    = fullfile(rootpath_data, 'pat_02599_1057', 'eeg');
config{1}.datasavedir               = fullfile(rootpath_analysis, 'data',   'seg');         % where to write data
config{1}.imagesavedir              = fullfile(rootpath_analysis, 'images', 'seg');       % where to print images
config{1}.labels.macro              = {'_HaT2_6'};                                                                  % at least one macro for calculating sample rate etc.
config{1}.labels.micro              = {'mHaT2_6'};                                                                  % at least one micro for calculating sample rate etc.
% config{1}.directory_searchstring    = '*';
config{1}.directorylist{1}          = {'02599_2018-04-24_18-36',...
                                       '02599_2018-04-24_20-36'};                                   
config{1}.directorylist{2}          = {'02599_2018-04-25_11-24',...
                                       '02599_2018-04-25_13-24'};                                   
config{1}.directorylist{3}          = {'02599_2018-04-25_13-24',...
                                       '02599_2018-04-25_15-24'};                                   
config{1}.directorylist{4}          = {'02599_2018-04-27_13-23',...
                                       '02599_2018-04-27_15-23'};                                
config{1}.directorylist{5}          = {'02599_2018-04-28_11-23',...
                                       '02599_2018-04-28_13-23'};
                                   
config{1}.circus.channel            = {'mHaT2_6'};
config{1}.circus.reref              = 'no';
config{1}.circus.refchan            = 'xxx';
config{1}.circus.outputdir          = fullfile(rootpath_analysis, 'data', 'seg', 'SpykingCircus');
config{1}.circus.hpfilter           = 'no'; % hp before writing data for SC, does not change the hp of SC
config{1}.epoch.toi{1}              = [-10,  60];           
config{1}.stats.actoi{1}            = [-10, 60];
config{1}.stats.bltoi{1}            = [-10, -5];
config{1}.stats.alpha               = 0.025;
config{1}.spike.slidestep           = [0.01];
config{1}.spike.toispikerate{1}     = [-0.1 0.1];           % for plotting spikerate
config{1}.spike.resamplefs          = 1000;
config{1}.spike.bltoi{1}            = [-10, -5];

% patient 2
config{2}.prefix                    = 'seg_2256_';
config{2}.datasavedir               = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/seg';
config{2}.imagesavedir              = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/images/seg';
config{2}.rawdir                    = '/network/lustre/iss01/epimicro/patients/raw/pat_02256_0700/eeg/';
config{2}.labels.macro              = {'_Cinp_1'};
config{2}.labels.micro              = {'mCinp_1'};
config{2}.directory_searchstring    = '*';

% patient 3
config{3}.prefix                    = 'seg_2379_';
config{3}.datasavedir               = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/seg';
config{3}.imagesavedir              = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/images/seg';
config{3}.rawdir                    = '/network/lustre/iss01/epimicro/patients/raw/pat_02379_0828/eeg/';
config{3}.labels.macro              = {'_Ha2d_1'};
config{3}.labels.micro              = {'mHa2d_1'};
config{3}.directory_searchstring    = '*';

% patient 4
config{4}.prefix                    = 'seg_2614_';
config{4}.datasavedir               = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/data/seg';
config{4}.imagesavedir              = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/images/seg';
config{4}.rawdir                    = '/network/lustre/iss01/epimicro/patients/raw/pat_02614_1073/eeg/';
config{4}.labels.macro              = {'_Casd_1'};
config{4}.labels.micro              = {'mCasd_1'};
config{4}.directory_searchstring    = '*';

end
