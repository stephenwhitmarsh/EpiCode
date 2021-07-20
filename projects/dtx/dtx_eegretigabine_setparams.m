function [config] = dtx_eegretigabine_setparams(config)
%The conversion from Deltamed to '.vhdr' loose information in the header. Do
%not use the header read by ft_read_header, but the header saved with
%the data (data_header, for each file) 

disp('setting parameters for EEG-retigabine');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/lgi1/DTX-EEG-RETIGABINE/';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-EEGRetigabine-Brainvision/';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-EEG-RETIGABINE';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-EEGRetigabine-Brainvision\';
    os                  = 'windows';
else
    error('Platform not supported')
end

%% Config common for all rats
datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis,'image');

configcommon.name                      = {''};
configcommon.datasavedir               = datasavedir;         

configcommon.seizuretimings.marker_start = 'Crise_Start';
configcommon.seizuretimings.marker_end   = 'Crise_End'; 
configcommon.seizuretimings.analysis_start = 'Analysis_Start';
configcommon.seizuretimings.analysis_end   = 'Analysis_End';
configcommon.seizuretimings.winsize        = 3600;%s
configcommon.seizuretimings.winstep        = 1300;%s

%% Rodent 1
config{1}                           = configcommon;
config{1}.prefix                    = 'Rat-2021_04_07-1-';
config{1}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_07-1');
config{1}.rawlabels.oldnames        = {'47','45','48','46','49'};
config{1}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{1}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_07-1');
config{1}.directorylist             = {}; %found automatically

config{1}.injectiontime             = datetime('06-Apr-2021 10:19:00');
config{1}.injectionretig            = datetime('06-Apr-2021 15:22:00');
config{1}.group                     = 'retigabine';

%% Rodent 2
config{2}                           = configcommon;
config{2}.prefix                    = 'Rat-2021_04_07-2-';
config{2}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_07-2');
config{2}.rawlabels.oldnames        = {'52','50','53','51','54'};
config{2}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{2}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_07-2');
config{2}.directorylist             = {}; %found automatically

config{2}.injectiontime             = datetime('06-Apr-2021 11:48:00');
config{2}.injectionretig            = datetime('06-Apr-2021 15:48:00');
config{2}.group                     = 'dmso';

%% Rodent 3
config{3}                           = configcommon;
config{3}.prefix                    = 'Rat-2021_04_07-3-';
config{3}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_07-3');
config{3}.rawlabels.oldnames        = {'57','55','58','56','59'};
config{3}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{3}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_07-3');
config{3}.directorylist             = {}; %found automatically

config{3}.injectiontime             = datetime('06-Apr-2021 14:45:00');
config{3}.injectionretig            = datetime('06-Apr-2021 18:55:00');
config{3}.group                     = 'retigabine';

%% Rodent 4
config{4}                           = configcommon;
config{4}.prefix                    = 'Rat-2021_04_09-1-';
config{4}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_09-1');
config{4}.rawlabels.oldnames        = {'52','50','53','51','54'};
config{4}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{4}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_09-1');
config{4}.directorylist             = {}; %found automatically

config{4}.injectiontime             = datetime('09-Apr-2021 9:32:00');
config{4}.injectionretig            = datetime('09-Apr-2021 14:13:00');
config{4}.group                     = 'dmso';

%% Rodent 5
config{5}                           = configcommon;
config{5}.prefix                    = 'Rat-2021_04_09-2-';
config{5}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_09-2');
config{5}.rawlabels.oldnames        = {'47','45','48','46','49'};
config{5}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{5}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_09-2');
config{5}.directorylist             = {}; %found automatically

config{5}.injectiontime             = datetime('09-Apr-2021 12:23:00');
config{5}.injectionretig            = datetime('09-Apr-2021 16:55:00');
config{5}.group                     = 'retigabine';

%% Rodent 6
config{6}                           = configcommon;
config{6}.prefix                    = 'Rat-2021_04_09-3-';
config{6}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_09-3');
config{6}.rawlabels.oldnames        = {'57','58'};
config{6}.rawlabels.newnames        = {'M1G','PtA'};
config{6}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_09-3');
config{6}.directorylist             = {}; %found automatically

config{6}.injectiontime             = datetime('09-Apr-2021 15:07:00');
config{6}.injectionretig            = datetime('09-Apr-2021 20:04:00');
config{6}.group                     = 'dmso';

%% Rodent 7
config{7}                           = configcommon;
config{7}.prefix                    = 'Rat-2021_04_16-1-';
config{7}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_16-1');
config{7}.rawlabels.oldnames        = {'52','50','53','51','54'};
config{7}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{7}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_16-1');
config{7}.directorylist             = {}; %found automatically

config{7}.injectiontime             = datetime('16-Apr-2021 12:33:00');
config{7}.injectionretig            = datetime('16-Apr-2021 17:35:00');
config{7}.group                     = 'retigabine';

%% Rodent 8
config{8}                           = configcommon;
config{8}.prefix                    = 'Rat-2021_04_16-2-';
config{8}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_16-2');
config{8}.rawlabels.oldnames        = {'57','55','58','56','59'};
config{8}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{8}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_16-2');
config{8}.directorylist             = {}; %found automatically

config{8}.injectiontime             = datetime('16-Apr-2021 15:29:00');
config{8}.injectionretig            = datetime('16-Apr-2021 20:08:00');
config{8}.group                     = 'retigabine';

%% Rodent 9
config{9}                           = configcommon;
config{9}.prefix                    = 'Rat-2021_04_16-3-';
config{9}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_16-3');
config{9}.rawlabels.oldnames        = {'47','45','48','46','49'};
config{9}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{9}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_16-3');
config{9}.directorylist             = {}; %found automatically

config{9}.injectiontime             = datetime('16-Apr-2021 16:54:00');
config{9}.injectionretig            = datetime('16-Apr-2021 21:07:00');
config{9}.group                     = 'dmso';

%% Rodent 10
config{10}                           = configcommon;
config{10}.prefix                    = 'Rat-2021_04_17-1-';
config{10}.rawdir                    = fullfile(rootpath_data,'Rat-2021_04_17-1');
config{10}.rawlabels.oldnames        = {'52','50','53','51','54'};
config{10}.rawlabels.newnames        = {'M1G','M1D','PtA', 'EMG1', 'EMG2'};
config{10}.imagesavedir              = fullfile(imagesavedir,'Rat-2021_04_17-1');
config{10}.directorylist             = {}; %found automatically

config{10}.injectiontime             = datetime('17-Apr-2021 10:50:00');
config{10}.injectionretig            = datetime('17-Apr-2021 17:26:00');
config{10}.group                     = 'dmso';

%% find data files
for irat = 1:size(config,2)
    filelist = dir(config{irat}.rawdir);
    filelist = natsort({filelist.name});
    i=0;
    for ifile = 1:length(filelist)
        [~,~,file_extension] = fileparts(filelist{ifile});
        if strncmp(file_extension,'.eeg',4)
            i=i+1;
            config{irat}.directorylist{1}{i} =  filelist{ifile}(1:end-4);
        end
    end
    clear filelist
end