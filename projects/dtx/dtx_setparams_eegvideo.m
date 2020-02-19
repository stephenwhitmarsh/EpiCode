
%% Setting parameters DTX project Paul Baudin
%The conversion from Deltamed to '.vhdr' loose infioramtion in the header. Do
%not use the header created by ft_read_header, but the header saved with
%the data (data_header, for each file) 

function [config] = dtx_setparams_eegvideo(config)

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh';
    rootpath_data       = '/network/lustre/iss01/epimicro/rodents/raw/DTX-raw/';
    os                  = 'unix'; 
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\lgi1\DTX-EEG-VIDEO\Analyses_Paul';
    rootpath_data       = '\\lexport\iss01.epimicro\rodents\raw\DTX-EEGRodents-Brainvision\';
    os                  = 'windows';
else
    error('Platform not supported')
end

datasavedir = fullfile(rootpath_analysis, 'data');
imagesavedir = fullfile(rootpath_analysis); 

%% Rodent 1
config{1}.os                        = os;
config{1}.format                    = 'brainvision';
config{1}.types                     = ['macro'];
config{1}.name                      = {'SlowWave','Seizure','InterIctal'};
config{1}.prefix                    = 'DTXEEG12-';
config{1}.muse.startend             = {'SlowWave','Crise_Start'; 'Crise_Start','Crise_End'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{1}.foldername                = 'DTXEEG12-100uM';
config{1}.rawdir                    = fullfile(rootpath_data, config{1}.foldername);
config{1}.datasavedir               = datasavedir;         % where to write data
config{1}.imagesavedir              = imagesavedir;       % where to print images
config{1}.rawlabels                 = {'47','45','48'};
config{1}.labels.macro              = {'M1G','M1D','PtA'};

%config{irat}.directorylist
filelist = dir(config{1}.rawdir);
i=0;

% changer directory list : ajouter filelist. pcq tous dans le même
% directory
for ifile = 1:length(filelist)
    if length(filelist(ifile).name)>3
        fileExtension = filelist(ifile).name(1,length(filelist(ifile).name)-3:length(filelist(ifile).name));
        if strncmp(fileExtension,'.eeg',4) 
            i=i+1;
            config{1}.directorylist{1}{i}          =  filelist(ifile).name;
        end
    end
end
clear filelist
                                    
% config{1}.preproc.bsfilter    = 'yes';
% config{1}.preproc.bsfreq      = [49 51];
% config{1}.preproc.inversedata = -1; 

% config{1}.align.name                = {'SlowWave'};
% config{1}.align.channel             = {'E12.ncs'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{1}.align.flip                = {'no'};
% config{1}.align.abs                 = {'no'};
% config{1}.align.method              = {'max'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{1}.align.filter              = {'lp'};
% config{1}.align.freq                = {5};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{1}.align.hilbert             = {'no'};
% config{1}.align.thresh              = [0];
% config{1}.align.toiplot{1}          = [-1,  1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{1}.align.toiactive{1}        = [-0.5, 0.5];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{1}.align.toibaseline{1}      = [-1, -0.5];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% 
% config{1}.LFP.name                  = {'SlowWave'};
% config{1}.LFP.hpfilter              = 'no';
% config{1}.LFP.hpfreq                = 1;
% config{1}.LFP.resamplefs            = 1000;
% config{1}.LFP.baseline              = 'yes';
% config{1}.LFP.baselinewindow{1}     = [-2, -1];
% config{1}.LFP.slidestep             = 0.01;
% config{1}.LFP.electrodeToPlot       = [6, 11];
% 
% %1, 2 and 3 associated with config.muse.startend 
% config{1}.epoch.toi{1}              = [-5, 25];  
% config{1}.epoch.toi{2}              = [-2, 1];  
% config{1}.epoch.toi{3}              = [1, -2];  
% config{1}.epoch.pad{1}              = 1;
% config{1}.epoch.pad{2}              = 0.5;
% config{1}.epoch.pad{3}              = 0.5;

%% Rodent 2
%% Rodent 1
config{2}.os                        = os;
config{2}.format                    = 'brainvision';
config{2}.types                     = ['macro'];
config{2}.name                      = {'SlowWave','Seizure','InterIctal'};
config{2}.prefix                    = 'DTXEEG_EMG_6-';
config{2}.muse.startend             = {'SlowWave','Crise_Start'; 'Crise_Start','Crise_End'; 'Crise_End','SlowWave'};   % start and end Muse marker
config{2}.foldername                = 'DTXEEG_EMG_6-100uM';
config{2}.rawdir                    = fullfile(rootpath_data, config{2}.foldername);
config{2}.datasavedir               = datasavedir;         % where to write data
config{2}.imagesavedir              = imagesavedir;       % where to print images
config{2}.rawlabels                 = {'52','50','53'};
config{2}.labels.macro              = {'M1G','M1D','PtA'};

%config{irat}.directorylist
filelist = dir(config{2}.rawdir);
i=0;
for ifile = 1:length(filelist)
    if length(filelist(ifile).name)>3
        fileExtension = filelist(ifile).name(1,length(filelist(ifile).name)-3:length(filelist(ifile).name));
        if strncmp(fileExtension,'vhdr',4) 
            i=i+1;
            config{2}.directorylist{1}{i}          =  filelist(ifile).name;
        end
    end
end
clear filelist
                                    
config{2}.preproc.bsfilter    = 'yes';
config{2}.preproc.bsfreq      = [49 51];
config{2}.preproc.inversedata = -1; 

% config{2}.align.name                = {'SlowWave'};
% config{2}.align.channel             = {'E12.ncs'};                                                                                    % pattern to identify channel on which to based peak detection                                                                        % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
% config{2}.align.flip                = {'no'};
% config{2}.align.abs                 = {'no'};
% config{2}.align.method              = {'max'};                                                              % whether to align to max, first-after-zero, or nearest-to-t-zero peak, maxabs {'max','first', 'nearest', 'maxabs'}
% config{2}.align.filter              = {'lp'};
% config{2}.align.freq                = {5};                                                                                  % lowpass filter freq to smooth peak detection (Hz)
% config{2}.align.hilbert             = {'no'};
% config{2}.align.thresh              = [0];
% config{2}.align.toiplot{1}          = [-1,  1];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% config{2}.align.toiactive{1}        = [-0.5, 0.5];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
% config{2}.align.toibaseline{1}      = [-1, -0.5];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
% 
% config{2}.LFP.name                  = {'SlowWave'};
% config{2}.LFP.hpfilter              = 'no';
% config{2}.LFP.hpfreq                = 1;
% config{2}.LFP.resamplefs            = 1000;
% config{2}.LFP.baseline              = 'yes';
% config{2}.LFP.baselinewindow{1}     = [-2, -1];
% config{2}.LFP.slidestep             = 0.01;
% config{2}.LFP.electrodeToPlot       = [6, 11];
% 
% %1, 2 and 3 associated with config.muse.startend 
% config{2}.epoch.toi{1}              = [-5, 25];  
% config{2}.epoch.toi{2}              = [-2, 1];  
% config{2}.epoch.toi{3}              = [1, -2];  
% config{2}.epoch.pad{1}              = 1;
% config{2}.epoch.pad{2}              = 0.5;
% config{2}.epoch.pad{3}              = 0.5;
end


