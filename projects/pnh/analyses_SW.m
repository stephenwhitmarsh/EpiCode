%% Analysis script for Valerio Frazzini's study on electrophysiology of pathoanatomy
%
% (c) Stephen Whitmarsh, stephen.whitmarsh@gmail.com
%
% requires bandpassFilter.m from Mario
% requires releaseDec2015 from Neuralynx website

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy/releaseDec2015/
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/

% addpath /Volumes/iss01.charpier/stephen.whitmarsh/data/02230_2015-02-25_14-36/m1pNs
% addpath /Volumes/iss01.charpier/stephen.whitmarsh/fieldtrip/

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx
maxNumCompThreads(8)

% Setting parameters

% Patient 1, perivetricular heterotopia #1
hpfilter{1}                = {'no','no','no','no'};
hpfreq{1}                  = {1,1,1,1};
label{1}                   = { 'VF1','RR','P','PP'};                                                                      % marker name
startend{1}                = { 'VF1__START__','VF1__END__';'RR__START__','RR__END__';'P','P';'PP__START__','PP__END__'};   % start and end Muse marker
prestim{1}                 = [2, 0.50,   0.1,    0.25];                                                                  % list of onset timing with respect to start-marker (s)
poststim{1}                = [2,   0.50,   0.2,    0.5];                                                                   % list of offset timing with respect to end-marker (s)
pad{1}                     = 0.25;
slidestep{1}               = [0.01,0.001,0.001,0.001];
lpfreq{1}                  = [4,4,40,40];                                                                                  % lowpass filter freq to smooth peak detection (Hz)
toiactive{1}               = [-0.1, 0.8;  -0.1,  0.3;  -0.1, 0.1;   0.0,  0.3];                                            % active period in which to search for peaks [ -0.1,  30;  0, 30;  -0.1, 0.1;0,  0.1];
toibaseline{1}             = [-1.0,-0.1;  -1.0, -0.1;  -1.0 -0.1;  -1.0,  0.3];                                            % baseline period in which to search for peaks [ -1,  0; -1,  0;  -1,  -0.1;  -1, -0.1];
resamplefs{1}              = 640;                                                                                          % resample data to speed up and reduce memory (Hz), and align micro and macro
channel{1}                 = {'m1pNs_4','m1pNs_4','m1pNs_4','m1pNs_4'};                                                                                    % pattern to identify channel on which to based peak detection
alignmethod{1}             = {'first','first','max','first'};                                                              % whether to align to max, first-after-zero, or nearest-to-zero peak {'max','first', or 'nearest'}
pkthresh{1}                = [0,0.25,0.25,0.25];                                                                           % peak threshold: fraction (0:inf) of mean peak amplitude in baseline period
datasavedir{1}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/data';         % where to write data
imagesavedir{1}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/analysis/images';       % where to print images
patient_directory{1}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02230_0674/eeg';
directory_searchstring{1}  = '02230_2015-02-*';
micro_searchstring{1}{1}   = '*m1pNs*.ncs';
macro_searchstring{1}{1}   = '*_1pNs*.ncs';
micro_labels{1,1}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,1}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,2}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,2}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,3}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,3}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
micro_labels{1,4}          = ["m1pNs_1","m1pNs_4","m1pNs_6","m1pNs_7","m1pNs_8"];
macro_labels{1,4}          = ["_1pNs_1","_1pNs_2","_1pNs_3","_1pNs_4","_1pNs_5","_1pNs_6","_1pNs_7","_1pNs_8"];
timwin{1,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,2}                = [-0.1 0.1];      % for plotting spikerate
timwin{1,3}                = [-0.0025 0.0025];  % for plotting spikerate
timwin{1,4}                = [-0.01 0.01];      % for plotting spikerate

% Patient 2, microgyria
hpfilter{2}                = {'no','no','no','no'};
hpfreq{2}                  = {1,1,1,1};
label{2}                   = { 'RR','other'};
startend{2}                = { 'RR__START__','RR__END__';'other__START__','other__END__'};
prestim{2}                 = [0.5,0.5];
poststim{2}                = [2,2];
pad{2}                     = 0.25;
slidestep{2}               = [0.01,0.01];
lpfreq{2}                  = [4,4];
toiactive{2}               = [-0.1, 0.8; -0.1, 0.8];
toibaseline{2}             = [-1.0,-0.1; -1.0,-0.1];
resamplefs{2}              = 640;
channel{2}                 = {'mPlan_2','mPlan_2','mPlan_2','mPlan_2'};
alignmethod{2}             = {'max','max'};
pkthresh{2}                = [0.25, 0.25];
datasavedir{2}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/analysis/data';
imagesavedir{2}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/analysis/images';
patient_directory{2}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02476_0929/eeg';
directory_searchstring{2}  = '02476_2017-04-*';
micro_searchstring{2}{1}   = '*mPlan*.ncs';
macro_searchstring{2}{1}   = '*_Plan*.ncs';
micro_labels{2,1}          = ["mPlan_1","mPlan_2","mPlan_4","mPlan_5","mPlan_6","mPlan_7","mPlan_8"];
macro_labels{2,1}          = ["_Plan_1","_Plan_2","_Plan_3","_Plan_4","_Plan_5","_Plan_6","_Plan_7","_Plan_8"];
micro_labels{2,2}          = ["mPlan_1","mPlan_2","mPlan_4","mPlan_5","mPlan_6","mPlan_7","mPlan_8"];
macro_labels{2,2}          = ["_Plan_1","_Plan_2","_Plan_3","_Plan_4","_Plan_5","_Plan_6","_Plan_7","_Plan_8"];
timwin{2,1}                = [-0.01 0.01]; % for plotting spikerate
timwin{2,2}                = [-0.01 0.01]; % for plotting spikerate

% Patient 3, periventricular heterotopia #2
hpfilter{3}                = {'no','yes','no','no'};
hpfreq{3}                  = {1,1,1,1};
label{3}                   = { 'H02614_01','H02614_02','H02614_03','H02614_04'};
startend{3}                = { 'H02614_01__START__','H02614_01__END__';'H02614_02','H02614_02';'H02614_03__START__','H02614_03__END__';'H02614_04__START__','H02614_04__END__'};
prestim{3}                 = [2, 0.2, 0.5, 0.5];
poststim{3}                = [2,   0.2, 5,   2.5];
pad{3}                     = 0.25;
slidestep{3}               = [0.01,0.001,0.01,0.01];
lpfreq{3}                  = [4,40,4,40];
toiactive{3}               = [-0.1, 0.8; -0.1, 0.1; -0.1, 0.8; -0.1, 0.2];
toibaseline{3}             = [-1.0,-0.1; -1.0,-0.1; -1.0,-0.1; -1.0,-0.2];
resamplefs{3}              = 640;
channel{3}                 = {'mCasd_2','mCasd_2','mCasd_2','mCasd_2'};
alignmethod{3}             = {'first','max','first','first'};
pkthresh{3}                = [0.25, 0.25, 0.25, 0.25];
datasavedir{3}             = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/data';
imagesavedir{3}            = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/analysis/images';
patient_directory{3}       = '/network/lustre/iss01/charpier/stephen.whitmarsh/data/pat_02614_1073/eeg';
directory_searchstring{3}  = '02614_2018-06-*';
micro_searchstring{3}{1}   = '*mCasd*.ncs';
macro_searchstring{3}{1}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*mCasd*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*mCasd*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_searchstring{3}{2}   = '*mCasd*.ncs';
macro_searchstring{3}{2}   = '*_Casd*.ncs';
micro_labels{3,1}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
macro_labels{3,1}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,2}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
macro_labels{3,2}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,3}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
macro_labels{3,3}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
micro_labels{3,4}          = ["mCasd_1","mCasd_2","mCasd_3","mCasd_4","mCasd_5","mCasd_6","mCasd_7","mCasd_8"];
macro_labels{3,4}          = ["_Casd_1","_Casd_2","_Casd_3","_Casd_4","_Casd_5","_Casd_6","_Casd_7","_Casd_8"];
timwin{3,1}                = [-0.1 0.1];      % for plotting spikerate
timwin{3,2}                = [-0.01 0.01];      % for plotting spikerate
timwin{3,3}                = [-0.1 0.1];  % for plotting spikerate
timwin{3,4}                = [-0.1 0.1];      % for plotting spikerate


ipatient = 3;
imarker = 1;


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis project Stephen Whitmarsh %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ipatient = 1 : 1
    

    % read muse markers
    cfg = [];
    cfg.force                   = false;
    cfg.patient_directory       = patient_directory{ipatient};
    cfg.directory_searchstring  = directory_searchstring{ipatient};
    cfg.channel_searchstring    = micro_searchstring{ipatient}{1};
    cfg.datasavedir             = datasavedir{ipatient};
    cfg.suffix                  = [num2str(ipatient),'_micro'];
    MuseStruct_micro            = readMuseMarkers(cfg);
    
    cfg.suffix                  = [num2str(ipatient),'_macro'];
    cfg.channel_searchstring    = macro_searchstring{ipatient}{1};
    MuseStruct_macro            = readMuseMarkers(cfg);
    
    % write concatinated data for spyking-circus
    cfg                                         = [];
    cfg.force                                   = true;
    cfg.forcereload                             = true;
    cfg.hpfilter                                = 'no';
    cfg.hpfreq                                  = 100;
    cfg.channel                                 = micro_labels{ipatient,imarker};
    cfg.datasavedir                             = datasavedir{ipatient};
    [sampleinfo,fnames_ncs,deadfile_ms,deadfile_samples] = writemicroforspykingcircus_alldata(cfg, MuseStruct_micro);     
%    
%     filename = fullfile(cfg.datasavedir,'all_data_artefacts_ms.dead');
%     fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
%     dlmwrite(filename,deadfile_ms,'delimiter','	','precision','%.4f');
%     
%     filename = fullfile(cfg.datasavedir,'all_data_artefacts_samples.dead');
%     fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
%     dlmwrite(filename,deadfile_samples,'delimiter','\t','precision','%.0f');
%     
%     filename = fullfile(cfg.datasavedir,'all_data_dirlist.txt');
%     fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
%     dlmwrite(filename,dirlist);
%     
%     
    
    
    
    
    
    
    
    
    
    % read and plot spike data from Spyking Circus - reads for all markers!
    cfg                     = [];
    cfg.force               = true;
    cfg.imagesavedir        = imagesavedir{ipatient};
    cfg.datasavedir         = datasavedir{ipatient};
    cfg.channel             = channel{ipatient}{imarker};
    cfg.resamplefs          = resamplefs{ipatient};
    cfg.sampleinfo           = sampleinfo;
    cfg.timwin              = [-0.1 0.1];
    cfg.deadfile_samples    = deadfile_samples;
    cfg.fnames_ncs          = fnames_ncs;
    cfg.suffix              = '-1';
    [SpikeRaw, SpikeTrials] = readSpykingCircus_alldata(cfg);          
        
        %% cross correlation
        cfg = [];
        cfg.binsize = 0.001;
        cfg.maxlag = 0.020;
        cfg.debias = 'yes';
        cfg.method = 'xcorr';
        stat_x = ft_spike_xcorr(cfg,SpikeTrials);
        
        
        %         cfg.method = 'shiftpredictor';
        %         stat_s      = ft_spike_xcorr(cfg,SpikeTrials{imarker});
        %         stat_s.xcorr = stat_s.shiftpredictor;
        %
        %         stat_diff = stat_x;
        %         stat_diff.xcorr = stat_x.xcorr - stat_s.shiftpredictor;
        
        % plot crosscorrelation
        cfg = [];
        cfg.imagesavedir        = imagesavedir{ipatient};
        cfg.prefix = label{ipatient}{imarker};
        plotxcorr(cfg,stat_x);
        
        
        figure;
        
        bar(stat_x.time,squeeze(stat_x.xcorr(4,7,:)),1,'facecolor',[0 0 0],'edgecolor',c);
        
    end
        
        %         cfg.prefix = '-shiftpredictor';
        %         plotxcorr(cfg,stat_s);
        
        
        cfgtemp                     = [];
        cfgtemp.force               = false;
        cfgtemp.datasavedir         = datasavedir{ipatient};
        cfgtemp.label               = label{ipatient}{imarker};
        dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);
        
        
        
        
        % add header - better to do it in readMicroFs
        temp                        = dir(fullfile(datasavedir{ipatient},['all_concatinated_',channel{ipatient}{imarker}(1:end-2),'_*.ncs']));
        hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
        dat_microFs.hdr             = hdr;
        %         trialinfo_allmarkers
        
        % recreate trl
        trl                         = trialinfo_single{imarker} + trialinfo_allmarkers(imarker,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
        trl(:,3)                    = -ones(size(trl,1),1) * prestim{ipatient}(imarker) * hdr.Fs;
        dat_microFs.cfg.trl         = trl;
        
        % combine spike and LFP
        data_all = ft_appendspike([],dat_microFs,SpikeTrials{imarker});
        
        % plot spikes together with lfp
        cfg                 = [];
        cfg.viewmode        = 'vertical';
        cfg.continuous      = 'yes';
        cfg.mychan          = ft_channelselection('temp*', data_all, 'all');
        cfg.mychanscale     = repmat(max(max(data_all.trial{1}))*2,size(cfg.mychan));   
        cfg.channelcolormap = [0,0,0;1,0,0];
        cfg.colorgroups     = [1,1,1,1,1,2,2,2,2,2]; 
        cfg.linewidth       = 1;        
        ft_databrowser(cfg,data_all);
     
        
%         cfg.channel = {'temp_4','m1pNs_8'};
%         cfg.viewmode = 'butterfly';

        ft_databrowser(cfg,data_all);
        
        % interpolate LFP around spikes
        cfg                         = [];
        cfg.method                  = 'nan';
        cfg.timwin                  = [-0.002 0.002];
        cfg.spikechannel            = 'temp_2';
        cfg.channel                 = {'m1pNs_1','m1pNs_4','m1pNs_6','m1pNs_7','m1pNs_8'};
        data_i                      = ft_spiketriggeredinterpolation(cfg,data_all);
        
        cfg                         = [];
        cfg.timwin                  = [-0.25 0.25];
        cfg.spikechannel            = 'temp_2';
        cfg.channel                 = cellstr(micro_labels{1}); % first four chans
        cfg.latency                 = [-prestim{ipatient}(imarker) poststim{ipatient}(imarker)];
        spikeavg                    = ft_spiketriggeredaverage(cfg, data_i);
        
        figure;
        plot(spikeavg.time,spikeavg.avg(3,:))
        
        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.foi = 20:10:100;
        cfg.t_ftimwin = 5./cfg.foi;
        cfg.taper = 'hanning';
        stsConcol = ft_spiketriggeredspectrum(cfg,data_all);
        
        
        
        
    end
    close all
    
end % ipatient


% ilabel = 1;
%
% cfgtemp                         = [];
% cfgtemp.trl                     = cfg.trialinfo_single{ilabel} + cfg.trialinfo_all(ilabel,1) - 1; % + cfg.prestim(ilabel) * hdr.Fs;
% cfgtemp.trl(:,3)                = -ones(size(cfgtemp.trl,1),1) * cfg.prestim(ilabel) * hdr.Fs;
% cfgtemp.trlunit                 = 'samples';
% cfgtemp.hdr                     = hdr;
% SpikeTrials                     = ft_spike_maketrials(cfgtemp,SpikeRaw);


cfgtemp                     = [];
cfgtemp.force               = false;
cfgtemp.datasavedir         = datasavedir{ipatient};
cfgtemp.label               = label{ipatient}{ilabel};
dat_microFs                 = readMicroFs(cfgtemp,MuseStruct_micro);

temp                        = dir(fullfile(cfg.datasavedir,['all_concatinated_',cfg.channel(1:end-2),'_*.ncs']));
hdr                         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data

dat_microFs.hdr = hdr;
data_all = ft_appendspike([],dat_microFs,SpikeTrials);


cfgtemp                             = [];
cfgtemp.timwin                      = [-0.25 0.25];
cfgtemp.spikechannel                = channel{ipatient}{ilabel}; % first unit
cfgtemp.channel                     = micro_labels{ipatient,ilabel}; % first four chans
cfgtemp.latency                     = [-prestim{ipatient}(ilabel) poststim{ipatient}(ilabel)];
staPost                             = ft_spiketriggeredaverage(cfg, SpikeTrials);

dat_micro
%% write deadfile for Spyking Circus
% writeSpykingCircusDeadFile(MuseStruct_micro_aligned)
% writeSpykingCircusDeadFile_concatinated(MuseStruct_micro_aligned)


