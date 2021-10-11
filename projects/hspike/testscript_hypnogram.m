cd \\lexport\iss01.charpier\analyses\stephen.whitmarsh\EpiCode\projects\hspike
hspike_setpaths;

% set all parameters
config = hspike_setparams;

% replace data pointer to temporary data
config{2}.rawdir = '\\lexport\iss01.epimicro\patients\shared\tmp\test\pat_02718_1201\eeg';
config{2}.directorylist{1}          =  {'02718_2019-05-14_20-31',... % very artefacted
    '02718_2019-05-14_22-31',... % MISSING DATA
    '02718_2019-05-15_00-31',...
    '02718_2019-05-15_02-31',...
    '02718_2019-05-15_04-31',...
    '02718_2019-05-15_06-31'};
config{2}.directorylist{2}          =  {'02718_2019-05-15_23-52',...
    '02718_2019-05-16_01-52',...
    '02718_2019-05-16_03-52',...
    '02718_2019-05-16_05-52',...
    '02718_2019-05-16_07-52'};   % very artefacted for most of the recording
config{2}.directorylist{3}          =  {'02718_2019-05-16_23-15',...
    '02718_2019-05-17_01-15',...
    '02718_2019-05-17_03-15',...
    '02718_2019-05-17_05-15',...
    '02718_2019-05-17_07-15'};

config{4}.rawdir = '\\lexport\iss01.epimicro\patients\shared\tmp\test\pat_02680_1158\eeg';
config{4}.directorylist{1}          = { '02680_2019-01-15_21-31'...
    '02680_2019-01-15_23-31'...
    '02680_2019-01-16_01-31'...
    '02680_2019-01-16_03-31'...
    '02680_2019-01-16_05-31'}; % ChWrongSize
config{4}.directorylist{2}          = { '02680_2019-01-16_23-32'... % no signal - flatline - looks hp filtered
    '02680_2019-01-17_01-32'... % no signal - flatline - looks hp filtered
    '02680_2019-01-17_03-32'... % no signal - flatline - looks hp filtered
    '02680_2019-01-17_05-32'... % no signal - flatline - looks hp filtered
    '02680_2019-01-17_07-32'};
config{4}.directorylist{3}          = { '02680_2019-01-18_00-01'...
    '02680_2019-01-18_02-01'...
    '02680_2019-01-18_04-01'...
    '02680_2019-01-18_06-01'...
    '02680_2019-01-18_08-01'};

% read markers from (temporary) data
MuseStruct{ipatient} = readMuseMarkers(config{ipatient}, true);

% create hypnogram
[~, hypnogram{ipatient}, ~] = hypnogramMuseStats(config{ipatient}, MuseStruct{ipatient}, true);

% plot hypnogram
cfg                 = [];
cfg.ipart           = 1;
cfg.filename        = fullfile(config{ipatient}.imagesavedir, 'spikeoverview', sprintf('spikeoverview_s%sp%d.jpg', config{ipatient}.prefix, ipart));
cfg.xlim            = [MuseStruct{ipatient}{cfg.ipart}{1}.starttime, MuseStruct{ipatient}{cfg.ipart}{end}.endtime];
cfg.orientation     = 'landscape';
cfg.type{1}         = 'hypnogram';
cfg.title{1}        = 'Hypnogram';
cfg.PSGcolors(1)    = false;
plotWindowedData(cfg, MuseStruct{ipatient}, hypnogram{ipatient})

        
% calculate differences etc.
ipatient = 2; %4
t = (hypnogram{ipatient}.endtime(1:end-1) - hypnogram{ipatient}.starttime(2:end));
t(t ~= 0)
hypnogram{ipatient}(t ~= 0, :)
