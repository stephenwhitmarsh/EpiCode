% testfunction
addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip;
addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/wod'));
addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/CEDMATLAB/CEDS64ML
addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/CEDMATLAB/CEDS64ML/x64
CEDS64LoadLib('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/CEDMATLAB/CEDS64ML');

ft_defaults
feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neuralynx

[config] = wod_setparams([]);
irat = 1;
ipart = 1;
[MuseStruct]               = readCEDMarkers(config{irat}, true);
read_CED_continuous(config{irat},1,1);

hdrtemp = ft_read_header('X:\rodents\raw\TestConversionDeltamed\191014R-BEX_0004.edf');
eventtest = ft_read_event('X:\rodents\raw\TestConversionDeltamed\191014R-BEX_0004.edf');

cfgtemp = [];
cfgtemp.hpfilter              = 'yes';
cfgtemp.hpfreq                = 1;
cfgtemp.hpfiltord             = [];%leave empty : default value
cfgtemp.hpfilttype         = [];%leave empty : default valueco
dat = ft_preprocessing(cfgtemp, data);

figure;
hold;
plot(data.time{1}, data.trial{1}(7,:));
plot(dat.time{1}, dat.trial{1}(7,:));
xlim([-2 2]);
