if ispc
    addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
    addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/fieldtrip/
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/projects/dtx'));
    addpath(genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/external'));
    addpath /network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
    
end

ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = dtx_setparams_probe_spikes([]);
cfg = config{1};
ipart = 1;

datadir             = fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)]);
phydir_temp         = fullfile(datadir,'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*',cfg.circus.postfix,'.GUI']);
temp                = dir(phydir_temp);
if isempty(temp)
    error('Could not find Phy-converted Spyking-Circus results: %s\n',phydir);
end

cfgtemp = [];
cfgtemp.phy_directory              = fullfile(temp.folder,temp.name);
datafile            = [temp.name(1 : end - length(cfg.circus.postfix) - length('.GUI')), '.ncs'];
hdr_fname                        = fullfile(datadir,datafile);
cfgtemp.hdr                              = ft_read_header(hdr_fname);

spikeraw_test = read_spike_Phy(cfgtemp);
