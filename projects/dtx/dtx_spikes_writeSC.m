function dtx_spikes_writeSC(ipatient)

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
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/EpiCode/shared'));
    addpath (genpath('/network/lustre/iss01/charpier/analyses/lgi1/Git-Paul/NeuralynxMatlabImportExport_v6.0.0'));
end

ft_defaults
%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig   = 'silent';
ft_default.checkpath     = 'once';
ft_default.showcallinfo  = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

config = dtx_spikes_setparams;

MuseStruct                       = readMuseMarkers(config{ipatient}, true);
if strcmp(config{ipatient}.type, 'dtx') 
    MuseStruct                   = dtx_remove_wrong_seizure(config{ipatient}, MuseStruct,true);
end

% write the 2 different dead files
writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct, true, '_seizures_not_removed');
MuseStruct = addMuseBAD(config{ipatient},MuseStruct);
writeSpykingCircusDeadfiles(config{ipatient}, MuseStruct, true, '_seizures_removed');

writeSpykingCircusParameters(config{ipatient});
writeSpykingCircus(config{ipatient}, true);
  