addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

% ft_read_event modifié : commenté avec %Paul
filename = '\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx\tests\read-nev_nlx_events\test_read_events_nlx.nev';
events = ft_read_event(filename);

