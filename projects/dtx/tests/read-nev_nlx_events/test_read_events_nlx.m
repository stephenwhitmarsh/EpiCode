addpath \\lexport\iss01.charpier\analyses\lgi1\Git-Paul\fieldtrip;
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\external'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\shared'));
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\NeuralynxMatlabImportExport_v6.0.0'));
ft_defaults


feature('DefaultCharacterSet', 'CP1252') % To fix bug for weird character problems in reading neurlynx

%read_neuralynx_nev was modified to add eventstrings. 2 lines are modified,
%I put 'Paul' in comment for each.
% read_neuralynx_nev is in fieldtrip/fileio/private
filename = 'Z:\analyses\lgi1\Git-Paul\EpiCode\projects\dtx\read-nev_nlx_events\test_read_events_nlx.nev';
events = ft_read_event(filename);

