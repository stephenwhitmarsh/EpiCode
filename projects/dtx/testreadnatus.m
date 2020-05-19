function testreadnatus()

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

cfgtemp = [];
if isunix
    cfgtemp.dataset = '/network/lustre/iss01/charpier/echanges/paul.baudin/Natus_DC/P005alex_0beccbb8-35f1-464c-95b1-97a80e2a9bbc_003.e';
elseif ispc
    cfgtemp.dataset = '\\lexport\iss01.charpier\echanges\paul.baudin\Natus_DC\P005alex_0beccbb8-35f1-464c-95b1-97a80e2a9bbc_003.e';
end
testreadNatus = ft_preprocessing(cfgtemp);

[a, b, ~] = fileparts(cfgtemp.dataset);
save(fullfile(a,[b,'.mat']),'testreadNatus','-v7.3');

