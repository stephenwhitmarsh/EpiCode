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

% if isunix
%     cfgtemp.dataset = '/network/lustre/iss01/charpier/echanges/paul.baudin/Natus_DC/P005alex_0beccbb8-35f1-464c-95b1-97a80e2a9bbc_003.e';
% elseif ispc
%     cfgtemp.dataset = '\\lexport\iss01.charpier\echanges\paul.baudin\Natus_DC\P005alex_0beccbb8-35f1-464c-95b1-97a80e2a9bbc_003.e';
% end

fname = 'C:\Users\paul.baudin\Pictures\Diapos LGI1\Schwa\SCHWARTZ__ROSE_B429F9B3_DEDB_11E0_B53F_0017317D9694\data.eeg';

cfgtemp = [];
cfgtemp.dataset = 'C:\Users\paul.baudin\Pictures\Diapos LGI1\Schwa\SCHWARTZ__ROSE_B429F9B3_DEDB_11E0_B53F_0017317D9694\data.eeg';
testreadNatus = ft_preprocessing(cfgtemp);

fid = fopen(fname);

 header.originfilename   = fname;
 header.filename         = fname;
 header.date             = datetime(2011,14,09,15,24,03);
 header.dateEnd          = datetime(2011,14,09,15,36,05);
 header.orig.PID         = 'pat_lgi1_020';
 header.nChans           = 1;
 header.chanunit         = 'uV';
 header.chantype         = 'eeg';
 header.NumberOfChannels = 1;
 header.DataFormat = 'binary';
 header.label = {'Fz'};
dat = read_brainvision_eeg(fname, header, 0, 10000, 1);

  fid = fopen(fname, 'rb', 'ieee-le');
  dat = fread(fid,'float32');
  dat = dat';
  plot(1:10000, dat(1:10000));
  
fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1), 'cof');
dat = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], sampletype);

[a, b, ~] = fileparts(cfgtemp.dataset);
save(fullfile(a,[b,'.mat']),'testreadNatus','-v7.3');

