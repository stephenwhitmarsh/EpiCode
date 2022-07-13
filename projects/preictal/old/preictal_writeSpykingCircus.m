function preictal_writeSpykingCircus(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'))
    addpath \\lexport\iss01.charpier\analyses\vn_preictal\scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/EpiCode/projects/preictal'))
    addpath /network/lustre/iss01/charpier/analyses/vn_preictal/scripts/fieldtrip-20200607
end


ft_defaults

config = preictal_setparams;

%% create Spyking circus files

for ipatient = slurm_task_id%1 %[2 4] %1:6

%read muse markers
MuseStruct = readMuseMarkers(config{ipatient}, true);

%remove the end of the file, after the seizure to analyse (because we do
%not use the post-ictal data, so better to remove it of the spike sorting)
% Paramètres indiqués dans preictal_setparams
MuseStruct                  = addMuseBAD(config{ipatient},MuseStruct);

%save deadfile without removing seizures
config{ipatient}.circus.deadfilesuffix = '_SeizuresNotRemoved';
writeSpykingCircusDeadFile(config{ipatient},MuseStruct);

%remove seizures for the Spyking-Circus clustering
cfgtemp                       = [];
cfgtemp.bad.markerStart       = 'CriseStart';
cfgtemp.bad.markerEnd         = 'CriseEnd';
MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);
 
%create multifile and dead file :
config{ipatient}.circus.deadfilesuffix = '_SeizuresRemoved';
writeSpykingCircus(config{ipatient},MuseStruct,true, true);

%create params file and prb file :
writeSpykingCircusParameters(config{ipatient});


end



