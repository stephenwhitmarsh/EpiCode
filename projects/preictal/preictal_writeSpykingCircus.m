function preictal_writeSpykingCircus(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\EpiCode\projects\preictal'))
    addpath \\lexport\iss01.charpier\analyses\louis.cousyn\Scripts\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/EpiCode/projects/preictal'))
    addpath /network/lustre/iss01/charpier/analyses/louis.cousyn/Scripts/fieldtrip-20200607
end


ft_defaults

config = preictal_setparams;

%% create Spyking circus files

for ipatient = slurm_task_id%1 %[2 4] %1:6

%read muse markers
MuseStruct = readMuseMarkers(config{ipatient}, true);

%remove the end of the file, after the seizure to analyse (because we do
%not use the post-ictal data, so better to remove it of the spike sorting)
cfgtemp                     = [];
cfgtemp.bad.markerStart     = 'CriseEnd'; %BAD à partir crise end
cfgtemp.bad.markerEnd       = 'end'; % BAD jusque fin du fichier
cfgtemp.bad.dir_list        = 'last'; %nouveau marqueur BAD sur dernier des 2 fichiers 
cfgtemp.bad.sample_list     = ft_getopt(config{ipatient},'seizure_index', 'last'); %index of the seizure to analyze, on the LAST dir. can be 'last' (default)
cfgtemp.bad.time_from_begin = 60; %début à +60s de crise END (pour rfaire spike sorting sur crise et post critique immédiat)
MuseStruct                  = addMuseBAD(cfgtemp,MuseStruct);

%save deadfile without removing seizures
writeSpykingCircusDeadFile(config{ipatient},MuseStruct);

%remove seizures for the Spyking-Circus clustering
cfgtemp                       = [];
cfgtemp.bad.markerStart       = 'CriseStart';
cfgtemp.bad.markerEnd         = 'CriseEnd';
MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);
 
%create multifile and dead file :
writeSpykingCircus(config{ipatient},MuseStruct,true, true);

%create params file and prb file :
writeSpykingCircusParameters(config{ipatient});


end



