restoredefaultpath
if ispc
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\shared'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\external'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\templates'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'))
    addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\SPIKY_apr_2021'))
    %     addpath (genpath('\\lexport\iss01.charpier\analyses\vn_preictal\scripts\SPIKY_Dec_2019'))

    addpath \\l2export\iss02.charpier\analyses\vn_preictal\scripts\fieldtrip-20200607

elseif isunix
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/external'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/EpiCode/projects/preictal'))
    addpath (genpath('/network/lustre/iss02/charpier/analyses/vn_preictal/scripts/SPIKY_apr_2021'))
    %     addpath (genpath('/network/lustre/iss01/charpier/analyses/vn_preictal/scripts/SPIKY_Dec_2019'))
    addpath /network/lustre/iss02/charpier/analyses/vn_preictal/scripts/fieldtrip-20200607
end
ft_defaults

config = preictal_setparams;

preictal_spikes_slurm_joblist


for elec = 1 %[1:4] => patients 1 2 3 et 4, [1 4] => 1 et 4, [1 4 7] => patients 1, 4 et 7

    % read muse markers
    % MuseStruct{PART}{DIRECTORY}, e:g MuseStruct{1}{1}
    MuseStruct = readMuseMarkers(config{elec}, true);

    % remove post ictal from the whole analysis,
    % according to config (some 'patients' will have
    % shorter postictal kept because of noise, see setparams)
	MuseStruct                  = addMuseBAD(config{ipatient},MuseStruct);

	%save deadfile without removing seizures
	writeSpykingCircusDeadFile(config{ipatient}, MuseStruct, '_SeizuresNotRemoved');

	%remove seizures for the Spyking-Circus clustering
	cfgtemp                       = [];
	cfgtemp.bad.markerStart       = 'CriseStart';
	cfgtemp.bad.markerEnd         = 'CriseEnd';
	MuseStruct                    = addMuseBAD(cfgtemp,MuseStruct);

	writeSpykingCircusDeadFile(config{ipatient}, MuseStruct, '_SeizuresRemoved');

    % write parameters file for spyking circus
    writeSpykingCircusParameters(config{elec});

    % write file list for spyking circus
    writeSpykingCircusFileList(config{elec}, true);

end
