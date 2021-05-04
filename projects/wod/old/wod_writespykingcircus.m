function wod_writespykingcircus(slurm_task_id)

if ispc
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\shared'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\external'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\templates'))
    addpath (genpath('\\lexport\iss01.charpier\analyses\wod\EpiCode\projects\wod'))
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
    
elseif isunix
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/shared'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/external'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/templates'))
    addpath (genpath('/network/lustre/iss01/charpier/analyses/wod/EpiCode/projects/wod'))
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end


ft_defaults

config = wod_setparams;


for irat = slurm_task_id

    %create multifile (concatenated nlx file) and dead file
    writeSpykingCircus(config{irat}, [], true, true);
    %create param and prb file
    writeSpykingCircusParameters(config{irat})

end



end