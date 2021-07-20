function dtx_spikes_slurm_joblist

%find config script :
addpath (genpath('\\lexport\iss01.charpier\analyses\lgi1\Git-Paul\EpiCode\projects\dtx'));
config = dtx_spikes_setparams;
sc_version = '1.0.1_dev';

fname = fullfile(config{1}.datasavedir,'dtx_slurm_job_list.txt');
fid = fopen(fname,'w');

for ipatient = 1:size(config,2)
    for ipart = 1:size(config{ipatient}.directorylist,2)
        for channame = string(unique(config{ipatient}.circus.channelname))
            
            %get patient's information
            subjdir     = config{ipatient}.prefix(1:end-1);
            partdir     = ['p',num2str(ipart)];
            chandir     = char(channame);
            chan_param_idx = find(strcmp(config{ipatient}.circus.channel,channame),1,'first');
            filename    = [config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{chan_param_idx},'.ncs'];
            dirname     = sprintf('/network/lustre/iss01/charpier/analyses/lgi1/DTX-PROBE/data/%s/%s/%s',subjdir,partdir,chandir);
            
            %code to launch SC
            load_sc         = sprintf('module load spyking-circus/%s;', sc_version);
            open_dir        = sprintf('cd %s;', dirname);
            change_dead_1   = 'cp SpykingCircus_artefacts_samples_seizures_removed.dead SpykingCircus_artefacts_samples.dead;';
            launch_sc_1     = sprintf('spyking-circus %s -c 14;', filename);
            
            %code to launch extracting
            change_dead_2   = 'cp SpykingCircus_artefacts_samples_seizures_not_removed.dead SpykingCircus_artefacts_samples.dead;';
            launch_sc_2      = sprintf('spyking-circus %s -m whitening,extracting,fitting -c 14;', filename);
            
            %code to convert
            convert_results = sprintf('spyking-circus %s -m converting;', filename);
            
            %pull all together in one line of the text file
            fprintf(fid,[load_sc,open_dir,change_dead_1,launch_sc_1,change_dead_2,launch_sc_2,convert_results,'\n']);            
            
        end
    end
end

close('all')    