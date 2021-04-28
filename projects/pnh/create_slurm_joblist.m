
%% Create slurm job list
config = pnh_setparams;

fname_slurm_joblist = fullfile(config{1}.datasavedir, 'slurm_job_list.txt');
delete(fname_slurm_joblist);
for ipatient = 1:4
    for ipart = 1 : size(config{ipatient}.directorylist, 2)
        subjdir     = config{ipatient}.prefix(1:end-1);
        partdir     = ['p', num2str(ipart)];
        if ~isfield(config{ipatient}.circus, 'channelname')
            fid = fopen(fname_slurm_joblist, 'a');
            if fid == -1
                error('Could not create/open %s', fname_slurm_joblist);
            end
            dirname     = fullfile(config{1}.datasavedir, subjdir, partdir);
            fprintf(fid,'cd %s; ~/.local/bin/spyking-circus spyking-circus.params -c 28; ~/.local/bin/spyking-circus spyking-circus.params -m converting -c 28; echo DONE!!! \n', dirname);
            fclose(fid);  
        else
            for chandir = unique(config{ipatient}.circus.channelname)
                fid = fopen(fname_slurm_joblist, 'a');
                if fid == -1
                    error('Could not create/open %s', fname_slurm_joblist);
                end
                temp        = strcmp(config{ipatient}.circus.channelname, chandir);
                firstchan   = string(config{ipatient}.circus.channel(find(temp,1,'first')));
                dirname     = fullfile(config{1}.datasavedir, subjdir, partdir, chandir);
                fprintf(fid,'cd %s; /.local/bin/spyking-circus spyking-circus.params -c 28; /.local/bin/spyking-circus spyking-circus.params -m converting -c 28; echo DONE!!! \n', dirname);
                fclose(fid);
            end
        end
    end
end

