function PET_slurm_joblist

%find config script :
addpath (genpath('\\l2export\iss02.charpier\analyses\vn_preictal\scripts\EpiCode\projects\preictal'));
config = PET_setparams;
sc_version = '1.0.8';

fname = fullfile(config{1}.datasavedir,'spyking_circus_slurm_job_list.txt');
fid = fopen(fname,'w');

for ielec = 1:size(config,2)
    
    for ipart = 1:size(config{ielec}.directorylist,2)
        
        %get patient's information
        subjdir     = config{ielec}.prefix(1:end-1);
        partdir     = ['p',num2str(ipart)];
%         filename    = [config{ipatient}.prefix,'p',num2str(ipart),'-multifile-',config{ipatient}.circus.channel{1},'.ncs'];
        filename    = 'SpykingCircus.params';
        
        dirname     = sprintf('/network/lustre/iss02/charpier/analyses/vn_preictal/data/%s/%s',subjdir,partdir);
        
        %code to launch SC
        load_sc         = sprintf('module load spyking-circus/%s;', sc_version);
        open_dir        = sprintf('cd %s;', dirname);
        launch_sc_1     = sprintf('spyking-circus %s -c 28;', filename);
        
        %code to convert
        convert_results = sprintf('spyking-circus %s -m converting;', filename);
        
        %pull all together in one line of the text file
        fprintf(fid,[load_sc,open_dir,launch_sc_1,convert_results,'\n']);
        
    end
end

close('all')
