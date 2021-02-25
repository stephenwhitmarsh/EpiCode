function wod_concatenateLFP(slurm_task_id)
% code copied from writespykingcircus.m
% Must used a modified version of ft_writedata to give the good channel
% name to Muse. On added line, which is commented with %Paul.

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

if slurm_task_id(1) > 0

for irat = slurm_task_id
    
    % loop through different parts
    for ipart = 1 : size(config{irat}.directorylist,2)
        
        % process channels separately
        for ichan = 1 : size(config{irat}.LFP.channel,2)
            
            if isempty(config{irat}.LFP.channel{ichan})
                continue
            end
            
            [~,foldername] = fileparts(config{irat}.rawdir);
            output_datapath = fullfile(config{irat}.datasavedir,'concatenated_LFP',foldername);
            
            if ~isfolder(output_datapath)
                mkdir(output_datapath);
            end
            
            clear dirdat
            % loop over all directories (time), concatinating channel
            for idir = 1 : size(config{irat}.directorylist{ipart},2)
                clear fname
                temp                                    = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},['*',config{irat}.LFP.channel{ichan},'.ncs']));
                fname{1}                                = fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},temp.name);
                cfgtemp                                 = [];
                cfgtemp.dataset                         = fname{1};
                fprintf('LOADING: %s\n',cfgtemp.dataset);
                clear fname
                fname{1}                                = cfgtemp.dataset;
                dirdat{idir}                            = ft_read_neuralynx_interp(fname);
            end
            
            % concatinate data over files
            chandat = dirdat{1};
            for idir = 2 : length(config{irat}.directorylist{ipart})
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
                chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
            end
            
            if isempty(idir) %correct idir if there is only 1 file
                idir = 1;
            end
            
            % create filename for concatinated data
            temp        = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},['*',config{irat}.LFP.channel{ichan},'.ncs']));
            hdrtemp     = ft_read_header(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir}, temp.name));
            clear fname
            
            fname = fullfile(output_datapath,[config{irat}.prefix,config{irat}.LFP.channel{ichan}, '.ncs']);
                        
            hdr                         = [];
            hdr.Fs                      = hdrtemp.Fs;
            hdr.nSamples                = size(chandat.trial{1},2);
            hdr.nSamplePre              = 0;
            hdr.nChans                  = 1;
            hdr.FirstTimeStamp          = 0;
            hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
            hdr.label{1}                = config{irat}.LFP.channel{ichan};
            ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
            
            clear chandat
            clear dirdat
            
            %create Muse Marker File
            copyfile(config{irat}.muse.templatemarker,fullfile(output_datapath,'Events.mrk'));
            
        end % ichan
    end %ipart
end %irat
end %slurm task id
end %wod_concatenate