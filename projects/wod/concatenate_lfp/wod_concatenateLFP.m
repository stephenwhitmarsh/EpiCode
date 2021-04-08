function wod_concatenateLFP(rat_list, configscript)
% code copied from writespykingcircus.m
% Must use a modified version of ft_writedata to give the good channel
% name to Muse. Each added line in ft_writedata is commented with %Paul.

% configscript : name of the setparams script, ie 'wod_setparams', or
% 'wod_setparams_32chans';

%retrouver le chemin du dossier EpiCode par rapport à ce script
scriptpath = matlab.desktop.editor.getActiveFilename;
epicodepath = fileparts(fileparts(fileparts(fileparts(scriptpath))));

addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects\wod']))
if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

%/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
%Attention au risque d'écraser des fichiers marqueurs déjà existants
overwriteMuseMarkerFile = false; %false : do not write a new muse marker if one is founed.
%/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

ft_defaults

config = eval(configscript);

feature('DefaultCharacterSet', 'CP1252');

for irat = rat_list
    
    % loop through different parts
    for ipart = 1 : size(config{irat}.directorylist,2)
        
        % process channels separately
        for ichan = 1 : size(config{irat}.LFP.channel,2)
            
            if isempty(config{irat}.LFP.channel{ichan})
                continue
            end
            
            [~,foldername] = fileparts(config{irat}.rawdir);
            output_datapath = fullfile(config{irat}.concatdata_path);
            
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
                header_orig{idir}                       = ft_read_header(fname{1});
            end
            
            % concatinate data over files
            chandat = dirdat{1};
            for idir = 2 : length(config{irat}.directorylist{ipart})
                fprintf('Concatinating directory %d, channel %d\n',idir, ichan);
                chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
                chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
                %concatenate events
            end
            
            if isempty(idir) %correct idir if there is only 1 file
                idir = 1;
            end
            
            % create filename for concatinated data
            temp        = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},['*',config{irat}.LFP.channel{ichan},'.ncs']));
            hdrtemp  = ft_read_header(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir}, temp.name));
            clear fname
            
            fname = fullfile(output_datapath,[config{irat}.prefix,config{irat}.LFP.channel{ichan}, '.ncs']);
                        
            hdr                         = [];
            hdr.Fs                      = hdrtemp.Fs;
            hdr.nSamples                = size(chandat.trial{1},2);
            hdr.nSamplesPre             = hdrtemp.nSamplesPre;
            hdr.nChans                  = 1;
            hdr.FirstTimeStamp          = hdrtemp.FirstTimeStamp;
            hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
            hdr.label{1}                = config{irat}.LFP.channel{ichan};
            hdr.chantype                = hdrtemp.chantype;
            hdr.chanunit                = hdrtemp.chanunit;
            
            ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
            
            clear chandat
            clear dirdat
            
            %create Muse Marker File
            
            add_nev = false;
            if overwriteMuseMarkerFile || ~exist(fullfile(output_datapath,'Events.mrk'), 'file')
                copyfile(config{irat}.muse.templatemarker,fullfile(output_datapath,'Events.mrk'));
                add_nev = true;
            end
            
        end % ichan
        
        %process events
        for idir = 1 : size(config{irat}.directorylist{ipart},2)
            temp                                    = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},'*.nev'));
           
            if size(temp,1) > 1
                error('there should be only one nev file')
            end
            if size(temp,1) == 0
                nev_data{idir} = [];
                continue
            end
            
            fname                                = fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},temp.name);
            nev_data{idir} = read_neuralynx_nev(fname,'eventformat','neuralynx_nev');
            for ievent = 1:size(nev_data{idir}, 1)
                nev_data{idir}(ievent).idir = idir;
            end
            
        end
        
        %concatenate events
        nev_all{1} = nev_data{1};
        for idir = 2 : length(config{irat}.directorylist{ipart})
            fprintf('Concatinating events %d\n',idir);
            if isempty(nev_data{idir})
                continue
            end
            for ievent = 1:size(nev_data{idir},1)
                idx = size(nev_all{1},1) + 1;
                for ifield = string(fieldnames(nev_data{idir})')
                    nev_all{1}(idx).(ifield) = nev_data{idir}(ievent).(ifield);
                end
            end
        end
        
        %save events' MATLAB structure
        fname = fullfile(output_datapath,'Events_nev.mat');
        save(fname, 'nev_all');
        
        %add events to MuseStruct
        if add_nev
            cfgtemp = config{irat};
            [~,dir_name]                  = fileparts(config{irat}.rawdir);
            cfgtemp.rawdir                = fullfile(config{irat}.datasavedir,'concatenated_LFP');
            cfgtemp.directorylist{ipart}  = {dir_name};
            MuseStruct = readMuseMarkers(cfgtemp, true);
            MuseStruct = read_nev_Muse(cfgtemp, MuseStruct, nev_all);
            
            fname = fullfile(cfgtemp.rawdir, cfgtemp.directorylist{1}{1}, 'Events.mrk');
            writeMuseMarkerfile(MuseStruct{1}{1}, fname);
        end
        
    end %ipart
end %irat
end %wod_concatenate