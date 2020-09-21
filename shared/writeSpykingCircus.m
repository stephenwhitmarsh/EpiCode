function [sampleinfo] = writeSpykingCircus(cfg, MuseStruct, force, varargin)

% WRITESPYKINGCIRCUS writes a concatinated data file for SpykingCircus and
% sampleinfo is returned for later (cutting results back up)
%
% use as
%   writeSpykingCircus(cfg, MuseStruct, force, varargin)
%
% third argument "force" is to force recalculation of samplinfo - needed in spike analysis
% pipelines
% fourth optional argument is to force rewriting the data (takes a lot of
% time)

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%   EpiCode is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   EpiCode is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with EpiCode. If not, see <http://www.gnu.org/licenses/>.


if isempty(varargin)
    write = true;
else
    write = varargin{1};
end

cfg.circus.writedeadfile    = ft_getopt(cfg.circus, 'writedeadfile', 'yes');
cfg.circus.hpfilter         = ft_getopt(cfg.circus, 'hpfilter', 'no');
cfg.circus.hpfreq           = ft_getopt(cfg.circus, 'hpfilter', 'no');
fname_output                = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_trialinfo_parts.mat']);

if exist(fname_output, 'file') && force == false
    fprintf('\nLoading sampleinfo: %s \n', fname_output);
    temp         = load(fname_output, 'cfg');
    sampleinfo   = temp.cfg.sampleinfo;
    %     cfg.deadfile_ms         = temp.cfg.deadfile_ms;
    %     cfg.deadfile_samples    = temp.cfg.deadfile_samples;
else
    
    % loop through different parts
    for ipart = 1 : size(cfg.directorylist, 2)
        
        fprintf('\n*** Starting on part %d ***\n',ipart)
        
        % process channels separately
        for ichan = 1 : size(cfg.circus.channel,2)
            
            sampleinfo{ipart}{ichan} = [];
            clear dirdat
            
            % loop over all directories (time), concatinating channel
            for idir = 1 : size(cfg.directorylist{ipart},2)
                
                %
                %                 % NEED SOMETHING LIKE THIS
                %                 % select MICRO files
                %                 micro_filenrs = [];
                %                 for ifile = 1 : size(MuseStruct_micro{ipart}{idir}.filenames,2)
                %                     for ilabel = 1 : size(cfg.labels.micro,2)
                %                         if ~isempty(strfind(MuseStruct_micro{ipart}{idir}.filenames{ifile},cfg.labels.micro{ilabel}))
                %                             micro_filenrs       = [micro_filenrs, ifile];
                %                             microlabel{ifile}   = cfg.labels.micro{ilabel};
                %                         end
                %                     end
                %                 end
                
                if write
                    clear fname
                    temp                            = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
                    fname{1}                        = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                    cfgtemp                         = [];
                    cfgtemp.dataset                 = fname{1};
                    fprintf('LOADING: %s\n',cfgtemp.dataset);
                    clear fname
                    fname{1}                        = cfgtemp.dataset;
                    dirdat{idir}                    = ft_read_neuralynx_interp(fname);
                    
                    if strcmp(cfg.circus.reref, 'yes')
                        temp                        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.refchan,'.ncs']));
                        cfgtemp                     = [];
                        cfgtemp.dataset             = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, temp.name);
                        fprintf('LOADING (reference): %s\n',cfgtemp.dataset);
                        clear fname
                        fname{1}                    = cfgtemp.dataset;
                        refdat                      = ft_read_neuralynx_interp(fname);
                        dirdat{idir}.trial{1}       = dirdat{idir}.trial{1} - refdat.trial{1};
                        clear refdat
                    end
                    
                    % Large offsets can create artefacts, dealt with by
                    % filtering. Should not happen with continuous data
                    if strcmp(cfg.circus.hpfilter, 'yes')
                        cfgtemp                   = [];
                        cfgtemp.hpfilter          = cfg.circus.hpfilter;
                        cfgtemp.hpfreq            = cfg.circus.hpfreq;
                        dirdat{idir}              = ft_preprocessing(cfgtemp,dirdat{idir});
                    end
                    
                    % truncate label to make them equal over files
                    dirdat{idir}.label{1} = dirdat{idir}.label{1}(end-6:end); % can be replaced by circus.channel
                    
                end % write
                
                % create sampleinfo (nr of samples in file)
                temp = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.channel{ichan},'.ncs']));
                hdr  = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, temp.name));
                
                % save sampleinfo to reconstruct data again after reading SC
                sampleinfo{ipart}{ichan}(idir,:) = [1 hdr.nSamples];
                
            end
            
            % concatinate data over files
            if write
                chandat = dirdat{1};
                for idir = 2 : length(cfg.directorylist{ipart})
                    fprintf('Concatinating directory %d, channel %d\n', idir, ichan);
                    chandat.trial{1}        = [chandat.trial{1} dirdat{idir}.trial{1}];
                    chandat.time{1}         = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
                end
            end
            
            % create filename for concatinated data
            temp        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
            hdrtemp     = ft_read_header(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name));
            clear fname
            subjdir     = cfg.prefix(1:end-1);
            partdir     = ['p', num2str(ipart)];
            filename    = [cfg.prefix, 'p', num2str(ipart), '-multifile-', cfg.circus.channel{ichan}, '.ncs'];
            fname       = fullfile(cfg.datasavedir, subjdir, partdir, filename);
            
            % save filenames to cfg (output of function)
            cfg.fnames_ncs{ipart}{ichan} = fname;
            
            % write data in .ncs format
            if write
                
                % if it doesnt exist, create subject and part directory
                if ~exist(fullfile(cfg.datasavedir, subjdir),'dir')
                    mkdir(fullfile(cfg.datasavedir, subjdir));
                end
                if ~exist(fullfile(cfg.datasavedir, subjdir, partdir),'dir')
                    mkdir(fullfile(cfg.datasavedir, subjdir, partdir));
                end
                
                hdr                         = [];
                hdr.Fs                      = hdrtemp.Fs;
                hdr.nSamples                = size(chandat.trial{1},2);
                hdr.nSamplePre              = 0;
                hdr.nChans                  = 1;
                hdr.FirstTimeStamp          = 0;
                hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
                hdr.label                   = chandat.label;
                ft_write_data(fname,chandat.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
                
                %
                %                 hdr                         = [];
                %                 hdr.Fs                      = hdrtemp.Fs;
                %                 hdr.nSamples                = size(dirdat{1}.trial{1},2);
                %                 hdr.nSamplePre              = 0;
                %                 hdr.nChans                  = 1;
                %                 hdr.FirstTimeStamp          = 0;
                %                 hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
                %                 hdr.label                   = dirdat{1}.label;
                %                 ft_write_data(fname{1},dirdat{1}.trial{1},'chanindx',1,'dataformat','neuralynx_ncs','header',hdr);
                
                %
                %                 st_FieldSelection  = ones(1,6);
                %
                %                 s_ExtractHeader    = 1;
                %                 s_ExtractMode      = 1;
                %                 v_ModeArray        = [];
                %                 s_AppendToFileFlag = 0; %1 if append data to file when creating ncs file
                %
                %                 s_ExportMode       = 1; %1 if a
                %                 v_ExportModeVector = [];
                %                 if ispc
                %                     Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, st_FieldSelection, ...
                %                         v_TimestampI, v_ChanNumI, v_SampleFrequency, v_NumValSamplesI, v_SamplesI, st_Header);
                %                 end
                %                 if isunix
                %                     v_NumRecs = length(v_TimestampI);
                %                     Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, v_NumRecs, st_FieldSelection, ...
                %                         v_TimestampI, v_ChanNumI, v_SampleFrequency, v_NumValSamplesI, v_SamplesI, st_Header);
                %                 end
            end
            
            clear chandat
            clear dirdat
            
        end % ichan
        
        if ~write
            continue
        end
        
        if ~strcmp(cfg.circus.writedeadfile, 'yes')
            continue
        end
        
        % write deadtime, i.e. artefact file for Spyking-Circus
        deadfile_ms         = [];
        deadfile_samples    = [];
        last_ms             = 0;
        last_samples        = 0;
        dirlist{ipart}      = [];
        
        for idir = 1 : size(cfg.directorylist{ipart},2)
            if isfield(MuseStruct{ipart}{idir},'markers')
                if isfield(MuseStruct{ipart}{idir}.markers,'BAD__START__')
                    if isfield(MuseStruct{ipart}{idir}.markers.BAD__START__,'synctime')
                        
                        % check if there is an equal amount of start and end markers
                        if size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) == 0
                            
                            temp    = [];
                            hdr     = [];
                            temp    = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.channel{1},'.ncs'])); %first chan. all chans have the same nr of samples
                            hdr     = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, temp.name));
                            
                            fprintf('Great, recovered same number of start and end markers \n')
                            deadfile_ms         = [deadfile_ms;      MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                            deadfile_samples    = [deadfile_samples; MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples, MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
                            last_samples        = last_samples  + hdr.nSamples;
                            last_ms             = last_ms       + hdr.nSamples/hdr.Fs * 1000;
                            
                        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2) - size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) > 0
                            
                            fprintf('ERROR! more start than end found in %s \n',MuseStruct{ipart}{idir}.directory);
                            
                        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime,2) - size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime,2) < 0
                            
                            fprintf('ERROR! more end than start found in %s - CORRECTING \n',MuseStruct{ipart}{idir}.directory)
                            for i = 1 : 10
                                for itrial = 1 : length(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime)
                                    start(itrial) = MuseStruct{ipart}{idir}.markers.BAD__START__.synctime(itrial);
                                    stop(itrial)  = MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(itrial);
                                end
                                
                                x = find(start > stop, 1, 'first');
                                if ~isempty(x)
                                    MuseStruct{ipart}{idir}.markers.BAD__END__.synctime(x) = [];
                                end
                            end
                            
                            deadfile_ms         = [deadfile_ms;      MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,        MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
                            deadfile_samples    = [deadfile_samples; MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples, MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
                            last_samples        = last_samples + hdr.nSamples;
                            last_ms             = last_ms + hdr.nSamples/hdr.Fs * 1000;
                        end
                    else
                        fprintf('Found no artefacts for file: %s\n', MuseStruct{ipart}{idir}.directory);
                    end
                else
                    fprintf('Found no artefacts for file: %s\n', MuseStruct{ipart}{idir}.directory);
                end
            end
            dirlist{ipart} = [dirlist{ipart}; MuseStruct{ipart}{idir}.directory];
            fprintf('%d\n',idir);
        end
        
        % return info
        cfg.deadfile_ms{ipart}         = deadfile_ms;
        cfg.deadfile_samples{ipart}    = deadfile_samples;
        
        % write artefacts to txt file for spyking circus
        subjdir  = cfg.prefix(1:end-1);
        partdir  = ['p', num2str(ipart)];
        
        filename = 'SpykingCircus_artefacts_ms.dead';
        fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
        dlmwrite(fullfile(cfg.datasavedir,subjdir,partdir,filename),deadfile_ms,'delimiter','	','precision','%.4f');
        
        filename = 'SpykingCircus_artefacts_samples.dead';
        fprintf('Writing artefacts for Spyking-Circus to: %s\n',filename);
        dlmwrite(fullfile(cfg.datasavedir,subjdir,partdir,filename),deadfile_samples,'delimiter','	','precision','%.4f');
        
        filename = 'SpykingCircus_dirlist.txt';
        fprintf('Writing list of directories for Spyking-Circus to: %s\n',filename);
        writematrix(dirlist{ipart},fullfile(cfg.datasavedir,subjdir,partdir,filename));
        
    end % ipart
    
    % save trialinfo stuff
    save(fname_output,'cfg');
end
