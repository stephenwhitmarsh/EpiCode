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
                
                if write
                    clear fname
                    temp                            = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
                    fname                           = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                    fprintf('LOADING: %s\n', fname);
                    
                    st_FieldSelection(1) = 1; %timestamps
                    st_FieldSelection(2) = 1; %Channel Numbers
                    st_FieldSelection(3) = 1; %sample freq
                    st_FieldSelection(4) = 1; %Number of Valid Samples
                    st_FieldSelection(5) = 1; %samples
                    st_FieldSelection(6) = 1; %header (for creating .ncs only)
                    s_ExtractHeader      = 1;
                    s_ExtractMode        = 1; %1 if all samples
                    v_ModeArray          = []; %[] if all, 2 elements if range
                    s_AppendToFileFlag   = 0;
                    s_ExportMode         = 1; %1 if a
                    v_ExportModeVector   = [];
                    
                    [v_Timestamp{idir}, v_ChanNum, v_SampleFrequency, v_NumValSamples, dirdat{idir}, st_Header] = Nlx2MatCSC_v3(fname, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);
                    
                    if strcmp(cfg.circus.reref, 'yes')
                        temp                        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*',cfg.circus.refchan, '.ncs']));
                        fnameref                    = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                        fprintf('LOADING (reference): %s\n', fnameref);
                        [~, ~, ~, ~, refdat, ~]     = Nlx2MatCSC_v3(fname, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);                          
                        dirdat{idir}                = dirdat{idir} - refdat;
                        clear refdat
                    end
                    
                end % write
              
                % save sampleinfo to reconstruct data again after reading SC
                sampleinfo{ipart}{ichan}(idir,:)            = [1 numel(dirdat{idir})];
            end
            
            % concatinate data over files
            if write
                chandat = dirdat{1};
                timestamps_concat = v_Timestamp{1};
                for idir = 2 %length(cfg.directorylist{ipart}) 8?
                    fprintf('Concatinating directory %d, channel %d\n', idir, ichan);
                    chandat = [chandat dirdat{idir}];
                    timestamps_concat = [timestamps_concat v_Timestamp{idir}];
                end
            end
            
            % create filename for concatinated data
            temp        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
            subjdir     = cfg.prefix(1:end-1);
            partdir     = ['p', num2str(ipart)];
            filename    = [cfg.prefix, 'p', num2str(ipart), '-multifile-', num2str(ichan), '.ncs'];
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
  
                if ispc
                    Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, st_FieldSelection, timestamps_concat, v_ChanNum, v_SampleFrequency, v_NumValSamples, chandat, st_Header);
                end
                if isunix
                    v_NumRecs = length(timestamps_concat);
                   
                    Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, v_NumRecs, st_FieldSelection, timestamps_concat, ...
                        ones(size(timestamps_concat)) * v_ChanNum(1), ...
                        ones(size(timestamps_concat)) * v_SampleFrequency(1), ...
                        ones(size(timestamps_concat)) * v_NumValSamples(1), ...
                        chandat, st_Header);
                end
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
