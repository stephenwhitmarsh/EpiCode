function [artefacts] = writeSpykingCircusDeadfiles(cfg, MuseStruct, force, suffix)

% WRITESPYKINGCIRCUS_DEADFILES creates artefact files (deadtime) for
% SpykingCircus in data directories for spike analysis
%
% use as
%   writeSpykingCircus_deadfiles(cfg, MuseStruct, force)
% or 
%   writeSpykingCircus_deadfiles(cfg, MuseStruct, force, suffix)
% 
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

if nargin < 4
    suffix = [];
end
fname_output           = fullfile(cfg.datasavedir, [cfg.prefix, 'SpykingCircus_artefacts', char(suffix), '.mat']);
cfg.circus.channelname = ft_getopt(cfg.circus, 'channelname', []);

if exist(fname_output, 'file') && force == false
    fprintf('\nLoading sampleinfo artefacts: %s \n', fname_output);
    load(fname_output, 'artefacts');
    return
end

% loop through different parts
for ipart = 1 : size(cfg.directorylist, 2)
    
    fprintf('\n*** Writing artefacts of part %d ***\n', ipart)
    
    % write deadtime, i.e. artefact file for Spyking-Circus
    deadfile_ms         = [];
    deadfile_samples    = [];
    deadfile_idir       = [];
    last_ms             = 0;
    last_samples        = 0;
    dirlist{ipart}      = [];
    
    for idir = 1 : size(cfg.directorylist{ipart}, 2)
        temp                            = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{1}, '.ncs']));
        hdr                             = ft_read_header(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name));
        
        if idir > 1
            last_samples                    = last_samples  + hdr.nSamples;
            last_ms                         = last_ms       + hdr.nSamples/hdr.Fs * 1000;
        end
        
        if ~isfield(MuseStruct{ipart}{idir}, 'markers')
            continue
        end
        if ~isfield(MuseStruct{ipart}{idir}.markers, 'BAD__START__')
            fprintf('Found no artefacts for file: %s\n', MuseStruct{ipart}{idir}.directory);
            continue
        end
        if ~isfield(MuseStruct{ipart}{idir}.markers.BAD__START__, 'synctime')
            fprintf('Found no artefacts for file: %s\n', MuseStruct{ipart}{idir}.directory);
            continue
        end
        
        % check if there is an equal amount of start and end markers
        if size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime, 2)-size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime, 2) == 0
            fprintf('Great, recovered same number of start and end markers \n')
            
            artefacts.ms{ipart}{idir}       = [MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,         MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
            artefacts.samples{ipart}{idir}  = [MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples,  MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
            deadfile_ms                     = [deadfile_ms; artefacts.ms{ipart}{idir}];
            deadfile_samples                = [deadfile_samples; artefacts.samples{ipart}{idir} ];
            deadfile_idir                   = [deadfile_idir; ones(size(artefacts.ms{ipart}{idir}))*idir];
            
        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime, 2) - size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime, 2) > 0
            fprintf('ERROR! more start than end found in %s \n', MuseStruct{ipart}{idir}.directory);
            
        elseif size(MuseStruct{ipart}{idir}.markers.BAD__START__.synctime, 2) - size(MuseStruct{ipart}{idir}.markers.BAD__END__.synctime, 2) < 0
            fprintf('ERROR! more end than start found in %s - CORRECTING \n', MuseStruct{ipart}{idir}.directory)
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
            artefacts.ms{ipart}{idir}       = [MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*1000+last_ms,         MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*1000+last_ms];
            artefacts.samples{ipart}{idir}  = [MuseStruct{ipart}{idir}.markers.BAD__START__.synctime'*hdr.Fs+last_samples,  MuseStruct{ipart}{idir}.markers.BAD__END__.synctime'*hdr.Fs+last_samples];
            deadfile_ms                     = [deadfile_ms; artefacts.ms{ipart}{idir}];
            deadfile_samples                = [deadfile_samples; artefacts.samples{ipart}{idir} ];       
            deadfile_idir                   = [deadfile_idir; ones(size(artefacts.ms{ipart}{idir}))*idir];
        end
        dirlist{ipart} = [dirlist{ipart}; MuseStruct{ipart}{idir}.directory];
        fprintf('%d\n', idir);
    end
    
    % write artefacts to txt file for spyking circus
    subjdir  = cfg.prefix(1:end-1);
    partdir  = ['p', num2str(ipart)];
    
    % if it doesnt exist, create directory
    if ~exist(fullfile(cfg.datasavedir, subjdir), 'dir')
        sprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir));
        mkdir(fullfile(cfg.datasavedir, subjdir));
    end
    
    if ~exist(fullfile(cfg.datasavedir, subjdir, partdir), 'dir')
        sprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir, partdir));
        mkdir(fullfile(cfg.datasavedir, subjdir, partdir));
    end
    
    if isempty(cfg.circus.channelname)
        filename = sprintf('SpykingCircus_artefacts_ms%s.dead', suffix);
        fprintf('Writing artefacts for Spyking-Circus to: %s\n', filename);
        dlmwrite(fullfile(cfg.datasavedir, subjdir, partdir, filename), deadfile_ms, 'delimiter', '	', 'precision', '%.4f');
        
        filename = sprintf('SpykingCircus_artefacts_samples%s.dead', suffix);
        fprintf('Writing artefacts for Spyking-Circus to: %s\n', filename);
        dlmwrite(fullfile(cfg.datasavedir, subjdir, partdir, filename), deadfile_samples, 'delimiter', '	', 'precision', '%.4f');    
    else
        for chandir = unique(cfg.circus.channelname)
            
            % if it doesnt exist, create directory
            if ~exist(fullfile(cfg.datasavedir, subjdir, partdir, string(chandir)), 'dir')
                sprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir, partdir, string(chandir)));
                mkdir(fullfile(cfg.datasavedir, subjdir, partdir, string(chandir)));
            end

            filename = sprintf('SpykingCircus_artefacts_ms%s.dead', suffix);
            fprintf('Writing artefacts for Spyking-Circus to: %s\n', filename);
            dlmwrite(fullfile(cfg.datasavedir, subjdir, partdir, string(chandir), filename), deadfile_ms, 'delimiter', '	', 'precision', '%.4f');
            
            filename = sprintf('SpykingCircus_artefacts_samples%s.dead', suffix);
            fprintf('Writing artefacts for Spyking-Circus to: %s\n', filename);
            dlmwrite(fullfile(cfg.datasavedir, subjdir, partdir, string(chandir), filename), deadfile_samples, 'delimiter', '	', 'precision', '%.4f');
        end
    end
    
    % return info
    artefacts.deadfile_ms{ipart}         = deadfile_ms;
    artefacts.deadfile_samples{ipart}    = deadfile_samples;
    artefacts.idir{ipart}                = deadfile_idir;
    save(fname_output, 'artefacts');

end