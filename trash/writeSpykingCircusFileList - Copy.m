function [cfg] = writeSpykingCircusFileList(cfg, force)

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

cfg.circus.channelname      = ft_getopt(cfg.circus, 'channelname', []);
fname_output                = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_filelist.mat']);

if exist(fname_output, 'file') && force == false
    fprintf('\n Loading sampleinfo: %s \n', fname_output);
    load(fname_output, 'cfg');
    continue
end

filelist = [];

% loop through different parts
for ipart = 1 : size(cfg.directorylist, 2)
    
    fprintf('\n*** Starting on part %d ***\n', ipart)
    
    % process channels separately
    for ichan = 1 : size(cfg.circus.channel, 2 )
        
        sampleinfo{ipart}{ichan} = [];
        clear chandat dirdat v_Timestamp
        
        % read one channel for all directories (time)
        for idir = 1 : size(cfg.directorylist{ipart},2)
       
            % get sampleinfo (nr of samples in file)
            temp = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan},'.ncs']));
            fname = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
            hdr  = ft_read_header(fname);
            
            % save sampleinfo to reconstruct data again after reading SC
            sampleinfo{ipart}{ichan}(idir,:) = [1, hdr.nSamples];
            
            filelist{ipart}(idir, ichan) = fname;         
        end
        
        % save filelist
        if ~isempty(cfg.circus.channelname); chandir = cfg.circus.channelname{ichan}; else; chandir = []; end
        temp        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
        hdrtemp     = ft_read_header(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name));
        temp        = [cfg.prefix, 'p', chandir, 'filelist.txt');
        subjdir     = cfg.prefix(1:end-1);
        partdir     = ['p', num2str(ipart)];
        fname       = fullfile(cfg.datasavedir, subjdir, partdir, chandir, temp);
        
        % if it doesnt exist, create directory
        if ~exist(fullfile(cfg.datasavedir, subjdir), 'dir')
            fprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir));
            mkdir(fullfile(cfg.datasavedir, subjdir));
        end
        
        if ~exist(fullfile(cfg.datasavedir, subjdir, partdir), 'dir')
            fprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir, partdir));
            mkdir(fullfile(cfg.datasavedir, subjdir, partdir));
        end
        
        if ~exist(fullfile(cfg.datasavedir, subjdir, partdir, chandir), 'dir')
            fprintf('Creating directory %s', fullfile(cfg.datasavedir, subjdir, partdir, chandir));
            mkdir(fullfile(cfg.datasavedir, subjdir, partdir, chandir));
        end
        
        tdlwrite(fname, filelist);

    end % ichan
end % ipart

% save trialinfo stuff
save(fname_output,'cfg');
