function [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(cfg, force)

% writeSpykingCircusFileList writes a textfile with filenames for
% SpykingCircus and returns sampleinfo for later (cutting results back up)
%
% use as
%   [filelist, sampleinfo, timestamps, hdr] = writeSpykingCircusFileList(cfg, force)
%
% "force" is to force recalculation

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


% TODO: deal with several channel combination in 'channelname'
cfg.circus.channelname      = ft_getopt(cfg.circus, 'channelname', []);
cfg.circus.timestamps       = ft_getopt(cfg.circus, 'timestamps', false);
fname_output                = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_filelist.mat']);

if nargin == 1
    if exist(fname_output, 'file')
        fprintf('Reading %s\n', fname_output);
        % repeat to deal with load errors
        count = 0;
        err_count = 0;
        while count == err_count
            try
                if cfg.circus.timestamps
                    load(fname_output, 'filelist', 'sampleinfo', 'timestamps', 'hdr');
                else
                    load(fname_output, 'filelist', 'sampleinfo', 'hdr');
                    timestamps = [];
                end
            catch ME
                err_count = err_count + 1;
                disp('Something went wrong loading the file. Trying again...')
            end
            count = count + 1;
        end
        return;
    else
        warning('No precomputed data is found, not enough input arguments to compute data');
        return
    end
end

if exist(fname_output, 'file') && force == false
    fprintf('Loading %s\n', fname_output);
    
    count = 0;
    err_count = 0;
    while count == err_count
        try
            if cfg.circus.timestamps
                load(fname_output, 'filelist', 'sampleinfo', 'timestamps', 'hdr');
            else
                load(fname_output, 'filelist', 'sampleinfo', 'hdr');
                timestamps = [];
            end
        catch ME
            err_count = err_count + 1;
            disp('Something went wrong loading the file. Trying again...')
        end
        count = count + 1;
    end
    return
end

filelist = [];

% loop through different parts
for ipart = 1 : size(cfg.directorylist, 2)
    
    fprintf('\n*** Starting on part %d ***\n', ipart)
    
    if isempty(cfg.circus.channelname)
        
        % process channels separately
        for ichan = 1 : size(cfg.circus.channel, 2 )
            
            % read one channel for all directories (time)
            for idir = 1 : size(cfg.directorylist{ipart},2)
                
                temp                            = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan},'.ncs']));
                fname                           = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                fname_filelist = fname;
                if ispc %correct file name for its use with linux
                    fname_filelist = strrep(fname_filelist, '\\lexport\iss01.', '/network/lustre/iss01/');
                    fname_filelist = strrep(fname_filelist, '\\l2export\iss02.', '/network/lustre/iss02/');
                    fname_filelist = strrep(fname_filelist, '\', '/');
                end
                filelist{ipart}(idir, ichan)    = string(fname_filelist);
                
                if ichan == 1
                    fprintf('Reading header & timestamps %s\n', fname)
                    hdr{ipart}{idir}            = ft_read_header(fname);
                    sampleinfo{ipart}(idir, 1)  = 1;
                    sampleinfo{ipart}(idir, 2)  = hdr{ipart}{idir}.nSamples;
                    if cfg.circus.timestamps
                        timestamps{ipart}{idir} = ft_read_data(fname, 'timestamp', 'true');  % take the first concatinated channel to extract the timestamps
                    else
                        timestamps{ipart}{idir} = [];
                    end
                end
            end
        end % ichan
        
        % save filelist
        if ~isempty(cfg.circus.channelname); chandir = cfg.circus.channelname{ichan}; else; chandir = []; end
        
        subjdir     = cfg.prefix(1:end-1);
        partdir     = ['p', num2str(ipart)];
        fname       = fullfile(cfg.datasavedir, subjdir, partdir, chandir, 'filelist.txt');
       
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
        fid = fopen(fname, 'w');
        for irow = 1 : size(filelist{ipart}, 1)
            for icol = 1 : size(filelist{ipart}, 2)
                fprintf(fid, '%s\t', filelist{ipart}(irow, icol));
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
        
    else
        
        % go through different channels
        for chandir = string(unique(cfg.circus.channelname))
            
            chanlist = find(contains(cfg.circus.channel, chandir));
            %             disp(chanlist);
            % process channels separately
            
            for ichan = 1 : length(chanlist)
                
                % read one channel for all directories (time)
                for idir = 1 : size(cfg.directorylist{ipart},2)
                    
                    temp                                    = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{chanlist(ichan)},'.ncs']));
                    fname                                   = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                    fname_filelist = fname;
                    if ispc %correct file name for its use with linux
                        fname_filelist = strrep(fname_filelist, '\\lexport\iss01.', '/network/lustre/iss01/');
                        fname_filelist = strrep(fname_filelist, '\\l2export\iss02.', '/network/lustre/iss02/');
                        fname_filelist = strrep(fname_filelist, '\', '/');
                    end
                    filelist{ipart}.(chandir)(idir, ichan) = string(fname_filelist);
                    
                    if ichan == 1
                        fprintf('Reading header & timestamps %s\n', fname)
                        hdr{ipart}{idir}            = ft_read_header(fname);
                        sampleinfo{ipart}(idir, 1)  = 1;
                        sampleinfo{ipart}(idir, 2)  = hdr{ipart}{idir}.nSamples;
                        if cfg.circus.timestamps
                            timestamps{ipart}{idir} = ft_read_data(fname, 'timestamp', 'true');  % take the first concatinated channel to extract the timestamps
                        else
                            timestamps{ipart}{idir} = [];
                        end
                    end
                end
            end % ichan
            
            % save filelist
            subjdir     = cfg.prefix(1:end-1);
            partdir     = ['p', num2str(ipart)];
            fname       = fullfile(cfg.datasavedir, subjdir, partdir, chandir, 'filelist.txt');
            
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
            fid = fopen(fname, 'w');
            for irow = 1 : size(filelist{ipart}.(chandir), 1)
                for icol = 1 : size(filelist{ipart}.(chandir), 2)
                    fprintf(fid, '%s\t', filelist{ipart}.(chandir)(irow, icol));
                end
                fprintf(fid, '\n');
            end
            fclose(fid);
            
        end
    end
end % ipart

% save trialinfo stuff
if cfg.circus.timestamps
    save(fname_output,'filelist', 'sampleinfo', 'timestamps', 'hdr', '-v7.3');
else
    save(fname_output,'filelist', 'sampleinfo', 'hdr', '-v7.3');
end