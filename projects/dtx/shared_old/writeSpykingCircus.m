function [cfg] = writeSpykingCircus(cfg, force)

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

cfg.circus.writedeadfile    = ft_getopt(cfg.circus, 'writedeadfile', 'yes');
cfg.circus.hpfilter         = ft_getopt(cfg.circus, 'hpfilter', 'no');
cfg.circus.hpfreq           = ft_getopt(cfg.circus, 'hpfilter', 'no');
cfg.circus.channelname      = ft_getopt(cfg.circus, 'channelname', []);
cfg.circus.version          = ft_getopt(cfg.circus, 'version', 'mex');
fname_output                = fullfile(cfg.datasavedir,[cfg.prefix,'SpykingCircus_trialinfo_parts.mat']);

if exist(fname_output, 'file') && force == false
    fprintf('\n Loading sampleinfo: %s \n', fname_output);
    load(fname_output, 'cfg');
else

    % loop through different parts
    for ipart = 1 : size(cfg.directorylist, 2)

        fprintf('\n*** Starting on part %d ***\n',ipart)

        % process channels separately
        for ichan = 1 : size(cfg.circus.channel,2)

            cfg.sampleinfo{ipart}{ichan} = [];
            clear chandat dirdat v_Timestamp

            % read one channel for all directories (time)
            for idir = 1 : size(cfg.directorylist{ipart},2)

                % get filename to read
                clear fname
                temp     = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
                fname{1} = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);

                % use updated fieldtrip MATLAB function
                if strcmp(cfg.circus.version, 'fieldtrip')

                    fprintf('Loading %s using Fieldtrip\n', fname{1});
                    fprintf('***************************************');
                    fprintf('*** CHECK POLARITY AS IT IS IGNORED ***');
                    fprintf('***************************************');

                    cfgtemp                         = [];
                    dirdat{idir}                    = ft_read_neuralynx_interp(fname);

                    if strcmp(cfg.circus.reref, 'yes')
                        temp                        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.refchan,'.ncs']));
                        fname_ref{1}                = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                        %fprintf('LOADING (reference): %s\n',cfgtemp.dataset);
                        refdat                      = ft_read_neuralynx_interp(fname_ref);
                        dirdat{idir}.trial{1}       = dirdat{idir}.trial{1} - refdat.trial{1};
                        clear refdat
                    end

                % use NeuraLynx MEX file
                elseif strcmp(cfg.circus.version, 'mex')

                    fprintf('Loading %s using NeuraLynx MEX\n', fname{1});
                    fprintf('****************************************************');
                    fprintf('*** CHECK POLARITY AS IT IS CORRECTED IF FLIPPED ***');
                    fprintf('****************************************************');

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

                    if ispc
                        [v_Timestamp{idir}, v_ChanNum, v_SampleFrequency, v_NumValSamples, dirdat{idir}.trial{1}, st_Header] = Nlx2MatCSC(fname{1}, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);
                    else
                        [v_Timestamp{idir}, v_ChanNum, v_SampleFrequency, v_NumValSamples, dirdat{idir}.trial{1}, st_Header] = Nlx2MatCSC_v3(fname{1}, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);
                    end

                    if strcmp(cfg.circus.reref, 'yes')
                        temp                        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.refchan,'.ncs']));
                        fname_ref{1}                = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                        fprintf('LOADING (reference): %s\n',cfgtemp.dataset);
                        if ispc
                            [~, ~, ~, ~, refdat.trial{1}, ~] = Nlx2MatCSC(fname{1}, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);
                        else
                            [~, ~, ~, ~, refdat.trial{1}, ~] = Nlx2MatCSC_v3(fname{1}, st_FieldSelection(1:5), s_ExtractHeader, s_ExtractMode, v_ModeArray);
                        end
                        dirdat{idir}.trial{1}       = dirdat{idir}.trial{1} - refdat.trial{1};
                        clear refdat
                    end

                else
                    error('Version not supported');
                end

                % Large offsets can create artefacts, dealt with by
                % filtering. Should not happen with continuous data
                if strcmp(cfg.circus.hpfilter, 'yes')
                    cfgtemp                   = [];
                    cfgtemp.hpfilter          = cfg.circus.hpfilter;
                    cfgtemp.hpfreq            = cfg.circus.hpfreq;
                    dirdat{idir}              = ft_preprocessing(cfgtemp, dirdat{idir});
                end

                % truncate label to make them the same over subsequent
                % directories
                dirdat{idir}.label{1} = cfg.circus.channel{ichan};

                % create sampleinfo (nr of samples in file)
                temp = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.channel{ichan},'.ncs']));
                hdr  = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir}, temp.name));

                % save sampleinfo to reconstruct data again after reading SC
                cfg.sampleinfo{ipart}{ichan}(idir,:) = [1 hdr.nSamples];

            end

            % concatinate data over files
            if strcmp(cfg.circus.version, 'fieldtrip')
                chandat = dirdat{1};
                for idir = 2 : length(cfg.directorylist{ipart})
                    fprintf('Concatinating directory %d, channel %d\n', idir, ichan);
                    chandat.trial{1}    = [chandat.trial{1} dirdat{idir}.trial{1}];
                    chandat.time{1}     = [chandat.time{1} (dirdat{idir}.time{1} + chandat.time{1}(end))];
                end

            elseif strcmp(cfg.circus.version, 'mex')
                chandat             = dirdat{1}.trial{1};
                timestamps_concat   = v_Timestamp{1};
                for idir = 2 : length(cfg.directorylist{ipart})
                    fprintf('Concatinating directory %d, channel %d\n', idir, ichan);
                    chandat             = [chandat dirdat{idir}.trial{1}];
                    timestamps_concat   = [timestamps_concat v_Timestamp{idir}];
                end
            end

            % create filename for concatinated data
            if ~isempty(cfg.circus.channelname); chandir = cfg.circus.channelname{ichan}; else; chandir = []; end
            temp        = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.circus.channel{ichan}, '.ncs']));
            hdrtemp     = ft_read_header(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name));
            filename    = [cfg.prefix, 'p', num2str(ipart), '-multifile-', cfg.circus.channel{ichan}, '.ncs'];
            subjdir     = cfg.prefix(1:end-1);
            partdir     = ['p', num2str(ipart)];
            fname       = fullfile(cfg.datasavedir, subjdir, partdir, chandir, filename);

            % save filenames to cfg (output of function)
            cfg.fnames_ncs{ipart}{ichan} = fname;

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

            if strcmp(cfg.circus.version, 'fieldtrip')

                hdr                         = [];
                hdr.Fs                      = hdrtemp.Fs;
                hdr.nSamples                = size(chandat.trial{1},2);
                hdr.nSamplePre              = 0;
                hdr.nChans                  = 1;
                hdr.FirstTimeStamp          = 0;
                hdr.TimeStampPerSample      = hdrtemp.TimeStampPerSample;
                hdr.label                   = chandat.label;
                sprintf('Writing concatinated data with MATLAB to %s\n', fname);
                delete(fname);
                ft_write_data(fname, chandat.trial{1}, 'chanindx', 1, 'dataformat', 'neuralynx_ncs', 'header', hdr);

            elseif strcmp(cfg.circus.version, 'mex')

                if ispc
                    fprintf('Writing concatinated data with MEX file to %s\n', fname);
                    delete(fname);
                    Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, st_FieldSelection, timestamps_concat, ...
                        ones(size(timestamps_concat)) * v_ChanNum(1), ...
                        ones(size(timestamps_concat)) * v_SampleFrequency(1), ...
                        ones(size(timestamps_concat)) * v_NumValSamples(1), ...
                        chandat, st_Header);
                end
                if isunix % this write a file that is too long somehow - second half without data, too many timestamps?
                    delete(fname);
                    v_NumRecs = length(timestamps_concat);
                    fprintf('Writing concatinated data with MEX file to %s\n', fname);
                    Mat2NlxCSC(fname, s_AppendToFileFlag, s_ExportMode, v_ExportModeVector, v_NumRecs, st_FieldSelection, timestamps_concat, ...
                        ones(size(timestamps_concat)) * v_ChanNum(1), ...
                        ones(size(timestamps_concat)) * v_SampleFrequency(1), ...
                        ones(size(timestamps_concat)) * v_NumValSamples(1), ...
                        chandat, st_Header);
                end
            else
                error('Version not supported');
            end

            clear chandat dirdat

        end % ichan
    end % ipart

    % save trialinfo stuff
    save(fname_output,'cfg');
end
