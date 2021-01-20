function [SpikeRaw, SpikeTrials] = readSpykingCircusPSG(cfg,MuseStruct,force,varargin)

% READSPYKINGCIRCUSPSG reads Spyking-Circus results: spike-times and templates.
%
% use as
%   [SpikeRaw, SpikeTrials] = readSpykingCircusPSG(cfg,MuseStruct,force,varargin)
%
% Necessary input:
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files
% cfg.circus.channel    = micro electrode names
% MuseStruct{ipart}     = info (e.g. events, files) of original data,
%                         used to segment the spikes into trials
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% Output:
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% SpikeTrials           = spike data epoched in FieldTrip trial data structure

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

if isempty(varargin) || strcmp(varargin{1},'all')
    parts_to_read = 1:size(cfg.directorylist,2);
else
    parts_to_read = varargin{1};
end


fname = fullfile(cfg.datasavedir,[cfg.prefix,'spikedataPSG.mat']);
if exist(fname,'file') && force == false
    load(fname,'SpikeRaw','SpikeTrials');
else

    for ipart = parts_to_read

        % find spiking-circus output path, which is based on the name of the
        % first datafile
        temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'_*.result',cfg.circus.postfix,'.hdf5']));
        if isempty(temp)
            fprintf('Could not find Spyking-Circus results: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'_*.result',cfg.circus.postfix,'.hdf5']));
            return
        else
            fname_spikes = fullfile(temp.folder,temp.name);
        end

        temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'_*.templates',cfg.circus.postfix,'.hdf5']));
        if isempty(temp)
            fprintf('Could not find Spyking-Circus templates: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'_*.templates',cfg.circus.postfix,'.hdf5']));
            return
        end
        fname_templates = fullfile(temp.folder,temp.name);

        if ~exist(fname_spikes,'file')
            disp('Can''t find spyking circus results.');
            break
        end

        % load spiking data
        datinfo                     = h5info(fname_spikes);
        temp                        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p' num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'_*.ncs']));
        hdr_fname                   = fullfile(temp(1).folder,temp(2).name);
        hdr                         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
        SpikeRaw{ipart}.trialtime   = [0 hdr.nSamples / hdr.Fs];

        fprintf('Loading timestamps from: %s. This can take a while...',hdr_fname);
        timestamps                  = ft_read_data(hdr_fname,'timestamp','true'); fprintf('Done\n'); % take the first concatinated file to extract the timestamps
%       timestamps  =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster
%       [timestamps, ~, ~] = Nlx2MatCSC_v3(hdr_fname, [0 0 0 1 0]);

        % read spiketimes of clusters
        fprintf('Loading spike data from: %s\n',hdr_fname);
        for i = 1 : size(datinfo.Groups,1)
            names(i) = string(datinfo.Groups(i).Name);
            if strfind(names(i),'spiketimes')
                spiketimes_indx = i;
            end
        end

        clear clusternr
        for i = 1 : size(datinfo.Groups(spiketimes_indx).Datasets,1) % number of templates
            SpikeRaw{ipart}.label{i} = datinfo.Groups(spiketimes_indx).Datasets(i).Name;
            temp = strsplit(SpikeRaw{ipart}.label{i},'_');
            clusternr(i) = str2double(temp{2});
        end

        for i = 1:numel(clusternr)

            disp(['Reading cluster number ',num2str(i)]);

            % read spike timings (in seconds)
            datasetname = char(strcat('/spiketimes/',SpikeRaw{ipart}.label{i}));
            SpikeRaw{ipart}.sample{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0

            % read amplitudes
            datasetname = char(strcat('/amplitudes/',SpikeRaw{ipart}.label{i}));
            SpikeRaw{ipart}.amplitude{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0

            % map samplenrs onto timestamps
            SpikeRaw{ipart}.timestamp{clusternr(i)+1} = timestamps(SpikeRaw{ipart}.sample{clusternr(i)+1});
            % SpikeRaw.timestamp{clusternr(i)+1} = int64(SpikeRaw.sample{clusternr(i)+1}) * int64(hdr.TimeStampPerSample) + int64(hdr.FirstTimeStamp);
        end

        % load templates
        templates_size  = double(h5read(fname_templates, '/temp_shape'));
        N_e             = templates_size(2);
        N_t             = templates_size(1);
        temp_x          = double(h5read(fname_templates, '/temp_x') + 1);
        temp_y          = double(h5read(fname_templates, '/temp_y') + 1);
        temp_z          = double(h5read(fname_templates, '/temp_data'));
        templates       = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
        templates_size  = [templates_size(1) templates_size(2) templates_size(3)/2];

        for itemp = 1:numel(clusternr)
            template = full(reshape(templates(:, itemp), templates_size(2), templates_size(1)))';
            [~,i] = max(mean(abs(template),2));
            SpikeRaw{ipart}.template(itemp,:,:) = template;
            SpikeRaw{ipart}.template_maxchan(itemp) = i;
        end

        % create index list for hypnogram
        hypstageindx = ones(1,hdr.nSamples) * -1;
        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}

            % counter to keep track of nr. of samples per directory/file
            dirOnset = 0;

            for idir = 1:length(MuseStruct{ipart})

                if ~isfield(MuseStruct{ipart}{idir},'markers')
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                    continue
                end

                if ~isfield(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']),'synctime')
                    continue
                end

                for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                    y1 = int64(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i) * hdr.Fs + hdr.nSamplesPre);
                    y2 = int64(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i) * hdr.Fs + hdr.nSamplesPre);
                    switch cell2mat(hyplabel)
                        case 'PHASE_1'
                            disp('Found Stage 1');
                            hypstageindx((y1:y2)+dirOnset) = 1;
                        case 'PHASE_2'
                            disp('Found Stage 2');
                            hypstageindx((y1:y2)+dirOnset) = 2;
                        case 'PHASE_3'
                            disp('Found Stage 3');
                            hypstageindx((y1:y2)+dirOnset) = 3;
                        case 'REM'
                            disp('Found Stage REM');
                            hypstageindx((y1:y2)+dirOnset) = 4;
                        case 'AWAKE'
                            disp('Found Stage AWAKE');
                            hypstageindx((y1:y2)+dirOnset) = 0;
                        case 'NO_SCORE'
                            disp('Found Stage NO_SCORE');
                            hypstageindx((y1:y2)+dirOnset) = 0;
                    end
                end

                % assuming all channels have same length in samples
                temp        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.circus.channel{1},'.ncs']));
                hdr_ncs     = ft_read_header(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name));
                dirOnset    = dirOnset + hdr_ncs.nSamples;
            end
        end

        % create trials of x seconds
        cfgtemp                         = [];
        cfgtemp.trl(:,1)                = 1 : hdr.Fs * (cfg.hyp.spikewindow * cfg.hyp.spikewindowoverlap) : hdr.nSamples - hdr.Fs * cfg.hyp.spikewindow;
        cfgtemp.trl(:,2)                = cfgtemp.trl(:,1) + hdr.Fs * cfg.hyp.spikewindow;
        cfgtemp.trl(:,3)                = zeros(size(cfgtemp.trl,1),1);
        cfgtemp.trlunit                 = 'samples';
        cfgtemp.hdr                     = hdr;
        SpikeTrials{ipart}              = ft_spike_maketrials(cfgtemp,SpikeRaw{ipart});
        SpikeTrials{ipart}.trialinfo    = table;

        for i = 1 : size(cfgtemp.trl,1)

            % add sleepstage to trialinfo
            SpikeTrials{ipart}.trialinfo.stage(i)       = median(hypstageindx(cfgtemp.trl(i,1):cfgtemp.trl(i,2)));
            fprintf('Trial %d of %d overlaps with sleepstage %d \n',i,size(SpikeTrials{ipart}.trialinfo,1), SpikeTrials{ipart}.trialinfo.stage(i));

            % clock time of each segment, based on first file of part
            SpikeTrials{ipart}.trialinfo.starttime(i)   = MuseStruct{ipart}{1}.starttime + (i-1) * seconds(cfg.hyp.spikewindow);
            SpikeTrials{ipart}.trialinfo.endtime(i)     = MuseStruct{ipart}{1}.starttime + (i) * seconds(cfg.hyp.spikewindow);

            % add trialnr., for later stratification
            SpikeTrials{ipart}.trialinfo.trialnr(i)     = i;

        end
    end

    save(fname,'SpikeRaw','SpikeTrials');

end % save / force
