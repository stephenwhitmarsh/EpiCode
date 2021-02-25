function [SpikeRaw] = readSpikeRaw_SpykingCircus(cfg,force,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [SpikeRaw] = readSpikeRaw_SpykingCircus(cfg,force,varargin)
% Read results from Spyking-Circus analysis; spike-times and templates.
%
% Necessary input:
%
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of results
% cfg.circus.outputdir  = directory delow datasavedir for spyking-circus
% cfg.circus.suffix     = from spyking-circus params files
% cfg.circus.channel    = micro electrode names
%
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
%
% Output:
%
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
%
% Stephen Whitmarsh (stephen.whitmarsh@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) || strcmp(varargin{1},'all')
    parts_to_read = 1:size(cfg.directorylist,2);
else
    parts_to_read = varargin{1};
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeRaw_SC', cfg.circus.postfix, '.mat']);
if exist(fname,'file') && force == false
    load(fname,'SpikeRaw');
else

    for ipart = parts_to_read

        % find spiking-circus output path, which is based on the name of the
        % first datafile
        temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.result',cfg.circus.postfix,'.hdf5']));
        if isempty(temp)
            fprintf('Could not find Spyking-Circus results: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.result',cfg.circus.postfix,'.hdf5']));
            return
        else
            fname_spikes = fullfile(temp.folder,temp.name);
        end

        temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.templates',cfg.circus.postfix,'.hdf5']));
        if isempty(temp)
            fprintf('Could not find Spyking-Circus templates: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.templates',cfg.circus.postfix,'.hdf5']));
            return
        end
        fname_templates = fullfile(temp.folder,temp.name);

                % load spiking data
        if exist(fname_spikes,'file')
            fprintf('Loading spike data from: %s\n',fname_spikes);
            datinfo     = h5info(fname_spikes);
            temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
            hdr_fname   = fullfile(temp(1).folder,temp(1).name);
            hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
            %             tempdataset{1} = hdr_fname;
            %             temp        = ft_read_neuralynx_interp(tempdataset);
%                         timestamps  = ft_read_data(hdr_fname,'timestamp','true');  % take the first concatinated channel to extract the timestamps
            %          	  timestamps  =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster

            clear names
            % read spiketimes of clusters
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
                % read spike timings (in seconds)
                datasetname = char(strcat('/spiketimes/',SpikeRaw{ipart}.label{i}));
                SpikeRaw{ipart}.sample{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0

                % read amplitudes
                datasetname = char(strcat('/amplitudes/',SpikeRaw{ipart}.label{i}));
                SpikeRaw{ipart}.amplitude{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0

                % map samplenrs onto timestamps
%                 SpikeRaw{ipart}.timestamp{clusternr(i)+1} = timestamps(SpikeRaw{ipart}.sample{clusternr(i)+1});
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
                SpikeRaw{ipart}.template{itemp} = template;
                SpikeRaw{ipart}.template_maxchan(itemp) = i;
            end
        end % if exists spike file.
    end
    save(fname,'SpikeRaw');

end % save / force
