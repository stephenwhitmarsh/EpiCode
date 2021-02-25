function [SpikeRaw, SpikeTrials] = readSpykingCircus_selected_stim(cfg)

addpath /network/lustre/iss01/charpier/stephen.whitmarsh/WhitmarshEpilepsy
addpath /network/lustre/iss01/charpier/stephen.whitmarsh/fieldtrip/
ft_defaults

%% use spiking circus output

% SpikeRaw{height(cfg.markers)} = [];
SpikeTrials{height(cfg.markers)} = [];

for stimrate = [50]
    
    % find spiking-circus path
    temp = dir(fullfile(cfg.datasavedir,'SpykingCircus',[num2str(stimrate),'Hz_concatinated_',cfg.channel{1}(1:end-2),'_*.result',cfg.suffix,'.hdf5']));
    fname_spikes = fullfile(temp.folder,temp.name);
    
    temp = dir(fullfile(cfg.datasavedir,'SpykingCircus',[num2str(stimrate),'Hz_concatinated_',cfg.channel{1}(1:end-2),'_*.templates',cfg.suffix,'.hdf5']));
    fname_templates = fullfile(temp.folder,temp.name);
    
    
    
    % load spiking data
    if exist(fname_spikes,'file')
        fprintf('Loading spike data from: %s\n',fname_spikes);
        datinfo     = h5info(fname_spikes);
        temp        = dir(fullfile(cfg.datasavedir,[num2str(stimrate),'Hz_concatinated_',cfg.channel{1}(1:end-2),'_*.ncs']));
        hdr         = ft_read_header(fullfile(temp(1).folder,temp(1).name)); % take the first file to extract the header of the data
        timestamps  = ft_read_data(fullfile(temp(1).folder,temp(1).name),'timestamp','true');  % take the first concatinated file to extract the timestamps
           
        % read spiketimes of clusters
        for i = 1 : size(datinfo.Groups,1)
            names(i) = string(datinfo.Groups(i).Name);
            if strfind(names(i),'spiketimes')
                spiketimes_indx = i;
            end
        end
        for i = 1 : size(datinfo.Groups(spiketimes_indx).Datasets,1)
            SpikeRaw.label{i} = datinfo.Groups(spiketimes_indx).Datasets(i).Name;
            temp = strsplit(SpikeRaw.label{i},'_');
            clusternr(i) = str2num(temp{2});
        end
        
        for i = 1:numel(clusternr)
            % read spike timings (in seconds)
            datasetname = char(strcat('/spiketimes/',SpikeRaw.label{i}));
            SpikeRaw.sample{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % read amplitudes
            datasetname = char(strcat('/amplitudes/',SpikeRaw.label{i}));
            SpikeRaw.amplitude{clusternr(i)+1} = h5read(fname_spikes,datasetname); % count from 1 instead of 0
            
            % map samplenrs onto timestamps
%             SpikeRaw.timestamp{i} = timestamps(SpikeRaw.sample{i});
            SpikeRaw.timestamp{i} = SpikeRaw.sample{i} * hdr.TimeStampPerSample;
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
            SpikeRaw.template(itemp,:,:) = template;
            SpikeRaw.template_maxchan(itemp) = i;
        end
        
        % create spiketrials timelocked to events
        cfgtemp                 = [];
        cfgtemp.trl             = cfg.trialinfo;
        cfgtemp.trl(:,3)        = -ones(size(cfgtemp.trl,1),1)*cfg.prestim*hdr.Fs;
        cfgtemp.trlunit         = 'samples';
        cfgtemp.hdr             = hdr;
        SpikeTrials             = ft_spike_maketrials(cfgtemp,SpikeRaw);
        
        % plot spike data
        cfgtemp                 = [];
        cfgtemp.timwin          = [-0.5 0.5];
%         cfgtemp.fsample         = cfg.resamplefs; % sample at 1000 hz
        sdf                     = ft_spikedensity(cfgtemp,SpikeTrials);
        
        for itemp = 1 : size(SpikeTrials.label,2)
            
            fig = figure;
            
            subplot(1,2,1);
            cfgtemp                 = [];
            cfgtemp.spikechannel    = itemp;
            cfgtemp.topplotsize     = 0.2;
            cfgtemp.topplotfunc     = 'line'; % plot as a line
            cfgtemp.latency         = [-cfg.prestim, cfg.poststim];
            cfgtemp.spikelength     = 1;
            cfgtemp.linewidth       = 1;
            cfgtemp.errorbars       = 'sem'; % plot with the standard deviation
            cfgtemp.interactive     = 'yes'; % toggle off interactive mode
            ft_spike_plot_raster(cfgtemp,SpikeTrials, sdf)
            title('Spike density and rasterplot');
            
            subplot(1,2,2); hold;
            maxheight = max(max(abs(SpikeRaw.template(itemp,:,:))));
            for ichan = 1 : size(SpikeRaw.template,2)
                if ichan == SpikeRaw.template_maxchan(itemp)
                    plot(squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'r');
                else
                    plot(squeeze(SpikeRaw.template(itemp,ichan,:))+(ichan)*maxheight,'k');
                    
                end
            end
            ax = axis;
            axis([ax(1), ax(2), 0, ax(4)]);
            title('Template');
            set(gca, 'YTickLabel', '');
            
            set(fig,'PaperOrientation','landscape');
            set(fig,'PaperUnits','normalized');
            set(fig,'PaperPosition', [0 0 1 1]);
            print(fig, '-dpdf', fullfile(cfg.imagesavedir,[num2str(stimrate),'Hz_raster_template_',num2str(itemp),'.pdf']),'-r600');
        end
    else
        fprintf('Could not find: %s or %s\n',fname_spikes,fname_templates);
    end
end




