function [SpikeRaw] = read_spike_Phy(cfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SpikeRaw] = read_spike_Phy(cfg.phy_directory, cfg.hdr)
%
% Read SpykingCircus analysis results (spike times, templates, amplitudes)
% from Neuralynx data.
% Data must have been converted for Phy GUI, with  parameter prelabelling = True.
% https://spyking-circus.readthedocs.io/en/latest/advanced/extras.html#converting
%
% If data were checked with Phy, neurons are filtered and only the 'good' or
% 'mua' are loaded. Otherwise (data never opened with Phy), all neurons are
% loaded.
%
% Input :
% cfg.hdr           : header of the raw data. Must have fields :
%                     FirstTimeStanmp, TimeStampPerSample, nSamples.
% cfg.phy_directory : directory where the Phy-results files are located
% 
% Output:
% SpikeRaw = raw spike data in FieldTrip raw spike data structure
%
% Dependencies :
% - npy-matlab (https://github.com/kwikteam/npy-matlab)
% - Fieldtrip
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check path
if ~isfolder(cfg.phy_directory), error('Cannot find Phy results directory.'); end

% check cfg.hdr
if ~isfield(cfg.hdr, 'nSamplesPre'),         error('cfg.hdr.nSamplesPre must be specified'); end
if ~isfield(cfg.hdr, 'TimeStampPerSample'),  error('cfg.hdr.TimeStampPerSample must be specified'); end
if ~isfield(cfg.hdr, 'nSamples'),            error('cfg.hdr.nSamples must be specified'); end

%% load spike data

fprintf('Loading spike data from %s\n', cfg.phy_directory);

phydata.cluster_group           = tdfread(fullfile(cfg.phy_directory,'cluster_group.tsv'));   %phy classification.
phydata.spike_times             = readNPY(fullfile(cfg.phy_directory,'spike_times.npy'));     %each timing of any spike, in samples
phydata.spike_templates         = readNPY(fullfile(cfg.phy_directory,'spike_templates.npy')); %for each timing, which (non merged) template. Include the garbage templates
phydata.whitening_mat           = readNPY(fullfile(cfg.phy_directory,'whitening_mat.npy'));   %whitening infos to recover real amplitude data
phydata.whitening_mat_inv       = readNPY(fullfile(cfg.phy_directory,'whitening_mat_inv.npy'));
phydata.templates               = readNPY(fullfile(cfg.phy_directory,'templates.npy'));       %templates waveforms, before merging (merging in phy change clusters but not templates)
phydata.amplitudes              = readNPY(fullfile(cfg.phy_directory,'amplitudes.npy'));      %amplitude of each spike, relative to template

%2 files created by Phy, not Spyking-Circus. Not here if data were not
%opened at least once with Phy :
if exist(fullfile(cfg.phy_directory,'cluster_info.tsv')) && exist(fullfile(cfg.phy_directory,'spike_clusters.npy'))
    ischecked                   = true;
    phydata.cluster_info        = tdfread(fullfile(cfg.phy_directory,'cluster_info.tsv'));%id, amp, ch, depth, fr, group, n_spikes, sh
    phydata.spike_clusters      = readNPY(fullfile(cfg.phy_directory,'spike_clusters.npy')); %for each timing, which (merged) cluster. Include garbage clusters
else
    ischecked                    = false;
    warning('Data were not checked on Phy : loading all templates');
end

%convert templates from 'whitened' to data units
for itemplate = 1:size(phydata.templates,1)
    phydata.templates(itemplate,:,:)= squeeze(phydata.templates(itemplate,:,:)) * phydata.whitening_mat_inv;
end
phydata.templates               = permute(phydata.templates, [1 3 2]);             %permute to be consistent to SpikeRaw struct from MATLAB data

%convert amplitudes from template-normalized to data units
for itemplate = 1:size(phydata.templates,1)
    timings_idx                     = phydata.spike_templates == itemplate;
    template_amplitude              = max(max(abs(phydata.templates(itemplate,:,:)))); %take the abs because some data have negative spikes and some other positive spikes.
    phydata.amplitudes(timings_idx) = phydata.amplitudes(timings_idx) .* template_amplitude;
    clear timings_idx
end

% timestamps  = ft_read_data(hdr_fname,'timestamp','true');
timestamps    =  (cfg.hdr.nSamplesPre : cfg.hdr.nSamples-1-cfg.hdr.nSamplesPre) * cfg.hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster

%% reorganize data to make a Fieldtrip structure

%filter data if checked on Phy. Otherwise load all templates
if ischecked
    cluster_list   = phydata.cluster_group.cluster_id(phydata.cluster_group.group(:,1) ~= 'n')'; %only good, or mua clusters (not noise)
else
    cluster_list   = unique(phydata.spike_templates)';
end

%go trough each cluster
for icluster = 1:size(cluster_list,2)
    
    %add cluster label
    SpikeRaw.label{icluster}                 = sprintf('cluster_%d',cluster_list(icluster));
    
    %find spike time indexes
    clear timings_idx
    if ischecked
        timings_idx                                 = phydata.spike_clusters == cluster_list(icluster);
    else
        timings_idx                                 = phydata.spike_templates == cluster_list(icluster);
    end
    
    %add template waveforms
    cluster_templates                               = unique(phydata.spike_templates(timings_idx))+1; %+1 because template begins at zero but index in phyinfos.templates begins at 1
    SpikeRaw.template{icluster}              = phydata.templates(cluster_templates,:,:);%all template waveforms merged in this cluster
    
    %add template maxchan
    if ischecked
        imaxchan                                    = phydata.cluster_info.ch(phydata.cluster_info.id == cluster_list(icluster)) +1;%+1 because chans idx begin at zero
    else
        [~,imaxchan] = max(mean(abs(SpikeRaw.template{icluster}),3));
    end
    SpikeRaw.template_maxchan(icluster)  = imaxchan - 1; %-1 because chans idx begin at zero with Phy
    
    %add amplitude, samples and timestamps at selected spike timings
    SpikeRaw.amplitude{icluster}             = phydata.amplitudes(timings_idx)';
    SpikeRaw.samples{icluster}               = phydata.spike_times(timings_idx)';
    SpikeRaw.timestamp{icluster}             = timestamps(SpikeRaw.samples{icluster});
    
    %add Phy group info (good, mua)
    if ischecked
        SpikeRaw.cluster_group{icluster}         = phydata.cluster_group.group(phydata.cluster_group.cluster_id == cluster_list(icluster),:);
    end
end
