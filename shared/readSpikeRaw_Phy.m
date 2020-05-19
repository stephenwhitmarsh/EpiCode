function [SpikeRaw] = readSpikeraw_Phy(cfg,force,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SpikeRaw] = readSpikeraw_Phy(cfg,force,varargin)
% 
% Read SpykingCircus analysis results (spike times, templates, amplitudes).
% Data must have been converted for Phy GUI, with  parameter prelabelling = True.
% https://spyking-circus.readthedocs.io/en/latest/advanced/extras.html#converting
% 
% If data were checked with Phy, neurons are filtered and only the 'good' or
% 'mua' are loaded. Otherwise all neurons are loaded.
%
% ### Necessary input:
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of Phy converted results
% cfg.circus.postfix    = postfix in the name of spike data, if need to
%                         separate several analysis.
% cfg.circus.channel    = analyzed-electrode names
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% ### Output:
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
% 
% Dependencies : 
% - npy-matlab (https://github.com/kwikteam/npy-matlab)
% - Fieldtrip
%
% Paul Baudin (paul.baudin@live.fr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(varargin) || strcmp(varargin{1},'all')
    parts_to_read = 1:size(cfg.directorylist,2);
else
    parts_to_read = varargin{1};
end

fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeRaw_Phy', cfg.circus.postfix, '.mat']);

if exist(fname,'file') && force == false
    fprintf('Load precomputed SpikeRaw data\n');
    load(fname,'SpikeRaw');
    return
end

for ipart = parts_to_read
    
    
    %% find spiking-circus output path, which is based on the name of the first datafile
    
    datadir             = fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)]);
    phydir_temp         = fullfile(datadir,'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*',cfg.circus.postfix,'.GUI']);
    temp                = dir(phydir_temp);
    if isempty(temp)
        error('Could not find Phy-converted Spyking-Circus results: %s\n',phydir);
    end
    phydir              = fullfile(temp.folder,temp.name);
    datafile            = [temp.name(1 : end - length(cfg.circus.postfix) - length('.GUI')), '.ncs'];
    
    %% load spike data
    
    fprintf('Loading spike data from %s\n', phydir);
    
    phydata.amplitudes              = readNPY(fullfile(phydir,'amplitudes.npy'));      %amplitude of each spike (relative to template ?)
    phydata.cluster_group           = tdfread(fullfile(phydir,'cluster_group.tsv'));   %phy classification.
    phydata.spike_times             = readNPY(fullfile(phydir,'spike_times.npy'));     %each timing of any spike, in samples
    phydata.spike_templates         = readNPY(fullfile(phydir,'spike_templates.npy')); %for each timing, which (non merged) template. Include the garbage templates : the last ones, one per electrode
    phydata.templates               = readNPY(fullfile(phydir,'templates.npy'));       %templates waveforms, before merging (merging in phy change clusters but not templates)
    phydata.templates               = permute(phydata.templates, [1 3 2]);             %permute to be consistent to SpikeRaw struct from MATLAB data
    
    %2 files created by Phy, not Spyking-Circus. Not here if data were not
    %opened at least once with Phy :
    if exist(fullfile(phydir,'cluster_info.tsv')) && exist(fullfile(phydir,'spike_clusters.npy'))
        ischecked                   = true;
        phydata.cluster_info        = tdfread(fullfile(phydir,'cluster_info.tsv'));%id, amp, ch, depth, fr, group, n_spikes, sh
        phydata.spike_clusters      = readNPY(fullfile(phydir,'spike_clusters.npy')); %for each timing, which (merged) cluster. Include garbage clusters
    else
        ischecked                    = false;
        warning('Data were not checked on Phy : loading all templates');
    end
    
    % take the first concatinated channel to extract the timestamps
    hdr_fname                        = fullfile(datadir,datafile);
    hdr                              = ft_read_header(hdr_fname);
    % timestamps  = ft_read_data(hdr_fname,'timestamp','true');
    timestamps                       =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster
    
    %% reorganize data to make a structure compatible with Fieldtrip
    
    %filter data if checked on Phy. Otherwise load all templates
    if ischecked     
        cluster_list   = phydata.cluster_group.cluster_id(phydata.cluster_group.group(:,1) ~= 'n')'; %only good, or mua clusters (not noise)
    else
        cluster_list   = unique(phydata.spike_templates)';
    end
    
    for icluster = 1:size(cluster_list,2)
        
        %find spike time indexes for each cluster
        timings_idx                                     = [];
        if ischecked
            timings_idx                                 = phydata.spike_clusters == cluster_list(icluster);
        else
            timings_idx                                 = phydata.spike_templates == cluster_list(icluster);
        end
        
        %add data of selected timings
        SpikeRaw{ipart}.label{icluster}                 = sprintf('cluster_%d',cluster_list(icluster));
        SpikeRaw{ipart}.samples{icluster}               = phydata.spike_times(timings_idx)';
        SpikeRaw{ipart}.amplitude{icluster}             = phydata.amplitudes(timings_idx)';
        SpikeRaw{ipart}.timestamp{icluster}             = timestamps(SpikeRaw{ipart}.samples{icluster});
        
        %Phy group info (good, mua)
        if ischecked
            SpikeRaw{ipart}.cluster_group{icluster}         = phydata.cluster_group.group(phydata.cluster_group.cluster_id == cluster_list(icluster),:); 
        end
        
        %add templates
        %Template dimensions : 1 cluster, 2 electrode, 3 y value
        cluster_templates                               = unique(phydata.spike_templates(timings_idx))+1; %+1 because template begins at zero but index in phyinfos.templates begins at 1
        SpikeRaw{ipart}.template{icluster}              = phydata.templates(cluster_templates,:,:);%all template waveforms merged in this cluster
        if ischecked
            SpikeRaw{ipart}.template_maxchan(icluster)  = phydata.cluster_info.ch(phydata.cluster_info.id == cluster_list(icluster));
        else
            [~,imaxchan] = max(mean(abs(SpikeRaw{ipart}.template{icluster}),3));
            SpikeRaw{ipart}.template_maxchan(icluster)  = imaxchan - 1; %-1 because chans idx begin at zero with Phy
        end
    
    end
    
    clear timestamps
    
end % ipart

save(fname,'SpikeRaw');

end
