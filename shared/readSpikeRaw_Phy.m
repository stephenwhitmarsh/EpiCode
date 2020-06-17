function [SpikeRaw] = readSpikeRaw_Phy(cfg,force,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [SpikeRaw] = readSpikeRaw_Phy(cfg,force,varargin)
% 
% Read SpykingCircus analysis results (spike times, templates, amplitudes).
% Data must have been converted for Phy GUI, with  parameter prelabelling = True.
% https://spyking-circus.readthedocs.io/en/latest/advanced/extras.html#converting
% 
% If data were checked with Phy, neurons are filtered and only the 'good' or
% 'mua' are loaded. Otherwise (data never opened with Phy), all neurons are 
% loaded.
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

%to avoid specificity of our analysis paths : input hdr, phy_datapath, do not save data
%include fieldtrip tools to check datatype and cfg

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
    phydir              = fullfile(datadir,'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1},cfg.circus.postfix,'.GUI']);
%     temp                = dir(phydir_temp);
    if ~isfolder(phydir)
        error('Could not find Phy-converted Spyking-Circus results: %s\n',phydir);
    end
%     phydir              = fullfile(temp.folder,temp.name);
    if ispc, temp=split(phydir, '\'); elseif isunix, temp=split(phydir, '/'); end
    datafile            = [temp{end}(1 : end - length(cfg.circus.postfix) - length('.GUI')), '.ncs'];
    
    %% load spike data
    
    fprintf('Loading spike data from %s\n', phydir);
    
    phydata.cluster_group           = tdfread(fullfile(phydir,'cluster_group.tsv'));   %phy classification.
    phydata.spike_times             = readNPY(fullfile(phydir,'spike_times.npy'));     %each timing of any spike, in samples
    phydata.spike_templates         = readNPY(fullfile(phydir,'spike_templates.npy')); %for each timing, which (non merged) template. Include the garbage templates 
    phydata.whitening_mat           = readNPY(fullfile(phydir,'whitening_mat.npy'));   %whitening infos to recover real amplitude data
    phydata.whitening_mat_inv       = readNPY(fullfile(phydir,'whitening_mat_inv.npy'));  
    phydata.templates               = readNPY(fullfile(phydir,'templates.npy'));       %templates waveforms, before merging (merging in phy change clusters but not templates)
    phydata.amplitudes              = readNPY(fullfile(phydir,'amplitudes.npy'));      %amplitude of each spike, relative to template

    %2 files created by Phy, not Spyking-Circus. Not here if data were not
    %opened at least once with Phy :
    if exist(fullfile(phydir,'cluster_info.tsv')) && exist(fullfile(phydir,'spike_clusters.npy'))
        ischecked                   = true;
        phydata.cluster_info        = tdfread(fullfile(phydir,'cluster_info.tsv'));%id, amp, ch, depth, fr, group, n_spikes, sh
        phydata.spike_clusters      = readNPY(fullfile(phydir,'spike_clusters.npy')); %for each timing, which (merged) cluster. Include garbage clusters
    else
        ischecked                    = false;
        warning('Data were not checked on Phy : loading of all templates');
    end
    
    %convert templates from 'whitened' to data units
    for itemplate = 1:size(phydata.templates,1)
        phydata.templates(itemplate,:,:)= squeeze(phydata.templates(itemplate,:,:)) * phydata.whitening_mat_inv;
    end
    phydata.templates               = permute(phydata.templates, [1 3 2]);             %permute to be consistent to SpikeRaw struct from MATLAB data
    
    %convert amplitudes from template-normalized to data units
    for itemplate = 1:size(phydata.templates,1)
        timings_idx                     = find(phydata.spike_templates == itemplate-1); %-1 because template numerotation starts at zero
        template_amplitude              = max(max(phydata.templates(itemplate,:,:)))-min(min(phydata.templates(itemplate,:,:))); 
        for i = timings_idx'
            phydata.amplitudes(i) = phydata.amplitudes(i) * template_amplitude;
        end
        clear timings_idx
    end
     
    % take the first concatinated channel to extract the timestamps
    hdr_fname                        = fullfile(datadir,datafile);
    hdr                              = ft_read_header(hdr_fname);
    % timestamps  = ft_read_data(hdr_fname,'timestamp','true');
    timestamps                       =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster
    
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
        SpikeRaw{ipart}.label{icluster}                 = sprintf('cluster_%d',cluster_list(icluster));
        
        %find spike time indexes
        clear timings_idx                                    
        if ischecked
            timings_idx                                 = phydata.spike_clusters == cluster_list(icluster);
        else
            timings_idx                                 = phydata.spike_templates == cluster_list(icluster);
        end
        
        %add template waveforms
        cluster_templates                               = unique(phydata.spike_templates(timings_idx))+1; %+1 because template begins at zero but index in phyinfos.templates begins at 1
        SpikeRaw{ipart}.template{icluster}              = phydata.templates(cluster_templates,:,:);%all template waveforms merged in this cluster
        
        %add template maxchan
        if ischecked 
            imaxchan                                    = phydata.cluster_info.ch(phydata.cluster_info.id == cluster_list(icluster)) +1;%+1 because chans idx begin at zero 
        else 
            [~,imaxchan] = max(mean(abs(SpikeRaw{ipart}.template{icluster}),3));
        end
        SpikeRaw{ipart}.template_maxchan(icluster)  = imaxchan - 1; %-1 because chans idx begin at zero with Phy
               
        %add amplitude, samples and timestamps at selected spike timings
        SpikeRaw{ipart}.amplitude{icluster}             = phydata.amplitudes(timings_idx)';
        SpikeRaw{ipart}.samples{icluster}               = phydata.spike_times(timings_idx)';
        SpikeRaw{ipart}.timestamp{icluster}             = timestamps(SpikeRaw{ipart}.samples{icluster});
        
        %add Phy group info (good, mua)
        if ischecked
            SpikeRaw{ipart}.cluster_group{icluster}         = phydata.cluster_group.group(phydata.cluster_group.cluster_id == cluster_list(icluster),:); 
        end
    
    end
    
    clear timestamps
    
end % ipart

save(fname,'SpikeRaw');

end
