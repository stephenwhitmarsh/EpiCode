function [SpikeRaw] = readSpikeRaw_Phy(cfg, force)

% function [SpikeRaw] = readSpikeRaw_Phy(cfg, force)
%
% SPIKERAW reads SpykingCircus results (spike times, templates, amplitudes).
% Data must have been converted for Phy GUI, with  parameter prelabelling = True
% https://spyking-circus.readthedocs.io/en/latest/advanced/extras.html#converting
%
% use as
%    [SpikeRaw] = readSpikeRaw_Phy(cfg, force, varargin)
%
% If data were checked with Phy, neurons are filtered and only the 'good' or
% 'mua' are read. The data labeled with noise, or not labelled, are not loaded. 
% Otherwise (data never opened with Phy), all neurons are read.
%
% Necessary input:
% cfg.prefix            = prefix to output files
% cfg.datasavedir       = data directory of Phy converted results
% cfg.circus.channel    = analyzed-electrode names
% cfg.circus.channelname = when using more than one electrode bundle
% force                 = whether to redo analyses or read previous save
%                         (true/false)
%
% ### Optional cfg fields :
% cfg.circus.postfix     = string postfix appended to spike data results.
%                         Default = [].
% cfg.circus.part_list  = list of parts to analyse. Can be an array of
%                         integers, or 'all'. Default = 'all'.
%
% # Output:
% SpikeRaw              = raw spike data in FieldTrip raw spike data structure
%
% Dependencies :
% - npy-matlab (https://github.com/kwikteam/npy-matlab)
% - Fieldtrip

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

% get the default cfg options
cfg.circus.postfix       = ft_getopt(cfg.circus, 'postfix', []);
cfg.circus.part_list     = ft_getopt(cfg.circus, 'part_list', 'all');
cfg.circus.channelname   = ft_getopt(cfg.circus, 'channelname', []);
cfg.circus.maxchan       = ft_getopt(cfg.circus, 'maxchan', 'phy');

% check if depencies is on path, if not add
w = which('readNPY');

% must be in path and not a variable
if isempty(w) || ~isequal(w, 'variable')
    p = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'external', 'npy-matlab-master/npy-matlab');
    fprintf('Adding path: %s\n', p);
    addpath(p);
end

if strcmp(cfg.circus.part_list, 'all')
    cfg.circus.part_list = 1:size(cfg.directorylist, 2);
end

fname = fullfile(cfg.datasavedir, [cfg.prefix, 'SpikeRaw_Phy', cfg.circus.postfix, '.mat']);

if exist(fname, 'file') && force == false
    fprintf('Load precomputed SpikeRaw data\n');
    load(fname, 'SpikeRaw');
    return
end

% concatinate timestamps & samples
[~, samples_separate, ~, hdr_in] = writeSpykingCircusFileList(cfg, false);

% loop through parts
for ipart = cfg.circus.part_list
    
    temp        = sum(samples_separate{ipart});
    samples     = [1 temp(2)];
    clear temp
    
    % use dummy timestamps, 1 per sample, because SC ignores timestamps
    for idir = 1 : size(hdr_in, 2)
        hdr_in{ipart}{idir}.TimeStampPerSample = 1;
        hdr_in{ipart}{idir}.FirstTimeStamp = 0;
    end
    
    if isempty(cfg.circus.channelname)
        channelname = 'none';
    else
        channelname = unique(cfg.circus.channelname);
    end
    
    % to correct maxchan for multiple channels
    channelcount(1) = 0; 
    if ~strcmp(channelname, 'none')    
        for i = 2 : size(channelname, 2)
            channelcount(i) = channelcount(i-1) + sum(channelname{i-1} == string(cfg.circus.channelname));
        end
    end
    
    for chandir = string(channelname)
        
        %% find spiking-circus output path, which is based on the name of the first datafile
        if strcmp(channelname, 'none')   
            datadir = fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)]);
        else
            datadir = fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], char(chandir));
        end
        
        temp = dir(fullfile(datadir, 'SpykingCircus', '*.GUI'));
        if isempty(temp)
            error('Could not find Phy-converted Spyking-Circus results in %s\n', fullfile(datadir, 'SpykingCircus', '*.GUI'));
        else
            phydir = fullfile(temp.folder, temp.name);
        end
        
        %% load spike data
        fprintf('Loading spike data from %s\n', phydir);
        phydata.spike_times             = readNPY(fullfile(phydir, 'spike_times.npy'));     %each timing of any spike, in samples
        phydata.spike_templates         = readNPY(fullfile(phydir, 'spike_templates.npy')); %for each timing, which (non merged) template. Include the garbage templates
        phydata.whitening_mat           = readNPY(fullfile(phydir, 'whitening_mat.npy'));   %whitening infos to recover real amplitude data
        phydata.whitening_mat_inv       = readNPY(fullfile(phydir, 'whitening_mat_inv.npy'));
        phydata.templates               = readNPY(fullfile(phydir, 'templates.npy'));       %templates waveforms, before merging (merging in phy change clusters but not templates)
        phydata.amplitudes              = readNPY(fullfile(phydir, 'amplitudes.npy'));      %amplitude of each spike, relative to template

        if exist(fullfile(phydir, 'spike_clusters.npy'), 'file')
            ischecked                   = true;
            phydata.cluster_group       = tdfread(fullfile(phydir, 'cluster_group.tsv'));   %phy classification.
            try %sometimes this file is not created, I do not understand why
                phydata.cluster_info        = tdfread(fullfile(phydir, 'cluster_info.tsv'));    %id, amp, ch, depth, fr, group, n_spikes, sh
            catch
                phydata.cluster_info = [];
            end
            phydata.spike_clusters      = readNPY(fullfile(phydir, 'spike_clusters.npy'));  %for each timing, which (merged) cluster. Include garbage clusters
            
            % correct difference in field names depending on the Spyking-Circus' version
            try
                phydata.cluster_info.cluster_id = phydata.cluster_info.id;
            catch
            end
        else
            ischecked = false;
            warning('Data were not checked on Phy : loading of all templates');
        end
        
        % convert templates from 'whitened' to data units
        if size(phydata.templates, 3) ~= size(phydata.whitening_mat_inv, 2)
            disp('Size of templates does not match size of whitening matrix - did you accidentally enable ''sparse_export''?');
        else
            for itemplate = 1:size(phydata.templates, 1)
                phydata.templates(itemplate, :, :) = squeeze(phydata.templates(itemplate, :, :)) * phydata.whitening_mat_inv;
            end
            phydata.templates = permute(phydata.templates, [1 3 2]);             %permute to be consistent to SpikeRaw struct from MATLAB data
        end
        
        % convert amplitudes from template-normalized to data units
        for itemplate = 1:size(phydata.templates, 1)
            timings_idx                     = find(phydata.spike_templates == itemplate-1); %-1 because template numerotation starts at zero
            template_amplitude              = max(max(phydata.templates(itemplate, :, :)))-min(min(phydata.templates(itemplate, :, :)));
            for i = timings_idx'
                phydata.amplitudes(i)       = phydata.amplitudes(i) * template_amplitude;
            end
            clear timings_idx
        end

        %% reorganize data to make a Fieldtrip structure
        
        % filter data if checked on Phy. Otherwise load all clusters
        if ischecked
            cluster_list   = phydata.cluster_group.cluster_id(phydata.cluster_group.group(:, 1) ~= 'n')'; %only good, or mua clusters (not noise)
        else
            cluster_list   = unique(phydata.spike_templates)';
        end
        
        % in case no clusters were included
        if isempty(cluster_list)
            if strcmp(channelname, 'none')
                SpikeRaw{ipart} = [];
                continue
            else
                SpikeRaw_chan{ipart}.(char(chandir)) = [];
                continue
            end
        end
        
        % go trough each cluster
        for icluster = 1:size(cluster_list, 2)
            
            % add cluster label
            SpikeRaw{ipart}.label{icluster}                 = sprintf('cluster_%d', cluster_list(icluster));
            SpikeRaw{ipart}.channelname{icluster}           = chandir;
  
            % find spike time indexes
            if ischecked
                cluster_idx                                 = phydata.spike_clusters == cluster_list(icluster);
            else
                cluster_idx                                 = phydata.spike_templates == cluster_list(icluster);
            end
            
            % add template waveforms
            cluster_templates                               = unique(phydata.spike_templates(cluster_idx))+1; %+1 because template begins at zero but index in phyinfos.templates begins at 1
            SpikeRaw{ipart}.template{icluster}              = phydata.templates(cluster_templates, :, :);%all template waveforms merged in this cluster
                        
            % add amplitude, samples and timestamps at selected spike timings
            SpikeRaw{ipart}.amplitude{icluster}             = phydata.amplitudes(cluster_idx)';
            SpikeRaw{ipart}.sample{icluster}                = phydata.spike_times(cluster_idx)';
            SpikeRaw{ipart}.timestamp{icluster}             = SpikeRaw{ipart}.sample{icluster}; % DUMMY TIMESTAMPS!

            % add Phy group info (good, mua)
            try SpikeRaw{ipart}.purity(icluster)        = phydata.cluster_info.purity(phydata.cluster_info.cluster_id   == cluster_list(icluster)); 
            catch; end %not always available, depending of the Spyking-Circus version
            if ischecked
                SpikeRaw{ipart}.cluster_group{icluster}     = phydata.cluster_group.group(phydata.cluster_group.cluster_id  == cluster_list(icluster), :);
            end
            if strcmp(cfg.circus.maxchan, 'phy') && ~isempty(phydata.cluster_info)
                SpikeRaw{ipart}.template_maxchan(icluster)  = phydata.cluster_info.ch(phydata.cluster_info.cluster_id       == cluster_list(icluster)) + channelcount(chandir == string(channelname));
                SpikeRaw{ipart}.template_maxchan_bundle(icluster)  = phydata.cluster_info.ch(phydata.cluster_info.cluster_id == cluster_list(icluster));
            else
                [~, imaxchan] = max(mean(abs(SpikeRaw{ipart}.template{icluster}), 2));
                SpikeRaw{ipart}.template_maxchan(icluster)  = imaxchan - 1 + channelcount(chandir == string(channelname)); % -1 because chans idx begin at zero with Phy
                SpikeRaw{ipart}.template_maxchan_bundle(icluster)  = imaxchan - 1; % -1 because chans idx begin at zero with Phy
            end   
        end
        
        SpikeRaw{ipart}.amplitudedimord = '{chan}_spike';
        SpikeRaw{ipart}.sampledimord    = '{chan}_spike';
        SpikeRaw{ipart}.timestampdimord = '{chan}_spike';
        clear timestamps
        
        %% Convert data into 1-trial Fieldtrip structure to be more consistent
        cfgtemp                     = [];
        cfgtemp.trl                 = [samples, 0];
        cfgtemp.trlunit             = 'samples';
        cfgtemp.hdr                 = rmfield(hdr_in{ipart}{1}, {'nSamples','orig'});
                
        if strcmp(channelname, 'none')         
            SpikeRaw{ipart}                     = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
            SpikeRaw{ipart}.trialinfo           = table;
            SpikeRaw{ipart}.trialinfo.begsample = samples(1);
            SpikeRaw{ipart}.trialinfo.endsample = samples(2);
            SpikeRaw{ipart}.trialinfo.offset    = 0;
            SpikeRaw{ipart}.hdr                 = hdr_in{ipart};             
        else
            SpikeRaw_chan{ipart}.(char(chandir))                     = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
            SpikeRaw_chan{ipart}.(char(chandir)).hdr                 = hdr_in{ipart}; 
            SpikeRaw_chan{ipart}.(char(chandir)).trialinfo           = table;
            SpikeRaw_chan{ipart}.(char(chandir)).trialinfo.begsample = samples(1);
            SpikeRaw_chan{ipart}.(char(chandir)).trialinfo.endsample = samples(2);
            SpikeRaw_chan{ipart}.(char(chandir)).trialinfo.offset    = 0;

            for ilabel = 1 : length(SpikeRaw_chan{ipart}.(char(chandir)).label)
                SpikeRaw_chan{ipart}.(char(chandir)).channelname{ilabel} = char(chandir);
            end
            
            % labels have to have unique names
            for ilabel = 1 : length(SpikeRaw_chan{ipart}.(char(chandir)).label)
                SpikeRaw_chan{ipart}.(char(chandir)).label{ilabel} = char(strcat(SpikeRaw_chan{ipart}.(char(chandir)).label{ilabel}, '_', chandir));
            end
            
            % make sure to clear
            SpikeRaw{ipart} = [];
        end
    end % channelname
    
    % combine different electrode bundles
    if ~strcmp(channelname, 'none')
        
        % check which electrode bundles are not empty
        f = fields(SpikeRaw_chan{ipart});
        fn = [];
        for field = f'
            disp(field)
             if ~isempty(SpikeRaw_chan{ipart}.(char(field)))
                 fn = [fn; field];
             end
        end
        
        % if none are found
        if isempty(fn)
            fprintf('No units found in part %d\n', ipart);
            continue
        end

        % start with first electrode bundle...
        SpikeRaw{ipart} = SpikeRaw_chan{ipart}.(fn{1});
        
        % ...and if that's it, continue...
        if size(fn, 1) == 1
            continue
        end

        % ...else add the rest
        for chandir = string(fn(2:end))'
            for field = string(fields(SpikeRaw{ipart}))'
                if ~any(strcmp(field,{'hdr', 'cfg', 'trialinfo', 'trialtime'})) && ~contains(field, 'dimord')
                    SpikeRaw{ipart}.(field) = [SpikeRaw{ipart}.(field), SpikeRaw_chan{ipart}.(char(chandir)).(field)];
                end
            end
        end
    end
    
end % ipart

save(fname, 'SpikeRaw', '-v7.3');
