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
% 'mua' are read. Otherwise (data never opened with Phy), all neurons are
% read.
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

for ipart = cfg.circus.part_list
    
    if isempty(cfg.circus.channelname)
        channelname = 'none';
    else
        channelname = unique(cfg.circus.channelname);
    end
    
    for chandir = string(channelname)
        
        %% find spiking-circus output path, which is based on the name of the first datafile
        if strcmp(channelname, 'none')   
            datadir = fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)]);
        else
            datadir = fullfile(cfg.datasavedir, cfg.prefix(1:end-1), ['p', num2str(ipart)], char(chandir));
        end
        
        temp    = dir(fullfile(datadir, 'SpykingCircus', '*.GUI'));
        if isempty(temp)
            error('Could not find Phy-converted Spyking-Circus results: %s\n', phydir);
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
        
        %2 files created by Phy, not Spyking-Circus. Not here if data were not
        %opened at least once with Phy :
        if exist(fullfile(phydir, 'cluster_info.tsv'), 'file') && exist(fullfile(phydir, 'spike_clusters.npy'), 'file')
            ischecked                   = true;
            phydata.cluster_group       = tdfread(fullfile(phydir, 'cluster_group.tsv'));   %phy classification.
            phydata.cluster_info        = tdfread(fullfile(phydir, 'cluster_info.tsv'));%id, amp, ch, depth, fr, group, n_spikes, sh
            phydata.spike_clusters      = readNPY(fullfile(phydir, 'spike_clusters.npy')); %for each timing, which (merged) cluster. Include garbage clusters
            
            %correct difference in field names depending on the Spyking-Circus' version
            try
                phydata.cluster_info.cluster_id = phydata.cluster_info.id;
            catch
            end
        else
            ischecked = false;
            warning('Data were not checked on Phy : loading of all templates');
        end
        
        %convert templates from 'whitened' to data units
        if size(phydata.templates, 3) ~= size(phydata.whitening_mat_inv, 2)
            error('Size of templates does not match size of whitening matrix - did you accidentally enable ''sparse_export''?');
        end
        for itemplate = 1:size(phydata.templates, 1)
            phydata.templates(itemplate, :, :) = squeeze(phydata.templates(itemplate, :, :)) * phydata.whitening_mat_inv;
        end
        phydata.templates = permute(phydata.templates, [1 3 2]);             %permute to be consistent to SpikeRaw struct from MATLAB data
        
        %convert amplitudes from template-normalized to data units
        for itemplate = 1:size(phydata.templates, 1)
            timings_idx                     = find(phydata.spike_templates == itemplate-1); %-1 because template numerotation starts at zero
            template_amplitude              = max(max(phydata.templates(itemplate, :, :)))-min(min(phydata.templates(itemplate, :, :)));
            for i = timings_idx'
                phydata.amplitudes(i)       = phydata.amplitudes(i) * template_amplitude;
            end
            clear timings_idx
        end
        
        %% take the first concatinated channel to extract the timestamps
        temp        = dir(fullfile(datadir, '*.ncs'));
        hdr_fname   = fullfile(datadir, temp(1).name);
        hdr         = ft_read_header(hdr_fname);  
        timestamps  = ft_read_data(hdr_fname, 'timestamp', 'true');
        
        %% The following is dangerous as missing data in e.g. Neuralynx files will result in unexpected offsets
        %     timestamps  =  (0:hdr.nSamples-1) * hdr.TimeStampPerSample; % calculate timemstamps myself, as this is much faster
        
        %% reorganize data to make a Fieldtrip structure
        
        % filter data if checked on Phy. Otherwise load all clusters
        if ischecked
            cluster_list   = phydata.cluster_group.cluster_id(phydata.cluster_group.group(:, 1) ~= 'n')'; %only good, or mua clusters (not noise)
        else
            cluster_list   = unique(phydata.spike_templates)';
        end
        
        % in case no clusters were included
        if isempty(cluster_list)
            continue
        end
        
        %go trough each cluster
        for icluster = 1:size(cluster_list, 2)
            
            % add cluster label
            SpikeRaw{ipart}.label{icluster}                 = sprintf('cluster_%d', cluster_list(icluster));
            
            % find spike time indexes
            clear timings_idx
            if ischecked
                timings_idx                                 = phydata.spike_clusters == cluster_list(icluster);
            else
                timings_idx                                 = phydata.spike_templates == cluster_list(icluster);
            end
            
            % add template waveforms
            cluster_templates                               = unique(phydata.spike_templates(timings_idx))+1; %+1 because template begins at zero but index in phyinfos.templates begins at 1
            SpikeRaw{ipart}.template{icluster}              = phydata.templates(cluster_templates, :, :);%all template waveforms merged in this cluster
            
            % add template maxchan
            if ischecked
                imaxchan                                    = phydata.cluster_info.ch(phydata.cluster_info.cluster_id == cluster_list(icluster)) +1; %+1 because chans idx begin at zero
            else
                [~, imaxchan] = max(mean(abs(SpikeRaw{ipart}.template{icluster}), 3));
            end
            SpikeRaw{ipart}.template_maxchan(icluster)      = imaxchan - 1; %-1 because chans idx begin at zero with Phy
            
            % add amplitude, samples and timestamps at selected spike timings
            SpikeRaw{ipart}.amplitude{icluster}             = phydata.amplitudes(timings_idx)';
            SpikeRaw{ipart}.sample{icluster}                = phydata.spike_times(timings_idx)';
            SpikeRaw{ipart}.timestamp{icluster}             = timestamps(SpikeRaw{ipart}.sample{icluster});
            
            % add Phy group info (good, mua)
            if ischecked
                SpikeRaw{ipart}.cluster_group{icluster}     = phydata.cluster_group.group(phydata.cluster_group.cluster_id == cluster_list(icluster), :);
                SpikeRaw{ipart}.template_maxchan(icluster)  = phydata.cluster_info.ch(phydata.cluster_info.cluster_id == cluster_list(icluster));
                %new field with purity index in 0.9.9 version of Spyking-Circus
                try SpikeRaw{ipart}.purity(icluster)        = phydata.cluster_info.purity(phydata.cluster_info.cluster_id == cluster_list(icluster)); end
            else
                [~, imaxchan] = max(mean(abs(SpikeRaw{ipart}.template{icluster}), 3));
                SpikeRaw{ipart}.template_maxchan(icluster)  = imaxchan - 1; %-1 because chans idx begin at zero with Phy
            end
            
        end
        
        clear timestamps
        
        %% Convert data into 1-trial Fieldtrip structure to be more consistent
        %if no need to have trials (otherwise some scripts written for trial
        %structures would not be compatible with this 'raw' structure, ie
        %ft_spike_isi.m).
        filebegin                   = 0-hdr.nSamplesPre;
        fileend                     = hdr.nSamples-hdr.nSamplesPre;
        cfgtemp                     = [];
        cfgtemp.trl                 = [filebegin, fileend, 0];
        cfgtemp.trlunit             = 'samples';
        cfgtemp.hdr                 = hdr;
        cfgtemp.timestampspersecond = hdr.TimeStampPerSample * hdr.Fs;
        
        if strcmp(channelname, 'none')         
            SpikeRaw{ipart}                     = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
            SpikeRaw{ipart}.hdr                 = hdr;
            SpikeRaw{ipart}.trialinfo           = table;
            SpikeRaw{ipart}.trialinfo.begsample = filebegin;
            SpikeRaw{ipart}.trialinfo.endsample = fileend;
            SpikeRaw{ipart}.trialinfo.offset    = 0;
        else
            SpikeRaw_temp{ipart}.(char(chandir))                     = ft_spike_maketrials(cfgtemp, SpikeRaw{ipart});
            SpikeRaw_temp{ipart}.(char(chandir)).hdr                 = hdr;
            SpikeRaw_temp{ipart}.(char(chandir)).trialinfo           = table;
            SpikeRaw_temp{ipart}.(char(chandir)).trialinfo.begsample = filebegin;
            SpikeRaw_temp{ipart}.(char(chandir)).trialinfo.endsample = fileend;
            SpikeRaw_temp{ipart}.(char(chandir)).trialinfo.offset    = 0;
            
            % labels have to have unique names
            for ilabel = 1 : length(SpikeRaw_temp{ipart}.(char(chandir)).label)
                SpikeRaw_temp{ipart}.(char(chandir)).label{ilabel} = char(strcat(SpikeRaw_temp{ipart}.(char(chandir)).label{ilabel},'_',chandir));
            end
            clear SpikeRaw
        end
                    
    end % channelname
    
    % combine diffrent electrodebundles - very smartly! :-)
    if ~strcmp(channelname, 'none')
        try % in case no clusters were selected
            SpikeRaw{ipart} = SpikeRaw_temp{ipart}.(channelname{1});
            for chandir = string(channelname(2:end))
                for field = string(fields(SpikeRaw{ipart}))'
                    if ~strcmp(field,{'hdr','cfg','trialinfo','trialtime'})
                        SpikeRaw{ipart}.(field) = [SpikeRaw{ipart}.(field), SpikeRaw_temp{ipart}.(char(chandir)).(field)];
                    end
                end
            end
        catch
        end
    end
    
end % ipart

save(fname, 'SpikeRaw');
