function [LFP] = readLFP(cfg, MuseStruct, force)

% READLFP reads data epoched according to markers from MuseStruct.
% Use as
%   [LFP] = readLFP(cfg, MuseStruct, force)
%
% Will also mark overlap with hypnogram markers
% Optionally downsample to same samplerate if data contains multiple
% samplerates to allow combination of channels in single data structure
% Markernames that contain a space (' ') or minus ('-') are
% replaced with an underscore ('_').

% Necessary LFP fields (with example):
% cfg.LFP.name                = {'HLFP'};
% cfg.LFP.channel             = {'_Ha2g_1', '_Ha2g_2', '_Ha2g_3', '_Ha2g_4'}
% cfg.LFP.write               = true (write LFP to disk, default = false)
%
% Necessary generic epoch fields (with example):
% cfg.epoch.toi{1}            = [-0.5  1];
% cfg.epoch.toi{2}            = [-0.5  1];
% cfg.epoch.pad               = {0.5, 0.5, 0.5};
%
% Optional fields (with example, defaults = false):
% cfg.LFP.lpfilter            = 'yes';
% cfg.LFP.hpfilter            = 'no';
% cfg.LFP.bpfilter            = 'no';
% cfg.LFP.bsfilter            = 'no';
% cfg.LFP.lpfreq              = 40;
% cfg.LFP.hpfreq              = [];
% cfg.LFP.bpfreq              = [];
% cfg.LFP.bsfreq              = [];
% cfg.LFP.resamplefs          = 1000;
% cfg.LFP.demean              = 'yes';
% cfg.LFP.baselinewindow      = [-0.5 0];
%
% Dependencies: dir2.m, recent FieldTrip version

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

write               = ft_getopt(cfg.LFP, 'write', true);
keepcfg             = ft_getopt(cfg.LFP, 'keepcfg', true);
cfg.LFP.reref       = ft_getopt(cfg.LFP, 'reref', 'no');
cfg.LFP.rerefmethod = ft_getopt(cfg.LFP, 'rerefmethod', []);
cfg.LFP.refchannel  = ft_getopt(cfg.LFP, 'refchannel', []);
cfg.LFP.postfix     = ft_getopt(cfg.LFP, 'postfix', []);
cfg.LFP.overlap     = ft_getopt(cfg.LFP, 'overlap', []);

% add markers to always look for overlap for
cfg.LFP.overlap                 = unique([cfg.LFP.overlap, "BAD", "PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "PRESLEEP", "POSTSLEEP"], 'stable');
cfg.muse.startmarker.BAD        = 'BAD__START__';
cfg.muse.endmarker.BAD          = 'BAD__END__';
cfg.muse.startmarker.PHASE_1    = 'PHASE_1__START__';
cfg.muse.endmarker.PHASE_1      = 'PHASE_1__END__';
cfg.muse.startmarker.PHASE_2    = 'PHASE_2__START__';
cfg.muse.endmarker.PHASE_2      = 'PHASE_2__END__';
cfg.muse.startmarker.PHASE_3    = 'PHASE_3__START__';
cfg.muse.endmarker.PHASE_3      = 'PHASE_3__END__';
cfg.muse.startmarker.REM        = 'REM__START__';
cfg.muse.endmarker.REM          = 'REM__END__';
cfg.muse.startmarker.AWAKE      = 'AWAKE__START__';
cfg.muse.endmarker.AWAKE        = 'AWAKE__END__';
cfg.muse.startmarker.PRESLEEP   = 'PRESLEEP__START__';
cfg.muse.endmarker.PRESLEEP     = 'PRESLEEP__END__';
cfg.muse.startmarker.POSTSLEEP  = 'POSTSLEEP__START__';
cfg.muse.endmarker.POSTSLEEP    = 'POSTSLEEP__END__';

% used later to determine sleepstage
hyplabels                       = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "PRESLEEP", "POSTSLEEP"];
hypindex                        = false(length(cfg.LFP.overlap), 1);
for i = 1 : length(cfg.LFP.overlap)
    if any(strcmp(cfg.LFP.overlap{i}, hyplabels))
        hypindex(i) = true;
    end
end

if nargin == 1
    for markername = string(cfg.LFP.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'LFP_', markername, cfg.LFP.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                    for ipart = 1 : size(cfg.directorylist, 2)
                        LFP{ipart}.(markername) = temp.LFP{ipart}.(markername);
                    end
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
        else
            fprintf('Will be (re-) computing LFP data for %s\n', markername);
        end
    end
    return
elseif ~force
    missing = [];
    for markername = string(cfg.LFP.name)
        fname = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'LFP_', markername, cfg.LFP.postfix, '.mat'));
        if exist(fname, 'file')
            fprintf('Reading %s\n', fname);
            count = 0;
            err_count = 0;
            while count == err_count
                try
                    temp = load(fname);
                    for ipart = 1 : size(cfg.directorylist, 2)
                        LFP{ipart}.(markername) = temp.LFP{ipart}.(markername);
                    end
                catch ME
                    err_count = err_count + 1;
                end
                count = count + 1;
            end
        else
            fprintf('Will be (re-) computing LFP data for %s\n', markername);
            missing = [missing, markername];
        end
    end
    cfg.LFP.name = missing;
end

% get file format
[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

% % initialize LFP, to return empty cell in case there is no LFP to load
% hyplabels = ["PHASE_1", "PHASE_2", "PHASE_3", "REM", "AWAKE", "NO_SCORE"];

% loop over markers
for markername = string(cfg.LFP.name)
    
    fname_out = fullfile(cfg.datasavedir, strcat(cfg.prefix, 'LFP_', markername, cfg.LFP.postfix, '.mat'));
    
    % loop over parts within subject
    for ipart = 1 : size(MuseStruct, 2)
        
        fprintf('For marker %s\n', cell2mat(markername));
        hasmarker = false(length(MuseStruct{ipart}), 1);
        
        % loop over directories
        for idir = 1 : length(MuseStruct{ipart})
            
            if ~isfield(MuseStruct{ipart}{idir}, 'markers')
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)), 'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime)
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.endmarker.(markername))
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)), 'synctime')
                continue
            end
            if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime)
                continue
            end
            
            if isNeuralynx
                nfile = size(cfg.LFP.channel, 2); % one file per channel
            elseif isMicromed
                nfile = 1; % only one file with all electrodes
                fname = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.TRC']);
            elseif isBrainvision
                nfile = 1; % only one file with all electrodes
                fname = fullfile(cfg.rawdir, [cfg.directorylist{ipart}{idir} '.eeg']);
            end
            
            % loop over files (electrodes)
            for ifile = 1 : nfile
                
                %load data
                if isNeuralynx
                    temp                = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', cfg.LFP.channel{ifile}, '.ncs']));
                    cfgtemp             = [];
                    cfgtemp.dataset     = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                    dat                 = ft_preprocessing(cfgtemp);
                    
                    % rereferencing
                    if strcmp(cfg.LFP.reref, 'yes') && strcmp(cfg.LFP.rerefmethod, 'bipolar')
                        b               = cfg.LFP.channel{ifile}(1:end-1);
                        n               = num2str(str2double(cfg.LFP.channel{ifile}(end)) + 1);
                        temp            = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', b, n, '.ncs']));
                        cfgtemp         = [];
                        cfgtemp.dataset = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
                        refdat          = ft_preprocessing(cfgtemp);
                        dat.trial{1}    = dat.trial{1} - refdat.trial{1};
                        clear refdat
                    end
                    
                else
                    cfgtemp             = [];
                    cfgtemp.dataset     = fname;
                    cfgtemp.channel     = cfg.LFP.channel;
                    dat                 = ft_preprocessing(cfgtemp);
                end
                
                % filtering
                cfgtemp                 = [];
                cfgtemp.lpfilter        = ft_getopt(cfg.LFP, 'lpfilter', 'no');
                cfgtemp.hpfilter        = ft_getopt(cfg.LFP, 'hpfilter', 'no');
                cfgtemp.bpfilter        = ft_getopt(cfg.LFP, 'bpfilter', 'no');
                cfgtemp.bsfilter        = ft_getopt(cfg.LFP, 'bsfilter', 'no');
                cfgtemp.dftfilter       = ft_getopt(cfg.LFP, 'dftfilter', 'no');
                
                cfgtemp.lpfreq          = ft_getopt(cfg.LFP, 'lpfreq', []);
                cfgtemp.hpfreq          = ft_getopt(cfg.LFP, 'hpfreq', []);
                cfgtemp.bpfreq          = ft_getopt(cfg.LFP, 'bpfreq', []);
                cfgtemp.bsfreq          = ft_getopt(cfg.LFP, 'bsfreq', []);
                dat                     = ft_preprocessing(cfgtemp, dat);
                
                % append EMG data (if any)
                if isMicromed || isBrainvision %not adapted for nlx data (1 file per electrode) for now
                    
                    cfg.EMG = ft_getopt(cfg, 'EMG', []);
                    
                    cfgtemp                   = [];
                    cfgtemp.channel           = ft_getopt(cfg.EMG, sprintf('%s',markername)); % load the emg associated with eeg marker, and the ref if any
                    
                    if ~isempty(cfgtemp.channel)
                        cfgtemp.dataset           = fname;
                        cfgtemp.reref             = ft_getopt(cfg.EMG, 'reref', 'no');
                        cfgtemp.rerefmethod       = ft_getopt(cfg.EMG, 'rerefmethod', []);
                        cfgtemp.refchannel        = ft_getopt(cfg.EMG, 'refchannel', []);
                        data_EMG                  = ft_preprocessing(cfgtemp);
                        
                        % filtering
                        cfgtemp                 = [];
                        cfgtemp.lpfilter        = ft_getopt(cfg.EMG, 'lpfilter', 'no');
                        cfgtemp.hpfilter        = ft_getopt(cfg.EMG, 'hpfilter', 'no');
                        cfgtemp.bpfilter        = ft_getopt(cfg.EMG, 'bpfilter', 'no');
                        cfgtemp.bsfilter        = ft_getopt(cfg.EMG, 'bsfilter', 'no');
                        cfgtemp.lpfreq          = ft_getopt(cfg.EMG, 'lpfreq', []);
                        cfgtemp.hpfreq          = ft_getopt(cfg.EMG, 'hpfreq', []);
                        cfgtemp.bpfreq          = ft_getopt(cfg.EMG, 'bpfreq', []);
                        cfgtemp.bsfreq          = ft_getopt(cfg.EMG, 'bsfreq', []);
                        cfgtemp.lpfiltord       = ft_getopt(cfg.EMG, 'lpfiltord', []);
                        cfgtemp.hpfiltord       = ft_getopt(cfg.EMG, 'hpfiltord', []);
                        cfgtemp.bpfiltord       = ft_getopt(cfg.EMG, 'bpfiltord', []);
                        cfgtemp.bsfiltord       = ft_getopt(cfg.EMG, 'bsfiltord', []);
                        data_EMG                = ft_preprocessing(cfgtemp,data_EMG);
                        
                        % append EMG to EEG data
                        cfgtemp                 = [];
                        cfgtemp.keepsampleinfo  = 'yes';
                        dat                     = ft_appenddata(cfgtemp, dat, data_EMG);
                    end
                end
                
                % downsample data and correct baseline
                if isfield(cfg.LFP, 'resamplefs')
                    
                    cfgtemp                         = [];
                    cfgtemp.resamplefs              = ft_getopt(cfg.LFP, 'resamplefs', []);
                    if strcmp(ft_getopt(cfg.LFP, 'baseline', 'no'), 'no')
                        cfgtemp.demean              = 'no';
                    else
                        cfgtemp.demean              = 'yes';
                        cfgtemp.baselinewindow      = cfg.LFP.baselinewindow.(markername);
                    end
                    dat                             = ft_resampledata(cfgtemp, dat);
                end
                
                % create trial segmentation common to resampled
                % data. Neuralynx : same markers for all files
                % of one dir.
                Startsample             = [];
                Endsample               = [];
                Offset                  = [];
                
                % add trialnr before sorting for reference & checks
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr     = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr       = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime, 2);
                
                % sort times in case markers have been combined
                [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(sidx);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(sidx);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr     = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr(sidx);
                
                [~, sidx] = sort(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime      = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(sidx);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock         = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(sidx);
                MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr       = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr(sidx);
                
                overlap_sec = zeros(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2), size(cfg.LFP.overlap, 2));
                overlap_cnt = zeros(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2), size(cfg.LFP.overlap, 2));
                
                % loop over events
                for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime, 2)
                    
                    ss = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) * dat.fsample);
                    if strcmp(cfg.muse.startmarker.(markername), cfg.muse.endmarker.(markername))
                        es = ss;
                        idx_end = ievent;
                    else
                        idx_end = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime * dat.fsample) >= ss, 1, 'first');
                        es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(idx_end) * dat.fsample);
                    end
                    
                    if isempty(es)
                        continue
                    end
                    
                    Startsample(ievent) = ss + cfg.epoch.toi.(markername)(1) * dat.fsample - cfg.epoch.pad.(markername) * dat.fsample;
                    Endsample(ievent)   = es + cfg.epoch.toi.(markername)(2) * dat.fsample + cfg.epoch.pad.(markername) * dat.fsample;
                    Offset(ievent)      = (cfg.epoch.toi.(markername)(1) - cfg.epoch.pad.(markername)) * dat.fsample;
                    
                    % extra info for trialinfo
                    Trialnr(ievent)         = ievent;
                    Filenr(ievent)          = ifile;
                    Directory(ievent,:)     = cfg.directorylist{ipart}{idir};
                    StartTrialnr(ievent)    = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).trialnr(ievent);
                    EndTrialnr(ievent)      = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).trialnr(idx_end);
                    Starttime(ievent)       = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).clock(ievent) + seconds(cfg.epoch.toi.(markername)(1));
                    Endtime(ievent)         = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).clock(idx_end)   + seconds(cfg.epoch.toi.(markername)(2));
                    
                    % will be used to find overlap between events
                    trlstart                = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(markername)).synctime(ievent) + cfg.epoch.toi.(markername)(1);
                    trlend                  = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(markername)).synctime(idx_end) + cfg.epoch.toi.(markername)(2);       

                    % find overlap
                    for iother = 1 : size(cfg.LFP.overlap, 2)
                        
                        % if the event is not present continue
                        if ~isfield(MuseStruct{ipart}{idir}.markers, cfg.muse.startmarker.(cfg.LFP.overlap{iother}))
                            continue
                        end
                        if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.LFP.overlap{iother})), 'synctime')
                            continue
                        end
                 
                        % loop over instances of overlap-event
                        for i = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.LFP.overlap{iother})).synctime, 2)
                            
                            other_start = MuseStruct{ipart}{idir}.markers.(cfg.muse.startmarker.(cfg.LFP.overlap{iother})).synctime(i);
                            other_end   = MuseStruct{ipart}{idir}.markers.(cfg.muse.endmarker.(cfg.LFP.overlap{iother})).synctime(i);
                            
                            % trial falls fully within other
                            if other_start < trlstart && other_end > trlend
                                overlap_sec(ievent, iother) = overlap_sec(ievent, iother) + (trlend - trlstart);
                                overlap_cnt(ievent, iother) = overlap_cnt(ievent, iother) + 1;
                                % other falls fully within trial
                            elseif other_start > trlstart && other_end < trlend
                                overlap_sec(ievent, iother) = overlap_sec(ievent, iother) + (other_end - other_start);
                                overlap_cnt(ievent, iother) = overlap_cnt(ievent, iother) + 1;
                                % trial overlaps with beginning of other
                            elseif other_start > trlstart && other_start < trlend
                                overlap_sec(ievent, iother) = overlap_sec(ievent, iother) + (trlend - other_start);
                                overlap_cnt(ievent, iother) = overlap_cnt(ievent, iother) + 1;
                                % trial overlaps with end of other
                            elseif other_start < trlstart && other_end > trlstart
                                overlap_sec(ievent, iother) = overlap_sec(ievent, iother) + (other_end - trlstart);
                                overlap_cnt(ievent, iother) = overlap_cnt(ievent, iother) + 1;
                            end
                            
                        end % i
                    end % iother
                    
                    % add sleep stage to trialinfo
                    sel = overlap_sec(ievent, hypindex);
                    [val, indx] = max(sel);
                    if val > 0
                        hyplabels_trl(ievent) = hyplabels(indx);
                    else
                        hyplabels_trl(ievent) = "NO_SCORE";
                    end
   
                end
                
                % create Fieldtrip trl
                full_trial  = Startsample > 0 & Endsample < length(dat.trial{1}); % don't read before BOF or after EOF
                cfgtemp     = [];
                cfgtemp.trl = round([Startsample; Endsample; Offset]');
                cfgtemp.trl = cfgtemp.trl(full_trial, :);
                
                if isempty(cfgtemp.trl); continue; end
                
                filedat{ifile} = ft_redefinetrial(cfgtemp, dat);
                clear dat

                if ifile == 1
                    trialinfo            = table;
                    trialinfo.begsample  = Startsample(full_trial)';
                    trialinfo.endsample  = Endsample(full_trial)';
                    trialinfo.offset     = Offset(full_trial)';
                    trialinfo.trialnr    = Trialnr(full_trial)';
                    trialinfo.directory  = Directory(full_trial, :);                    
                    trialinfo.idir       = idir*ones(size(Startsample(full_trial)'));
                    trialinfo.hyplabel   = hyplabels_trl(full_trial)';
                    trialinfo.starttime  = Starttime(full_trial)';
                    trialinfo.endtime    = Endtime(full_trial)';
                end
                
                for ioverlap = 1 : size(cfg.LFP.overlap, 2)
                    othermarkername = cfg.LFP.overlap{ioverlap};
                    trialinfo.([othermarkername '_sec']) = overlap_sec(full_trial, ioverlap);
                    trialinfo.([othermarkername '_cnt']) = overlap_cnt(full_trial, ioverlap);
                end
                
                if isNeuralynx
                    % same label over files
                    filedat{ifile}.label{1} = cfg.LFP.channel{ifile};
                end
                
                % flag for averaging
                hasmarker(idir) = true;
                
            end % ifile (/electrodes)
            
            % concatinate channels
            if exist('filedat','var')
                cfgtemp                             = [];
                cfgtemp.keepsampleinfo              = 'no';
                dirdat{idir}                        = ft_appenddata(cfgtemp, filedat{:});
                clear filedat*
            else % if there are no events at all in this directory
                dirdat{idir} = [];
            end
            
            % additional trialinfo
            dirdat{idir}.trialinfo = trialinfo;
            
        end % idir
        
        if ~exist('dirdat', 'var') % in case there is no marker in the data
            LFP{ipart}.(markername) = [];
            fprintf('%s part %d : No data with marker ''%s''\n', cfg.prefix(1:end-1), ipart, markername);
            continue
        end
        
        % concatinate data of different datasets (over trials)
        LFP{ipart}.(markername)         = ft_appenddata([], dirdat{hasmarker});
        LFP{ipart}.(markername).fsample = dirdat{find(hasmarker, 1)}.fsample;
        clear dirdat*
        
        %remove cfg to save space on disk, if required
        if ~istrue(keepcfg)
            LFP{ipart}.(markername) = rmfield(LFP{ipart}.(markername),'cfg');
        end
        
    end % ipart
    
    if write
        fprintf('Saving LFP data for %s\n', markername);
        saveMarker_LFP(LFP, markername, fname_out)
    end
    
end % markername