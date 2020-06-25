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
% cfg.LFP.name                = {'Hspike'};
% cfg.LFP.channel             = {'_Ha2g_1','_Ha2g_2','_Ha2g_3','_Ha2g_4'}
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

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'LFP.mat']);
write     = ft_getopt(cfg.LFP, 'write', false);

if write
    fprintf('*** Will write LFP to: %s ***\n',fname_out);
else
    fprintf('*** Will NOT write LFP to disk ***\n');
end

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed LFP data ***\n');
    fprintf('************************************\n\n');
    load(fname_out,'LFP');
    return
else
    fprintf('********************************\n');
    fprintf('*** (re-) computing LFP data ***\n');
    fprintf('********************************\n\n');
end

    % get file format
    [isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

    % initialize LFP, to return empty array in case there is no LFP to load
    LFP = {};

    % loop over parts within subject
    for ipart = 1:length(MuseStruct)

        % loop over markers
        for imarker = 1 : size(cfg.LFP.name,2)

            fprintf('For marker %s\n',cfg.LFP.name{imarker});
            hasmarker = false(length(MuseStruct{ipart}),1);

            % loop over directories
            for idir = 1:length(MuseStruct{ipart})

                if ~isfield(MuseStruct{ipart}{idir},'markers')
                    continue
                end
                if ~isfield(MuseStruct{ipart}{idir}.markers,(cfg.muse.startend{imarker,1}))
                    continue
                end
                if ~isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}),'synctime')
                    continue
                end
                if isempty(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime)
                    continue
                end

                if isNeuralynx
                    nfile = size(cfg.LFP.channel,2); % one file per channel
                elseif isMicromed
                    nfile = 1; % only one file with all electrodes
                    fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
                elseif isBrainvision
                    nfile = 1; % only one file with all electrodes
                    fname = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
                end

                % loop over files
                for ifile = 1 : nfile

                    %load data
                    if isNeuralynx
                        temp                = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*',cfg.LFP.channel{ifile},'.ncs']));
                        fname{1}            = fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},temp.name);
                        dat                 = ft_read_neuralynx_interp(fname);
                    else
                        cfgtemp.dataset     = fname;
                        cfgtemp.channel     = cfg.labels.macro';
                        dat                 = ft_preprocessing(cfgtemp);
                    end

                    % filtering
                    cfgtemp                 = [];
                    cfgtemp.lpfilter        = ft_getopt(cfg.LFP, 'lpfilter', 'no');
                    cfgtemp.hpfilter        = ft_getopt(cfg.LFP, 'hpfilter', 'no');
                    cfgtemp.bpfilter        = ft_getopt(cfg.LFP, 'bpfilter', 'no');
                    cfgtemp.bsfilter        = ft_getopt(cfg.LFP, 'bsfilter', 'no');
                    cfgtemp.lpfreq          = ft_getopt(cfg.LFP, 'lpfreq', []);
                    cfgtemp.hpfreq          = ft_getopt(cfg.LFP, 'hpfreq', []);
                    cfgtemp.bpfreq          = ft_getopt(cfg.LFP, 'bpfreq', []);
                    cfgtemp.bsfreq          = ft_getopt(cfg.LFP, 'bsfreq', []);
                    dat                     = ft_preprocessing(cfgtemp,dat);

                    % append EMG data (if any)
                    if isMicromed || isBrainvision %not adapted for nlx data (1 file per electrode) for now

                        cfg.EMG = ft_getopt(cfg, 'EMG', []);

                        cfgtemp                   = [];
                        cfgtemp.channel           = {cfg.LFP.emg{imarker}, ft_getopt(cfg.EMG, 'refchannel', [])}; % load the emg associated with eeg marker, and the ref if any
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
                        data_EMG                = ft_preprocessing(cfgtemp,data_EMG);

                        % append EMG to EEG data
                        cfgtemp                 = [];
                        cfgtemp.keepsampleinfo  = 'yes';
                        dat                     = ft_appenddata(cfgtemp, dat, data_EMG);
                    end


                    % downsample data and correct baseline
                    if isfield(cfg.LFP,'resamplefs')

                        cfgtemp                         = [];
                        cfgtemp.resamplefs              = ft_getopt(cfg.LFP, 'resamplefs', []);
                        if strcmp(cfg.LFP.baseline, 'no')
                            cfgtemp.demean              = 'no';
                        else
                            cfgtemp.demean              = 'yes';
                            cfgtemp.baselinewindow      = cfg.LFP.baselinewindow{imarker};
                        end
                        dat                             = ft_resampledata(cfgtemp,dat);
                    end

                    fsample = dat.fsample; %store this info for output LFP

                    % create trial segmentation common to resampled
                    % data. Neuralynx : same markers for all files
                    % of one dir.
                    Startsample             = [];
                    Endsample               = [];
                    Stage                   = [];
                    Offset                  = [];
                    trialnr                 = [];
                    for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime,2)

                        ss = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,1}).synctime(ievent) *  dat.fsample);
                        if strcmp(cfg.muse.startend{imarker,1},cfg.muse.startend{imarker,2})
                            es = ss;
                        else
                            idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime * dat.fsample) >= ss,1,'first');
                            es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{imarker,2}).synctime(idx) * dat.fsample);
                        end

                        if isempty(es)
                            continue
                        end

                        Startsample(ievent) = ss + cfg.epoch.toi{imarker}(1) * dat.fsample - cfg.epoch.pad{imarker} * dat.fsample;
                        Endsample(ievent)   = es + cfg.epoch.toi{imarker}(2) * dat.fsample + cfg.epoch.pad{imarker} * dat.fsample;
                        Offset(ievent)      = (cfg.epoch.toi{imarker}(1) - cfg.epoch.pad{imarker}) * dat.fsample;
                        trialnr(ievent)     = ievent;
                        Stage(ievent)       = -1;

                        % find overlap with hypnogram markers
                        for hyplabel = {'PHASE_1','PHASE_2','PHASE_3','REM','AWAKE','NO_SCORE'}
                            if ~isfield(MuseStruct{ipart}{idir}.markers,[cell2mat(hyplabel),'__START__'])
                                continue
                            end
                            for i = 1 : size(MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime,2)
                                y1 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__START__']).synctime(i) * dat.fsample;
                                y2 = MuseStruct{ipart}{idir}.markers.([cell2mat(hyplabel),'__END__']).synctime(i) * dat.fsample;
                                if (y1 < ss) && (ss < y2)
                                    fprintf('Found "%s" overlapping with "%s" : adding to trialinfo: ',cfg.name{imarker},cell2mat(hyplabel));
                                    switch cell2mat(hyplabel)
                                        case 'PHASE_1'
                                            fprintf('%d\n',1);
                                            Stage(ievent) = 1;
                                        case 'PHASE_2'
                                            fprintf('%d\n',2);
                                            Stage(ievent) = 2;
                                        case 'PHASE_3'
                                            fprintf('%d\n',3);
                                            Stage(ievent) = 3;
                                        case 'REM'
                                            fprintf('%d\n',4);
                                            Stage(ievent) = 4;
                                        case 'AWAKE'
                                            fprintf('%d\n',0);
                                            Stage(ievent) = 0;
                                        case 'NO_SCORE'
                                            fprintf('%d\n',0);
                                            Stage(ievent) = 0;
                                        otherwise
                                            error('Unexpected label name in Hypnogram\n');
                                    end
                                end
                            end
                        end
                    end

                    % create Fieldtrip trl
                    cfgtemp                         = [];
                    cfgtemp.trl                     = round([Startsample; Endsample; Offset]');
                    cfgtemp.trl(:,4)                = trialnr;
                    cfgtemp.trl(:,6)                = idir;
                    cfgtemp.trl(:,7)                = Stage;
                    cfgtemp.trl(:,8)                = Startsample;
                    cfgtemp.trl(:,9)                = Endsample;
                    cfgtemp.trl(:,10)               = Offset;
                    cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < length(dat.trial{1}),:); % don't read before BOF or after EOF
                    filedat{ifile}                  = ft_redefinetrial(cfgtemp,dat);
                    clear dat

                    if isNeuralynx
                        % same label over files
                        filedat{ifile}.label{1}     = cfg.LFP.channel{ifile};
                    end

                    % flag for averaging
                    hasmarker(idir)                 = true;
                end

                % concatinate channels
                cfgtemp                             = [];
                cfgtemp.keepsampleinfo              = 'no';
                dirdat{idir}                        = ft_appenddata(cfgtemp,filedat{:});
                clear filedat*

            end % idir

            if exist('dirdat','var') % in case there is no marker in the data

                % concatinate data of different datasets (over trials)
                LFP{ipart}{imarker} = ft_appenddata([],dirdat{find(hasmarker)});
                LFP{ipart}{imarker}.fsample = fsample;
                clear dirdat*
            else
                LFP{ipart}{imarker} = [];
                fprintf('%s part %d : No data with marker ''%s''\n',cfg.prefix(1:end-1), ipart, cfg.LFP.name{imarker});
            end

        end % imarker

    end % ipart


if write
    save(fname_out,'LFP','-v7.3');
end
