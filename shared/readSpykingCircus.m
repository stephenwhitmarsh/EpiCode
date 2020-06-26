function [SpikeTrials] = readSpykingCircus(cfg,MuseStruct,SpikeRaw,force,varargin)


% READSPYKINGCIRCUS reads Spyking-Circus results: spike-times and templates.
%
% use as
% [SpikeRaw, SpikeTrials] = readSpykingCircus(cfg,MuseStruct,force,varargin)
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


fname = fullfile(cfg.datasavedir,[cfg.prefix,'SpikeTrials_Muse.mat']);
if exist(fname,'file') && force == false
    load(fname,'SpikeTrials');
    return;
else

    for ipart = parts_to_read

        % find spiking-circus output path, which is based on the name of the
        % first datafile
%         temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.result',cfg.circus.postfix,'.hdf5']));
%         if isempty(temp)
%             fprintf('Could not find Spyking-Circus results: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.result',cfg.circus.postfix,'.hdf5']));
%             return
%         else
%             fname_spikes = fullfile(temp.folder,temp.name);
%         end
        
%         temp = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.templates',cfg.circus.postfix,'.hdf5']));
%         if isempty(temp)
%             fprintf('Could not find Spyking-Circus templates: %s\n',fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],'SpykingCircus',[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.templates',cfg.circus.postfix,'.hdf5']));
%             return
%         end
%         fname_templates = fullfile(temp.folder,temp.name);
        
        temp        = dir(fullfile(cfg.datasavedir,cfg.prefix(1:end-1),['p',num2str(ipart)],[cfg.prefix,'p',num2str(ipart),'-multifile-',cfg.circus.channel{1}(1:end-2),'*.ncs']));
        hdr_fname   = fullfile(temp(1).folder,temp(1).name);
        hdr         = ft_read_header(hdr_fname); % take the first file to extract the header of the data
          

            % create trials
            clear Trials
            for ilabel = 1 : size(cfg.name,2)

                % clock time of each event
                clocktimes = [];
                for ifile = 1 : size(MuseStruct{ipart},2)
                    if isfield(MuseStruct{ipart}{ifile}.markers,cfg.muse.startend{ilabel})
                        if isfield(MuseStruct{ipart}{ifile}.markers.(cfg.muse.startend{ilabel}),'clock')
                            clocktimes = [clocktimes, MuseStruct{ipart}{ifile}.markers.(cfg.muse.startend{ilabel}).clock];
                        end
                    end
                end

                % create Fieldtrip trl based on concatinated files by adding nr of
                % samples of each file
                Startsample = [];
                Endsample   = [];
                Offset      = [];
                Trialnr     = [];
                Trialdir    = [];
                Filenr      = [];
                FileOffset  = [];
                
                dirOnset(1) = 0;

                trialcount = 1;
                for idir = 1 : size(MuseStruct{ipart},2)
                    %compute dir onset of the next dir based of the length
                    %of this dir. Do it before searching for events to avoid
                    %error if one dir has no event
                     temp        = dir(fullfile(cfg.rawdir,cfg.directorylist{ipart}{idir},['*', cfg.circus.channel{1}(1:end-2),'*.ncs']));
                     hdrtemp     = ft_read_header(fullfile(temp(1).folder,temp(1).name));
                     dirOnset(idir+1)    = dirOnset(idir) + hdrtemp.nSamples; % assuming all channels have same sampleinfo
                    
                    if isfield(MuseStruct{ipart}{idir}.markers,cfg.muse.startend{ilabel})
                        if isfield(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}),'synctime')
                            if ~isempty(size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}).synctime,2))
                                trialcount_dir = 1;
                                
                                for ievent = 1 : size(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel}).synctime,2)
                                    % FIXME
                                    % has to be fixezd for unequalnumber of
                                    % start-end, i.e. for Paul's project's
                                    % Interictal period
                                    %end of trial : take the following marker, no need to have the same index has begining marker
                                    ss  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,1}).synctime(ievent) * hdr.Fs);
                                    idx = find(round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime * hdr.Fs) >= ss,1,'first');
                                    es  = round(MuseStruct{ipart}{idir}.markers.(cfg.muse.startend{ilabel,2}).synctime(idx) * hdr.Fs);
                                    
                                    if ~isempty(es)
                                        Startsample  = [Startsample; ss + cfg.epoch.toi{ilabel}(1) * hdr.Fs + dirOnset(idir)];
                                        Endsample    = [Endsample;   es + cfg.epoch.toi{ilabel}(2) * hdr.Fs + dirOnset(idir)];
                                        Offset       = [Offset; cfg.epoch.toi{ilabel}(1) * hdr.Fs];
                                        Trialdir     = [Trialdir; trialcount_dir];
                                        Trialnr      = [Trialnr; trialcount];
                                        Filenr       = [Filenr; idir];
                                        FileOffset   = [FileOffset; dirOnset(idir)];
                                        trialcount   = trialcount + 1;
                                        trialcount_dir = trialcount_dir+1;
                                    end
                                end
                            else
                                fprintf('No events starting with %s found in filenr %d\n',cfg.muse.startend{ilabel},idir);
                            end
                        end
                    end

                end

                cfgtemp                         = [];
                cfgtemp.trl                     = [Startsample, Endsample, Offset];
                cfgtemp.trl(:,4)                = ones(size(cfgtemp.trl,1),1) * idir;
                cfgtemp.trl(:,5)                = Trialdir;                % trialnr. to try to find trials that are missing afterwards
                cfgtemp.trl(:,6)                = Startsample;                          % startsample
                cfgtemp.trl(:,7)                = Endsample;                            % endsample
                cfgtemp.trl(:,8)                = Offset;                               % offset
                cfgtemp.trl(:,9)                = Endsample-Startsample+1;              % duration in samples
                cfgtemp.trl(:,10)               = Trialnr;                              % trialnr. to try to find trials that are missing afterwards
                cfgtemp.trl(:,11)               = Filenr;                               % trialnr. to try to find trials that are missing afterwards
                cfgtemp.trl(:,12)               = FileOffset;                           % trialnr. to try to find trials that are missing afterwards

                cfgtemp.trl                     = cfgtemp.trl(Startsample > 0 & Endsample < hdr.nSamples,:); % so not to read before BOF or after EOFs

                % create spiketrials timelocked to events
                cfgtemp.trlunit                         = 'samples';
                cfgtemp.hdr                             = hdr;
                SpikeTrials{ipart}{ilabel}              = ft_spike_maketrials(cfgtemp,SpikeRaw{ipart});
                SpikeTrials{ipart}{ilabel}.clocktimes   = clocktimes;
                SpikeTrials{ipart}{ilabel}.hdr          = hdr;
                SpikeTrials{ipart}{ilabel}.label        = cfg.name{ilabel};

                % commented out on 9-8-2019 after looking with Zoe
                %             SpikeRaw.time{ilabel}           = SpikeRaw.samples{ilabel} / hdr.Fs;
                %             SpikeRaw.trial{ilabel}          = ones(size(SpikeRaw.samples{ilabel}));

            end % patterns
            SpikeRaw{ipart}.trialtime = [0 hdr.nSamples / hdr.Fs];

        end % if exists spike file.

    end
    
    save(fname,'SpikeTrials');
        
end % save / force

    %     save(fname,'SpikeRaw');

end % save / force
