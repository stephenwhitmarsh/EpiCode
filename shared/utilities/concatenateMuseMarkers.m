function MuseStruct_concat = concatenateMuseMarkers(cfg, MuseStruct, force)
% 
% Use as : 
%       MuseStruct_concat = concatenateMuseMarkers(cfg, MuseStruct, force)
% 
% concatenateMuseMarkers concatenate Muse markers over all the directories.

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

fname_out = fullfile(cfg.datasavedir, sprintf('%sMuseStruct_concatenated.mat', cfg.prefix));

if exist(fname_out,'file') && force == false
    fprintf('Reading %s\n',fname_out);
    load(fname_out,'MuseStruct_concat');
    return
end

[isNeuralynx, isMicromed, isBrainvision] = get_data_format(cfg);

for ipart = 1:size(MuseStruct,2)
    
    MuseStruct_concat{ipart}.starttime  = MuseStruct{ipart}{1}.starttime;
    MuseStruct_concat{ipart}.endtime    = MuseStruct{ipart}{end}.endtime;
    MuseStruct_concat{ipart}.directory  = cfg.directorylist{ipart};
    
    %find hdr for each dir
    ft_progress('init','text');
    for idir = 1:size(MuseStruct{ipart},2)
        ft_progress(idir/size(MuseStruct{ipart},2),'reading header for dir %d from %d',idir, size(MuseStruct{ipart},2));
        if isNeuralynx
            if ~isempty(cfg.LFP.channel)
                chan_list = cfg.LFP.channel;
            else
                chan_list = cfg.circus.channel;
            end
            temp  	 = dir(fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, ['*', chan_list{1}, '.ncs']));
            datapath = fullfile(cfg.rawdir, cfg.directorylist{ipart}{idir}, temp.name);
        elseif isMicromed
            datapath = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.TRC']);
        elseif isBrainvision
            datapath = fullfile(cfg.rawdir,[cfg.directorylist{ipart}{idir} '.eeg']);
        end
        hdr{ipart}{idir}     = ft_read_header(datapath);
    end
    ft_progress('close');
    
    %find all markers, even if not present in all dirs
    markerlist = [];
    for idir = 1:size(MuseStruct{ipart},2)
        markerlist = [markerlist, string(fieldnames(MuseStruct{ipart}{idir}.markers)')];
        markerlist = unique(markerlist);
    end
    
    %concatenate each marker
    for markername = markerlist
        
        MuseStruct_concat{ipart}.markers.(markername).clock    = datetime.empty;
        MuseStruct_concat{ipart}.markers.(markername).synctime = [];
        MuseStruct_concat{ipart}.markers.(markername).dir      = [];
        length_previous = 0;
        
        %concatenate clock times
        for idir = 1:size(MuseStruct{ipart},2)
            if isfield(MuseStruct{ipart}{idir}.markers, markername)
                if isfield(MuseStruct{ipart}{idir}.markers.(markername), 'clock')
                    
                    clock_temp = MuseStruct{ipart}{idir}.markers.(markername).clock;
                    
                    for isample = 1:size(clock_temp,2)
                        MuseStruct_concat{ipart}.markers.(markername).clock(isample+length_previous) = clock_temp(isample);
                        MuseStruct_concat{ipart}.markers.(markername).dir(isample+length_previous) = idir;
                    end
                    
                    length_previous = length_previous + size(clock_temp,2);
                    
                else
                    %fprintf('Marker %s exists but is empty in %s\n', imarker, MuseStruct{ipart}{idir}.directory);
                end
                
                
            else
                %fprintf('Marker %s does not exist in %s\n',imarker, MuseStruct{ipart}{idir}.directory);
            end
            
        end
        
        % concatenate synctime
        dirOnset = 0;
        for idir = 1:size(MuseStruct{ipart},2)
            if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
                dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
                continue
            end
            if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'synctime')
                dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
                continue
            end
            
            synctime_temp = MuseStruct{ipart}{idir}.markers.(markername).synctime;
            
            for isample = 1:size(synctime_temp,2)
                MuseStruct_concat{ipart}.markers.(markername).synctime(end+1) = synctime_temp(isample)+dirOnset/hdr{ipart}{idir}.Fs;
            end
            
            dirOnset    = dirOnset + hdr{ipart}{idir}.nSamples;
            
        end
        
        ft_info('%d times found for %s\n',size(MuseStruct_concat{ipart}.markers.(markername).clock,2),markername);
    end
end

% save output
save(fname_out,'MuseStruct_concat');

end