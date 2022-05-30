function [t_trial, t_muse, idir] = find_muse_marker_in_trial(data, MuseStruct, itrial, markername, time_from_start, time_from_end, ipart)
% find the first marker after the start of the trial. If no marker is
% found before the end of the trial, return empty values.
% 
% Required :
% - data.trialinfo.idir
% - data.trialinfo.starttime
% - data.trialinfo.endtime
% - data.trialinfo.offset
% - data.fsample or data.hdr.Fs
% - time_from_start and time_from_end : time (in seconds) from the start 
%          end the end of the trial during which the marker is searched 
%          (optional, default = 0)
% - ipart : (optional, default = 1)
% 

%set default values
if nargin < 7 
    ipart = 1;
end
if nargin < 6
    time_from_end = 0;
end
if nargin < 5
    time_from_start = 0;
end

t_trial = [];
t_muse  = [];

if isfield(data, 'hdr')
    fsample = data.hdr.Fs;
elseif isfield(data, 'fsample')
    fsample = data.fsample;
else
    error('There should be a field data.hdr or data.fsample');
end

idir        = data.trialinfo.idir(itrial);

while idir <= size(MuseStruct{ipart}, 2)
    
    if ~isfield(MuseStruct{ipart}{idir}.markers, markername)
        %in case the trial is overlaping 2 directories
        idir = idir +1;
        continue
    end
    
    if ~isfield(MuseStruct{ipart}{idir}.markers.(markername), 'clock')
        %in case the trial is overlaping 2 directories
        idir = idir +1;
        continue
    end
    
    starttrial  = data.trialinfo.starttime(itrial) + seconds(time_from_start);
    idx         = find(MuseStruct{ipart}{idir}.markers.(markername).clock >= starttrial, 1, 'first');
    t_muse      = MuseStruct{ipart}{idir}.markers.(markername).clock(idx);
    
    if t_muse > data.trialinfo.endtime(itrial) + seconds(time_from_end) %remove event if it is not in the corresponding window
        t_trial = [];
        t_muse  = [];
        break
    end
    
    t_trial     = t_muse - starttrial + seconds(data.trialinfo.offset(itrial) / fsample);
    
    %convert form datetime to real in seconds
    t_muse      = MuseStruct{ipart}{idir}.markers.(markername).synctime(idx);
    t_trial     = seconds(t_trial);
    break
end
