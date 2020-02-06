function number = GetMarkersNumber( z, beginpos, endpos )
% int = getMarkersNumber(z)
% returns number of Markers
% if <0 then error occurred
% INPUT  : z 
%          begin (optional): 1 by default
%          end   (optional): lastsample by default

if ~exist('beginpos','var')
    beginpos = 1;
end

if ~exist('endpos','var')
    endpos= z.FileInfo.duration * z.FileInfo.frequency;
end



z = ensureload(z);

% Get Number
res =calllib(z.libname,'Eeg3_GetMarkersNumber', beginpos, endpos);

number = res;