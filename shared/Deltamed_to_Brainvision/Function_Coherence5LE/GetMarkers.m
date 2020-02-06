function markers = GetMarkers( z, beginpos, endpos)
% GetMarkers(z, begin, end)
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

markersnumber = GetMarkersNumber(z);
if markersnumber<0
    error('Could not get Markersnumber');
end


sizeofMarker = 23; % size in ints: 80 char + 3 int
sizeofMarkerinBytes = sizeofMarker*4; % in bytes
tmarkers = libpointer('int32Ptr',1:sizeofMarker*markersnumber*2); % don't ask why *2 ... needed with linux 

% we have different function names in LinuxWrapper and WindowsDLL
os = getenv('OS');
if (strfind(os,'Windows')>0)
    functionname = 'Eeg3_GetMarkers2';
else
    functionname = 'Eeg3_GetMarkers';
end

[res,t] = calllib(z.libname,functionname,beginpos, endpos,tmarkers);

if (strfind(os,'Windows')>0) % not so platform independent
    t = tmarkers.Value;
end    

markers(1:markersnumber)=struct('pos',0,'duration',0,'evttype',0,'text','');
for i=0:markersnumber-1    
    markers(i+1) = setfield(markers(i+1),'pos',t(i*sizeofMarker+1));
    markers(i+1) = setfield(markers(i+1),'duration',t(i*sizeofMarker+2));
    markers(i+1) = setfield(markers(i+1),'evttype',t(i*sizeofMarker+3));
    textfield = t(i*sizeofMarker+4:i*sizeofMarker+13);
    markers(i+1) = setfield(markers(i+1),'text',char(typecast(textfield,'uint8')));
    
end


